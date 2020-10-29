import sys
import os
import re
import sys
import pathlib
import click
import numpy as np
from collections import Counter
from typing import TextIO, Iterable, Tuple, List
from itertools import chain
from functools import partial
from multiprocessing import Pool, cpu_count
from rdkit.Chem import AllChem as rdk
from rdkit import RDLogger
from tokenization import tokenize_smiles
from reaction_filter import MixedReactionFilter
from crc64iso.crc64iso import crc64
from tabulate import tabulate

RDLogger.DisableLog("rdApp.*")

DOCKER = os.getenv("RUNNING_IN_DOCKER")


def to_file(
    input: Iterable[str],
    suffix: str,
    out_path: str,
    tokenize: bool = True,
    split: bool = True,
) -> None:
    """Save a data set to a file, tokenize and split reactants and products into separate files.

    Args:
        input (Iterable[str]): The input Iterable.
        suffix (str): The suffix to apply to the file name (e.g. "train", "valid", or "test").
        out_path (str): The path to the output directory.
        tokenize (bool, optional): Whether to tokenize the SMILES. Defaults to True.
        split (bool, optional): Whether to split the reactants and products into separate files. Defaults to True.
    """
    if tokenize:
        for i in range(len(input)):
            input[i] = tokenize_smiles(input[i])

    if split:
        with open(
            "{0}/{1}-{2}.txt".format(out_path, "precursors", suffix),
            "w+",
        ) as f_r:
            with open(
                "{0}/{1}-{2}.txt".format(out_path, "product", suffix),
                "w+",
            ) as f_p:
                for reaction in input:
                    reaction_split = reaction.split(">>")
                    f_r.write(reaction_split[0] + "\n")
                    f_p.write(reaction_split[1] + "\n")
    else:
        with open("{0}/{1}.txt".format(out_path, suffix), "w+") as f:
            for reaction in input:
                f.write(reaction + "\n")


def export_invalids(
    input: Iterable[Tuple[str, List[str]]], suffix: str, out_path: str
) -> None:
    """Save the invalid reactions with the reasons for invalidity to a file.

    Args:
        input (Iterable[Tuple[str, List[str]]]): The input Iterable.
        suffix (str): The suffix to apply to the file name (e.g. "train", "valid", or "test").
        out_path (str): The path to the output directory.
    """
    with open("{0}/{1}.txt".format(out_path, suffix), "w+") as f:
        for reaction, invalid_reasons in input:
            f.write(reaction + "," + '"' + ",".join(invalid_reasons) + '"' "\n")


def summarize_invalids(invalid: Iterable[Tuple[str, List[str]]]) -> None:
    """Prints a table summarazing the results of the process.

    Args:
        invalid (Iterable[Tuple[str, List[str]]]): The list of invalid reactions with reasons.
    """
    reasons = []
    for inv in invalid:
        reasons.extend(inv[1])

    print(
        f"\033[93m- {str(len(invalid))} reactions were removed for the following reasons:"
    )
    headers = ["Reason", "Number of Reactions"]
    print(
        tabulate(list(Counter(reasons).items()), headers, tablefmt="fancy_grid")
        + "\033[0m"
    )


def stable_split(
    input: Iterable, split_ratio: float = 0.05
) -> (Iterable, Iterable, Iterable):
    """Generate stable splits of the input iterable into train, validation, and test sets.

    Args:
        input (Iterable): The input array to be split.
        split_ratio (split_ratio): The ratio of between the train set and the two other sets. Defaults to 0.05.

    Returns:
        (Iterable, Iterable, Iterable): The train, validation, and test set in this order.
    """
    intput_dict = {}
    for element in input:
        intput_dict[int(crc64(element), 16)] = element

    test = [v for k, v in intput_dict.items() if k < split_ratio * 2 ** 64]
    valid = [
        v
        for k, v in intput_dict.items()
        if k >= split_ratio * 2 ** 64 and k < split_ratio * 2 * 2 ** 64
    ]
    train = [v for k, v in intput_dict.items() if k >= split_ratio * 2 * 2 ** 64]

    return (train, valid, test)


def filter_sort_molecules(
    input: str,
    sep: str = ".",
    remove_single_atoms: bool = False,
    substructure: str = None,
) -> Iterable:
    """Filters and sorts a side of a reaction where molecules are separated by sep.

    Args:
        input (str): The input SMILES string.
        sep (str, optional): The delimiter used to separate molcules. Defaults to ".".
        remove_single_atoms (bool, optional): Whether to remove single atoms. Defaults to False.
        substructure (str, optional): A SMARTS pattern to filter for substructure. Defaults to None.

    Returns:
        Iterable: The filtered input string.
    """

    pattern = None

    if substructure:
        pattern = rdk.MolFromSmarts(substructure)

    result = []
    for s in input.split("."):
        m = rdk.MolFromSmiles(s)

        if pattern and len(list(m.GetSubstructMatch(pattern))) < 1:
            continue

        if m.GetNumAtoms() < 2 and remove_single_atoms:
            continue

        result.append(rdk.MolToSmiles(m))

    # Deduplicate
    result = list(dict.fromkeys(result))

    return sorted(
        result,
        key=str.lower,
    )


def process(lines: Iterable) -> Iterable:
    """Process the reaction SMARTs. See comments for details.

    Args:
        lines (Iterable): Lines containing reaction SMARTS (one per line).

    Returns:
        [Iterable]: A list of processed SMARTS.
    """

    # The regex removes isotopes
    lines = [
        re.sub(
            r"(?<=\[)([0-9]+)(?=[A-Za-z])",
            "",
            line.strip()
            .replace("(|", "(")  # Hack for bad input data
            .replace(")|", ")")  # Hack for bad input data
            .replace("~", ".")
            .replace("|", ".")
            .replace("..", ".")
            .replace(">.", ">")
            .replace(".>", ">"),
        )
        for line in lines
    ]

    original_lines = lines.copy()

    invalid_lines = []
    lines = [line.split(">") for line in lines]
    # Add the reagents to the reactants, validate and canonicalize all molecules, and
    # order the resulting SMILES alphabetically. Finally concat into reaction SMARTS.
    for i in range(len(lines) - 1, -1, -1):
        lines[i][0] += "." + lines[i].pop(1)
        lines[i][0] = lines[i][0].strip(".")

        # Remove duplicate reactants and products, canonicalize SMILES, and sort
        try:
            lines[i][0] = filter_sort_molecules(lines[i][0])
            lines[i][1] = filter_sort_molecules(
                lines[i][1], remove_single_atoms=True, substructure="[S+,s+](*)(*)*"
            )
        except:
            # Error will be printed by RDKit. Drop the line with the error
            lines.pop(i)
            invalid_lines.append((original_lines[i], ["rdkit_parser_error"]))
            continue

        # Remove products which are also reactants
        lines[i][1] = [m for m in lines[i][1] if m not in lines[i][0]]

        # Drop single reactant reaction or reactions without a product
        if len(lines[i][0]) < 2 or len(lines[i][1]) < 1:
            lines.pop(i)
            invalid_lines.append((original_lines[i], ["single_reactant_or_no_product"]))
            continue

        lines[i][0] = ".".join(lines[i][0])
        lines[i][1] = ".".join(lines[i][1])
        lines[i] = ">>".join(lines[i])

        # Apply the filter by ato
        mrf = MixedReactionFilter(lines[i])
        result = mrf.apply_all_filters()
        if not result[0]:
            lines.pop(i)
            invalid_lines.append((original_lines[i], result[1]))

    return (lines, invalid_lines)


@click.command()
@click.argument("input", type=click.File("r"), required=False)
@click.argument("output", nargs=1, required=False)
@click.option("--split/--no-split", default=True)
@click.option("--tokenize/--no-tokenize", default=True)
def cli(input: TextIO, output: str, split: bool, tokenize: bool) -> None:
    """The entry point for this cli script.

    Args:
        input (TextIO):  The input file (one reaction SMARTS per line).
        output (TextIO): The output file (one reaction SMARTS per line).
        split (bool): Whether to split the reaction into two separate reactants and products files.
        tokenize (bool): Whether to tokenize the files.
    """

    # If not running in docker, require intput and output file.
    # In docker these will be supplied by the volumes
    if not DOCKER:
        if not input:
            print("Please specify an input file.")
            sys.exit(1)
        if not output:
            print("Please specify an output directory.")
            sys.exit(1)
        output = output.rstrip("/")
    else:
        input = open("/data/input.txt", "r")
        output = "/data/output"

    # Going with multiprocessing rather than numpy vectorization due to the latter's
    # poor performance on strings
    result = []
    n_cores = cpu_count() - 1 or 1
    with Pool(n_cores) as pool:
        result = list(
            pool.map(
                partial(process),
                np.array_split(input.readlines(), n_cores),
            )
        )

    valid = []
    invalid = []
    for val, inval in result:
        valid.extend(val)
        invalid.extend(inval)

    valid_len = len(valid)
    invalid_len = len(invalid)

    # Remove duplicates
    valid = list(dict.fromkeys(valid))

    print(f"- {valid_len + invalid_len} total reactions.")
    print(f"\033[92m- {len(valid)} valid reactions.\033[0m")
    print(f"\033[93m- {valid_len - len(valid)} valid duplicates removed.\033[0m")

    summarize_invalids(invalid)

    # Write out the invalid reactions
    export_invalids(invalid, "invalid", output)

    # Split and save
    train, valid, test = stable_split(valid)

    print(f"- Training set size: {len(train)}.")
    print(f"- Validation set size: {len(valid)}.")
    print(f"- Test set size: {len(test)}.")
    print("")

    to_file(train, "train", output, tokenize=tokenize, split=split)
    to_file(valid, "valid", output, tokenize=tokenize, split=split)
    to_file(test, "test", output, tokenize=tokenize, split=split)


if __name__ == "__main__":
    cli()