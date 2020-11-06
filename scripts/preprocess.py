import sys
import os
import re
import sys
import pathlib
import click
import numpy as np
import pandas as pd
import data_preprocessor as dp
from collections import Counter
from typing import TextIO, Iterable, Tuple, List
from itertools import chain
from functools import partial
from multiprocessing import Pool, cpu_count
from rdkit.Chem import AllChem as rdk
from rdkit import RDLogger
from crc64iso.crc64iso import crc64
from tabulate import tabulate

RDLogger.DisableLog("rdApp.*")

DOCKER = os.getenv("RUNNING_IN_DOCKER")
TOKENIZER = dp.SmilesTokenizer()


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
            input[i] = TOKENIZER.tokenize(input[i])

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
        for reaction, processed_reaction, invalid_reasons in input:
            f.write(
                reaction
                + ","
                + processed_reaction
                + ","
                + '"'
                + ",".join(invalid_reasons)
                + '"'
                "\n"
            )


def summarize_invalids(invalid: Iterable[Tuple[str, str, List[str]]]) -> None:
    """Prints a table summarazing the results of the process.

    Args:
        invalid (Iterable[Tuple[str, List[str]]]): The list of invalid reactions with reasons.
    """
    reasons = []
    for inv in invalid:
        reasons.extend(inv[2])

    print(
        f"\033[93m- {str(len(invalid))} reactions were removed for the following reasons:"
    )
    headers = ["Reason", "Number of Reactions"]
    print(
        tabulate(list(Counter(reasons).items()), headers, tablefmt="fancy_grid")
        + "\033[0m"
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

    mrf = dp.MixedReactionFilter()
    invalid_reactions = []
    processed_reactions = []

    for i in range(len(lines) - 1, -1, -1):
        invalid_reasons = []
        reaction = dp.Reaction(lines[i], remove_duplicates=True)

        # If reactions contains None values for molecules, skip it
        if reaction.has_none():
            invalid_reasons.append("rdkit_molfromsmiles_failed")

        # Move agents to reactants
        reaction.reactants.extend(reaction.agents)
        reaction.agents = []

        reaction.filter(([], [], reaction.find("[S+,s+](*)(*)*")[2]))

        # Remove single atoms
        reaction.reactants = [
            m for m in reaction.reactants if m and m.GetNumAtoms() > 1
        ]
        reaction.agents = [m for m in reaction.agents if m and m.GetNumAtoms() > 1]
        reaction.products = [m for m in reaction.products if m and m.GetNumAtoms() > 1]

        # Remove products that are also reactants
        reaction.remove_precursors_from_products()
        reaction.sort()

        _, mrf_reasons = mrf.validate_reasons(reaction)
        invalid_reasons.extend(mrf_reasons)

        if len(invalid_reasons) > 0:
            invalid_reactions.append((lines[i], str(reaction), invalid_reasons))
        else:
            processed_reactions.append(str(reaction))

    return (processed_reactions, invalid_reactions)


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
    export_invalids(invalid, "invalid", output)

    sdd = dp.StableDataSplitter()

    # # Split and save
    df = pd.DataFrame({"rxn": valid})
    train, validation, test = sdd.split(df, index_column="rxn")
    df_train = df[train]
    df_validation = df[validation]
    df_test = df[test]

    print(f"- Training set size: {len(df_train)}.")
    print(f"- Validation set size: {len(df_validation)}.")
    print(f"- Test set size: {len(df_test)}.")
    print("")

    to_file(df_train["rxn"].tolist(), "train", output, tokenize=tokenize, split=split)
    to_file(
        df_validation["rxn"].tolist(), "valid", output, tokenize=tokenize, split=split
    )
    to_file(df_test["rxn"].tolist(), "test", output, tokenize=tokenize, split=split)


if __name__ == "__main__":
    cli()