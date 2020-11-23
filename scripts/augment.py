import os
import sys
import click
import data_preprocessor as dp
from typing import TextIO
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

DOCKER = os.getenv("RUNNING_IN_DOCKER")
TOKENIZER = dp.SmilesTokenizer()


@click.command()
@click.argument("input", type=click.File("r"), required=False)
@click.argument("output", nargs=1, required=False)
@click.option("--tokenize/--no-tokenize", default=True)
@click.option("--random_type", 'r', default='unrestricted',
              help='the randomization type for the SMILES. Available: unrestricted, restricted, rotated, molecules')
@click.option("--permutations", '-n', default=1, help='number of generated permutations per SMILES')
@click.option("--fragment_bond", default="~", help='fragment bond token in the SMILES')
def cli(input: TextIO, output: str, tokenize: bool, random_type: str, permutations: int, fragment_bond: str) -> None:
    """The entry point for this cli script.

    Args:
        input (TextIO):  The input file (one SMILES per line).
        output (TextIO): The output file (one SMILES per line).
        tokenize (bool): Whether to tokenize the augmented SMILES.
        random_type (str): Type of randomization for the SMILES. Available: unrestricted, restricted, rotated, molecules
        permutations (int): Number of permutations to retain for each of the SMILES
        fragment_bond (str): Token for fragment bond in the SMILES
    """

    # If not running in docker, require input and output file.
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

    # Create a instance of the Augmenter.
    ag = dp.Augmenter.read_csv(input.name, fragment_bond)

    # Perform augmentation
    augm = ag.augment(random_type, permutations, tokenize)

    # Exporting augmented samples
    # augm[f"{random_type}"].to_csv(f"{random_type}.csv")


if __name__ == "__main__":
    cli()
