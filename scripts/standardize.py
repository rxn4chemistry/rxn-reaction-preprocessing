import os
import sys
import click
import data_preprocessor as dp
from typing import TextIO
from rdkit import RDLogger
from data_preprocessor.utils import data_directory

RDLogger.DisableLog("rdApp.*")

DOCKER = os.getenv("RUNNING_IN_DOCKER")
TOKENIZER = dp.SmilesTokenizer()


@click.command()
@click.argument("input", type=click.File("r"), required=False)
@click.argument("output", nargs=1, required=False)
@click.option("--valid_column", default="rxn", help='column name for the SMILES of the reactions to process')
@click.option("--fragment_bond", default=".", help='fragment bond token in the SMILES of the reactions to process')
def cli(input: TextIO, output: str, valid_column: str, fragment_bond: str) -> None:
    """The entry point for this cli script.

    Args:
        input (TextIO):  The input file (one SMILES per line).
        output (str): The output file (one SMILES per line).
        valid_column (str): The column name for the SMILES of the reactions to process
        fragment_bond (str): The fragment bond token used in the files to be standardized
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

    # Create a instance of the Patterns.
    # for now jsonpath and fragment_bond (the one present in the jsonfile) fixed
    jsonfilepath = f"{str(data_directory())}/standardization-files/pistachio-200302"
    pt = dp.Patterns(jsonfilepath, fragment_bond='~')

    # Create an instance of the Standardizer
    std = dp.Standardizer.read_csv(input.name, pt, valid_column, fragment_bond)
    # Perform standardization
    std.standardize()

    # Exporting standardized samples
    std.df[f"std_{valid_column}"].to_csv(os.path.join(output,f"std_{valid_column}.csv"))


if __name__ == "__main__":
    cli()
