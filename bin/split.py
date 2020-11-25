import os
import sys
import click
import pandas as pd
import data_preprocessor as dp
from typing import TextIO

DOCKER = os.getenv("RUNNING_IN_DOCKER")


@click.command()
@click.argument("input", type=click.File("r"), required=False)
@click.argument("output", nargs=1, required=False)
@click.option("--split-ratio", default=0.05)
def cli(input: TextIO, output: str, split_ratio: float) -> None:
    """The entry point for this cli script.

    Args:
        input (TextIO):  The input file (one reaction SMARTS per line).
        output (TextIO): The output file (one reaction SMARTS per line).
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

    df = pd.read_csv(input.name)

    # Split into train, validation, and test sets, but do not export yet
    train, validation, test = dp.StableDataSplitter.split(
        df, "rxn", split_ratio=split_ratio
    )

    # Example of exporting one of the sets
    df.rxn[train].to_csv(os.path.join(output, "processed.training.csv"))
    df.rxn[validation].to_csv(os.path.join(output, "processed.validation.csv"))
    df.rxn[test].to_csv(os.path.join(output, "processed.test.csv"))


if __name__ == "__main__":
    cli()