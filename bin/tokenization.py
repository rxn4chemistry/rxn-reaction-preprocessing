import os
import sys
import click
import pandas as pd
import data_preprocessor as dp
from typing import TextIO

DOCKER = os.getenv("RUNNING_IN_DOCKER")


@click.command()
@click.argument("input", nargs=1, required=False)
def cli(input: str) -> None:
    """The entry point for this cli script.

    Args:
        input (str):  The input folder.
    """

    # If not running in docker, require intput and output file.
    # In docker these will be supplied by the volumes
    if not DOCKER:
        if not input:
            print("Please specify an input file.")
            sys.exit(1)
    else:
        input = "/data/input"

    # Tokenize the reactions
    tokenizer = dp.SmilesTokenizer()

    for filename in os.listdir(input):
        filename = os.path.join(input, filename)
        if filename[-4:] != ".csv":
            continue

        df = pd.read_csv(filename)

        if "rxn" not in df.columns:
            continue

        df["rxn_precursors"] = df.rxn.str.split(">>").str[0]
        df["rxn_products"] = df.rxn.str.split(">>").str[1]

        df.rxn_precursors = df.rxn_precursors.apply(tokenizer.tokenize)
        df.rxn_products = df.rxn_products.apply(tokenizer.tokenize)

        df[["rxn_precursors"]].to_csv(
            f"{filename}.precursors_tokens", header=False, index=False
        )
        df[["rxn_products"]].to_csv(
            f"{filename}.products_tokens", header=False, index=False
        )


if __name__ == "__main__":
    cli()