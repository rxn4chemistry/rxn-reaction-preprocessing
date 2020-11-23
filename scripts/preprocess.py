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

    # This is for the special SMILES extension where agents are separated by pipe.
    def clean_func(rxn: str) -> str:
        return re.sub(
            r"(?<=\[)([0-9]+)(?=[A-Za-z])",  # Remove isotopes
            "",
            rxn.strip()
            .replace("(|", "(")  # Hack for bad input data
            .replace(")|", ")")  # Hack for bad input data
            .replace("~", ".")
            .replace("|", ".")
            .replace("..", ".")
            .replace(">.", ">")
            .replace(".>", ">"),
        )

    # This is the function that is applied to each reaction.
    def apply_func(reaction: dp.Reaction) -> dp.Reaction:
        # Move agents to reactants
        reaction.reactants.extend(reaction.agents)
        reaction.agents = []
        reaction.filter(
            ([], [], reaction.find_in("[S+,s+](*)(*)*", dp.ReactionPart.products))
        )

        # Remove products that are also reactants
        reaction.remove_precursors_from_products()

        return reaction.sort()

    # Create a instance of the mixed reaciton filter with default values.
    # Make arguments for all properties in script
    mrf = dp.MixedReactionFilter()

    pp = dp.Preprocessor.read_csv(input.name, "rxn")

    # Remove duplicate reactions (useful for large dataset, this step is repeated later)
    pp.remove_duplicates()

    # In a first step, let's clean the data using the cleaning function
    # defined above
    pp.df.rxn = pp.df.rxn.apply(clean_func)

    # Apply the function above to all reactions, the remove_duplicate_molecules argument
    # is set to true to remove duplicate molecules within each reaction part
    pp.apply(apply_func, remove_duplicate_molecules=True)

    # Remove duplicate reactions
    pp.remove_duplicates()

    # Apply the mixed reaction filter instance defined above, enable verbose mode
    pp.filter(mrf, True)

    # Print the detailed stats
    pp.print_stats()

    # Drop the invalid ractions
    pp.remove_invalids()

    # After dropping invalid columns, display stats again (as an example)
    pp.print_stats()

    # Tokenize the reactions
    tokenizer = dp.SmilesTokenizer()
    pp.df.rxn = pp.df.rxn.apply(tokenizer.tokenize)

    # Split into train, validation, and test sets, but do not export yet
    train, validation, test = dp.StableDataSplitter.split(pp.df, "rxn")

    # Example of exporting one of the sets
    # pp.df.rxn[train].to_csv("training.csv")


if __name__ == "__main__":
    cli()