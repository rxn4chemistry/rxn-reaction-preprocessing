import os
import sys
from typing import TextIO

import click
from rdkit import RDLogger

import rxn_reaction_preprocessing as rrp

RDLogger.DisableLog('rdApp.*')

DOCKER = os.getenv('RUNNING_IN_DOCKER')
TOKENIZER = rrp.SmilesTokenizer()


@click.command()
@click.argument('input_file', type=click.File('r'), required=False)
@click.argument('output_directory', nargs=1, required=False)
@click.option('--tokenize/--no-tokenize', default=True)
@click.option(
    '--random_type',
    'r',
    default='unrestricted',
    help=
    'the randomization type for the SMILES. Available: unrestricted, restricted, rotated, molecules'
)
@click.option(
    '--permutations', '-n', default=1, help='number of generated permutations per SMILES'
)
@click.option('--fragment_bond', default='~', help='fragment bond token in the SMILES')
def cli(
    input_file: TextIO, output_directory: str, tokenize: bool, random_type: str, permutations: int,
    fragment_bond: str
) -> None:
    """The entry point for this cli script.

    Args:
        input (TextIO):  The input file (one SMILES per line).
        output (str): The output file (one SMILES per line).
    """

    # If not running in docker, require input and output file.
    # In docker these will be supplied by the volumes
    if not DOCKER:
        if not input_file:
            SystemExit('Please specify an input file.')
        if not output_directory:
            SystemExit('Please specify an output directory.')
        output = output_directory.rstrip('/')
    else:
        input_file = open('/data/input.txt')
        output_directory = '/data/output'

    # Create a instance of the Augmenter.
    ag = rrp.Augmenter.read_csv(input_file.name, fragment_bond)

    # Perform augmentation
    augm = ag.augment(random_type, permutations, tokenize)

    # Exporting augmented samples
    # augm[f"{random_type}"].to_csv(f"{random_type}.csv")


if __name__ == '__main__':
    cli()
