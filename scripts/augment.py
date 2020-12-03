import os
import sys
from typing import TextIO

import click

import rxn_reaction_preprocessing as rrp


@click.command()
@click.argument('input_file_path', type=click.Path(exists=True), required=True)
@click.argument('outupt_file_path', type=click.Path(), required=True)
@click.option('--tokenize/--no-tokenize', default=True)
@click.option(
    '--random_type',
    'r',
    default='unrestricted',
    type=click.Choice(
        ['unrestricted', 'restricted', 'restricted', 'rotated', 'molecules'], case_sensitive=True
    ),
    help=
    'the randomization type for the SMILES. Available: unrestricted, restricted, rotated, molecules'
)
@click.option(
    '--permutations', '-n', default=1, help='number of generated permutations per SMILES'
)
@click.option('--fragment_bond', default='~', help='fragment bond token in the SMILES')
def cli(
    input_file_path: str, output_file_path: str, tokenize: bool, random_type: str,
    permutations: int, fragment_bond: str
) -> None:
    """The entry point for this cli script.

    Args:
        input (TextIO):  The input file (one SMILES per line).
        output (str): The output file (one SMILES per line).
    """

    # Create a instance of the Augmenter.
    ag = rrp.Augmenter.read_csv(input_file_path, fragment_bond)

    # Perform augmentation
    augm = ag.augment(rrp.RandomType[random_type], permutations, tokenize)

    # Exporting augmented samples
    augm[f'{random_type}'].to_csv(f'{random_type}.csv')


if __name__ == '__main__':
    cli()
