import os
import sys
from pathlib import Path

import click

import rxn_reaction_preprocessing as rrp
import rxn_reaction_preprocessing.utils as utils


@click.command()
@click.argument('input_file_path', type=click.Path(exists=True), required=True)
@click.argument('output_file_path', type=click.Path(), required=True)
@click.option(
    '--fragment_bond',
    default='.',
    help='fragment bond token in the SMILES of the reactions to process'
)
def cli(input_file_path: str, output_file_path: str, fragment_bond: str) -> None:
    """The entry point for this cli script.

    Args:
        input_file_path (str):  The input file path (one SMILES per line).
        output_file_path (str): The output file path (one SMILES per line).
        fragment_bond (str): The fragment bond token used in the files to be standardized
    """

    # Create a instance of the Patterns.
    # for now jsonpath and fragment_bond (the one present in the jsonfile) fixed
    json_file_path = str(
        Path(utils.data_directory(), 'standardization-files/pistachio-200302.json')
    )
    print(json_file_path)
    pt = rrp.Patterns(json_file_path, fragment_bond='~')

    # Create an instance of the Standardizer
    std = rrp.Standardizer.read_csv(
        input_file_path,
        pt,
        reaction_column_name='rxn',
        fragment_bond=fragment_bond,
        kwargs={'lineterminator': '\n'}
    )
    # Perform standardization
    std.standardize()

    # Exporting standardized samples
    std.df.to_csv(output_file_path)


if __name__ == '__main__':
    cli()
