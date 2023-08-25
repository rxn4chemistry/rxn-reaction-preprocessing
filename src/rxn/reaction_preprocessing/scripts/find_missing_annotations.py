import click
from rxn.utilities.csv import iterate_csv_column

from rxn.reaction_preprocessing.annotations.molecule_annotation import (
    load_annotations_multiple,
)
from rxn.reaction_preprocessing.config import DEFAULT_ANNOTATION_FILES
from rxn.reaction_preprocessing.standardizer import Standardizer


@click.command()
@click.option("--csv_file", required=True)
@click.option(
    "--column_name", required=True, help="Column containing the reaction SMILES"
)
@click.option("--fragment_bond", default="~")
def main(csv_file: str, column_name: str, fragment_bond: str) -> None:
    """Find the missing annotation in a set of reactions.

    The missing annotations will be printed to standard output, and optionally
    a CSV will be created with the missing annotations in a new column.
    """

    # NB: in the future, the script may be updated to allow to change the annotation files.
    annotations = load_annotations_multiple(DEFAULT_ANNOTATION_FILES)

    # To find the missing annotations, we mis-use the standardizer, which anyway
    # looks for the missing annotations if `discard_unannotated_metals` is True.
    standardizer = Standardizer(
        annotations=annotations,
        discard_unannotated_metals=True,
        reaction_column_name=column_name,
        fragment_bond=fragment_bond,
        remove_stereo_if_not_defined_in_precursors=False,
    )

    missing_annotations = set()
    input_reactions = iterate_csv_column(csv_file, column_name)
    for rxn_smiles in input_reactions:
        result = standardizer.standardize_one(rxn_smiles)
        missing_annotations.update(result.missing_annotations)

    for missing_annotation in sorted(missing_annotations):
        print(missing_annotation)


if __name__ == "__main__":
    main()
