from typing import Optional

import click

from rxn.reaction_preprocessing import Standardizer
from rxn.reaction_preprocessing.annotations.molecule_annotation import (
    load_annotations_multiple,
)
from rxn.reaction_preprocessing.config import DEFAULT_ANNOTATION_FILES


@click.command()
@click.option("--csv_file", required=True)
@click.option(
    "--output_csv",
    help="(optional) Where to save the CSV augmented with column for missing annotations.",
)
@click.option(
    "--column_name", required=True, help="Column containing the reaction SMILES"
)
@click.option("--fragment_bond", default="~")
def main(
    csv_file: str, output_csv: Optional[str], column_name: str, fragment_bond: str
):
    """Find the missing annotation in a set of reactions.

    The missing annotations will be printed to standard output, and optionally
    a CSV will be created with the missing annotations in a new column.
    """

    # NB: in the future, the script may be updated to allow to change the annotation files.
    annotations = load_annotations_multiple(DEFAULT_ANNOTATION_FILES)

    # To find the missing annotations, we mis-use the standardizer, which anyway
    # looks for the missing annotations if `discard_unannotated_metals` is True.
    standardizer = Standardizer.read_csv(
        csv_file,
        annotations=annotations,
        discard_unannotated_metals=True,
        reaction_column_name=column_name,
        fragment_bond=fragment_bond,
        remove_stereo_if_not_defined_in_precursors=False,
    )

    standardizer.standardize()

    # Save csv if required
    if output_csv is not None:
        standardizer.df.to_csv(output_csv, index=False)

    missing_annotations = set()
    for missing_annotation_list in standardizer.df[
        standardizer.missing_annotations_column
    ]:
        for missing_annotation in missing_annotation_list:
            missing_annotations.add(missing_annotation)

    for missing_annotation in sorted(missing_annotations):
        print(missing_annotation)


if __name__ == "__main__":
    main()
