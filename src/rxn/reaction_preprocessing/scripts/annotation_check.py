import csv
from typing import Generator

import click

from rxn.reaction_preprocessing.annotations.annotation_info import AnnotationInfo
from rxn.reaction_preprocessing.annotations.missing_annotation_detector import (
    MissingAnnotationDetector,
)
from rxn.reaction_preprocessing.annotations.molecule_annotation import (
    load_annotations_multiple,
)
from rxn.reaction_preprocessing.config import DEFAULT_ANNOTATION_FILES


def iterate_rxn_smiles(csv_file: str, column_name: str) -> Generator[str, None, None]:
    with open(csv_file) as f:
        r = csv.reader(f)
        header = next(r)
        try:
            smiles_index = header.index(column_name)
        except ValueError as e:
            raise RuntimeError(f'No "{column_name}" column in {csv_file}') from e
        for row in r:
            yield row[smiles_index]


@click.command()
@click.option("--csv_file", required=True)
@click.option(
    "--column_name", required=True, help="Column containing the reaction SMILES"
)
def main(csv_file: str, column_name: str) -> None:
    """Check for missing annotations: what is already annotated (accepted /
    rejected), what still needs to be annotated."""

    iterator = iterate_rxn_smiles(csv_file, column_name)

    missing_annotation_detector = MissingAnnotationDetector(set())
    molecules_requiring_annotation = list(
        missing_annotation_detector.missing_in_reaction_smiles(
            iterator, fragment_bond="~"
        )
    )

    annotations = load_annotations_multiple(DEFAULT_ANNOTATION_FILES)
    annotation_info = AnnotationInfo(annotations)

    not_annotated = [
        m for m in molecules_requiring_annotation if not annotation_info.is_annotated(m)
    ]
    annotated = [
        m for m in molecules_requiring_annotation if annotation_info.is_annotated(m)
    ]
    accepted = [m for m in annotated if annotation_info.is_accepted(m)]
    rejected = [m for m in annotated if annotation_info.is_rejected(m)]

    to_print = [
        ("requiring annotation", molecules_requiring_annotation),
        ("not annotated", not_annotated),
        ("annotated", annotated),
        ("accepted", accepted),
        ("rejected", rejected),
    ]

    # Print summary
    for label, smiles_list in to_print:
        print(label, len(smiles_list), len(set(smiles_list)))

    # Print details
    for label, smiles_list in to_print:
        print()
        print(label)
        print("=" * len(label))
        for smiles in sorted(set(smiles_list)):
            print(smiles.replace("~", "."))


if __name__ == "__main__":
    main()
