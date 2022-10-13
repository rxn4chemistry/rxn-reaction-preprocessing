from pathlib import Path

import pytest

from rxn.reaction_preprocessing.annotations.molecule_annotation import (
    AnnotationDecision,
    MoleculeAnnotation,
    load_annotations,
    load_annotations_multiple,
)


def test_molecule_annotation() -> None:
    # basic instantiation check
    annotation = MoleculeAnnotation(
        original_smiles="CC~O~N",
        updated_smiles="CC.O~N",
        decision="accept",
        categories=[],
        dummy_a=2,
        dummy_b="a",
    )
    assert annotation.original_smiles == "CC~O~N"
    assert annotation.original_without_fragment_bond == "CC.O.N"
    assert annotation.updated_smiles == "CC.O~N"
    assert annotation.updated_without_fragment_bond == ["CC", "O.N"]
    assert annotation.decision is AnnotationDecision.ACCEPT
    assert annotation.extra_info["dummy_a"] == 2
    assert annotation.extra_info["dummy_b"] == "a"

    # raises if invalid decision
    with pytest.raises(ValueError):
        _ = MoleculeAnnotation(
            original_smiles="CC~O~N",
            updated_smiles="CC.O~N",
            decision="invalid_decision",
            categories=[],
        )

    # Raises if getting updated SMILES list when there is no updated SMILES
    annotation.updated_smiles = None
    with pytest.raises(RuntimeError):
        _ = annotation.updated_without_fragment_bond


def test_import_from_file() -> None:
    annotation_file = str(Path(__file__).parent / "test_molecule_annotations.json")
    annotations = load_annotations(annotation_file)

    # Check a few loaded values
    first = annotations[0]
    assert first.original_smiles == "CC[Zn]CC~Cc1ccccc1"
    assert first.updated_smiles == "CC[Zn]CC.Cc1ccccc1"
    assert first.categories == ["Zn"]

    third = annotations[2]
    assert third.original_smiles == "CCO"
    assert third.updated_smiles is None
    assert third.categories == []
    assert third.decision is AnnotationDecision.ACCEPT

    fourth = annotations[3]
    assert fourth.decision is AnnotationDecision.REJECT


def test_import_from_multiple_files() -> None:
    annotation_file = str(Path(__file__).parent / "test_molecule_annotations.json")

    # Loading 1 file only delivers the same if called from path or from list of paths
    load_single = load_annotations(annotation_file)
    load_multiple_1 = load_annotations_multiple([annotation_file])
    assert load_single == load_multiple_1

    # No removal of duplicates; loading from two files is the same as loading separately
    load_multiple_2 = load_annotations_multiple([annotation_file, annotation_file])
    assert load_multiple_2 == load_single + load_single

    # Empty list
    assert load_annotations_multiple([]) == []
