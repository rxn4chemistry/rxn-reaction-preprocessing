import pytest

from rxn.reaction_preprocessing.annotations.annotation_info import (
    AnnotationInfo,
    MoleculeNotAnnotated,
)
from rxn.reaction_preprocessing.annotations.molecule_annotation import (
    MoleculeAnnotation,
)

annotations = [
    MoleculeAnnotation("O", None, "accept", []),
    MoleculeAnnotation("(C)(C)", "CC", "reject", []),
    MoleculeAnnotation("[C]~CCCC~[Pd]", "[C]~[Pd].CCCC", "accept", []),
]


def test_is_annotated() -> None:
    annotation_info = AnnotationInfo(annotations)

    # The thee given molecules
    assert annotation_info.is_annotated("O")
    assert annotation_info.is_annotated("(C)(C)")
    assert annotation_info.is_annotated("[C].CCCC.[Pd]")

    # updated SMILES - not directly annotated
    assert not annotation_info.is_annotated("CC")
    # Requires molecules without fragment bond
    assert not annotation_info.is_annotated("[C]~CCCC~[Pd]")


def test_is_accepted() -> None:
    annotation_info = AnnotationInfo(annotations)

    # The thee given molecules
    assert annotation_info.is_accepted("O")
    assert not annotation_info.is_accepted("(C)(C)")
    assert annotation_info.is_accepted("[C].CCCC.[Pd]")

    # updated SMILES - not directly annotated
    with pytest.raises(MoleculeNotAnnotated):
        _ = annotation_info.is_accepted("CC")
    # Requires molecules without fragment bond
    with pytest.raises(MoleculeNotAnnotated):
        _ = annotation_info.is_accepted("[C]~CCCC~[Pd]")


def test_is_rejected() -> None:
    annotation_info = AnnotationInfo(annotations)

    # The thee given molecules
    assert not annotation_info.is_rejected("O")
    assert annotation_info.is_rejected("(C)(C)")
    assert not annotation_info.is_rejected("[C].CCCC.[Pd]")

    # updated SMILES - not directly annotated
    with pytest.raises(MoleculeNotAnnotated):
        _ = annotation_info.is_rejected("CC")
    # Requires molecules without fragment bond
    with pytest.raises(MoleculeNotAnnotated):
        _ = annotation_info.is_rejected("[C]~CCCC~[Pd]")
