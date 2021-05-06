from pathlib import Path

import pytest

from rxn_reaction_preprocessing.annotations.molecule_annotation import AnnotationDecision
from rxn_reaction_preprocessing.annotations.molecule_annotation import load_annotations
from rxn_reaction_preprocessing.annotations.molecule_annotation import MoleculeAnnotation


def test_molecule_annotation():
    # basic instantiation check
    annotation = MoleculeAnnotation(
        original_smiles='CC~O~N',
        updated_smiles='CC.O~N',
        decision='accept',
        categories=[],
        dummy_a=2,
        dummy_b='a'
    )
    assert annotation.original_smiles == 'CC~O~N'
    assert annotation.original_without_fragment_bond == 'CC.O.N'
    assert annotation.updated_smiles == 'CC.O~N'
    assert annotation.updated_without_fragment_bond == ['CC', 'O.N']
    assert annotation.decision is AnnotationDecision.ACCEPT
    assert annotation.extra_info['dummy_a'] == 2
    assert annotation.extra_info['dummy_b'] == 'a'

    # raises if invalid decision
    with pytest.raises(ValueError):
        _ = MoleculeAnnotation(
            original_smiles='CC~O~N',
            updated_smiles='CC.O~N',
            decision='invalid_decision',
            categories=[],
        )

    # Raises if getting updated SMILES list when there is no updated SMILES
    annotation.updated_smiles = None
    with pytest.raises(RuntimeError):
        _ = annotation.updated_without_fragment_bond


def test_import_from_file():
    annotation_file = str(Path(__file__).parent / 'test_molecule_annotations.json')
    annotations = load_annotations(annotation_file)

    # Check a few loaded values
    first = annotations[0]
    assert first.original_smiles == 'CC[Zn]CC~Cc1ccccc1'
    assert first.updated_smiles == 'CC[Zn]CC.Cc1ccccc1'
    assert first.categories == ['Zn']

    third = annotations[2]
    assert third.original_smiles == 'CCO'
    assert third.updated_smiles is None
    assert third.categories == []
    assert third.decision is AnnotationDecision.ACCEPT

    fourth = annotations[3]
    assert fourth.decision is AnnotationDecision.REJECT
