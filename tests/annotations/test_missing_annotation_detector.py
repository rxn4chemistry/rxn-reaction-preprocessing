from typing import List, Tuple

from rxn_chemutils.reaction_equation import ReactionEquation

from rxn_reaction_preprocessing.annotations.missing_annotation_detector import (
    MissingAnnotationDetector,
)
from rxn_reaction_preprocessing.annotations.molecule_annotation import (
    MoleculeAnnotation,
)


def mock_requires_annotation(smiles: str) -> bool:
    """Mock callback to use as an annotation criterion, returns True only
    when "1" is in the SMILES string."""
    return "1" in smiles


def test_molecule_needs_annotation():
    annotated_molecules = {
        "CCCC",
        "C1CC1",
        "c1ccccc1",
    }
    mad = MissingAnnotationDetector(annotated_molecules, mock_requires_annotation)

    # random molecules not fulfilling the annotation criterion
    assert not mad.molecule_needs_annotation("C[Zn]C")
    assert not mad.molecule_needs_annotation("CO")

    # random molecules fulfilling the annotation criterion
    assert mad.molecule_needs_annotation("O1CCC1")
    assert mad.molecule_needs_annotation("c1cc(CC)ccc1")

    # random molecules fulfilling the annotation criterion but already annotated
    assert not mad.molecule_needs_annotation("C1CC1")
    assert not mad.molecule_needs_annotation("c1ccccc1")


def test_missing_annotations_in_reaction_equations():
    annotated_molecules = {
        "CCCC",
        "C1CC1",
        "c1ccccc1",
    }
    mad = MissingAnnotationDetector(annotated_molecules, mock_requires_annotation)

    reactions_and_expected_annotations: List[Tuple[str, List[str]]] = [
        # nothing
        ("C.CC>>CCC", []),
        # already annotated
        ("C1CC1>>CCC", []),
        # one reactant and product must be annotated
        ("C.C1NC1>O>CCC1CC1", ["C1NC1", "CCC1CC1"]),
        # one molecule with fragment bond must be annotated
        ("C.N~CCC1NC1>O>CCCC", ["N.CCC1NC1"]),
    ]

    # Get the required annotations for reactions, one at a time
    for reaction_smiles, expected in reactions_and_expected_annotations:
        reaction = ReactionEquation.from_string(reaction_smiles, fragment_bond="~")
        assert list(mad.missing_in_reaction_equation(reaction)) == expected

    # Get the required annotations for a set of reactions, all at the same time
    reactions = [
        ReactionEquation.from_string(reaction_smiles, "~")
        for reaction_smiles, _ in reactions_and_expected_annotations
    ]
    expected_molecules = [
        expected_molecule
        for _, expected_list in reactions_and_expected_annotations
        for expected_molecule in expected_list
    ]
    assert list(mad.missing_in_reaction_equations(reactions)) == expected_molecules


def test_missing_annotations_in_reaction_smiles():
    annotated_molecules = {
        "CCCC",
        "C1CC1",
        "c1ccccc1",
    }
    mad = MissingAnnotationDetector(annotated_molecules, mock_requires_annotation)

    # nothing to annotate
    assert list(mad.missing_in_reaction_smiles("C.CC>>CCC")) == []
    # already annotated
    assert list(mad.missing_in_reaction_smiles("C1CC1>>CCC")) == []
    # one reactant and product must be annotated
    assert list(mad.missing_in_reaction_smiles("C.C1NC1>O>CCC1CC1")) == [
        "C1NC1",
        "CCC1CC1",
    ]
    # one molecule with fragment bond must be annotated
    assert list(mad.missing_in_reaction_smiles("C.N~CCC1NC1>O>CCCC", "~")) == [
        "N.CCC1NC1"
    ]

    # Multiple reaction SMILES at a time
    assert list(
        mad.missing_in_reaction_smiles(
            ["C1CC1>>CCC", "C.C1NC1>O>CCC1CC1", "C.N~CCC1NC1>O>CCCC"], "~"
        )
    ) == ["C1NC1", "CCC1CC1", "N.CCC1NC1"]


def test_default_required_annotations():
    # If the callback is not provided, falls back to the extended transition metals
    annotated_molecules = {
        "C[Pd]C",
    }
    mad = MissingAnnotationDetector(annotated_molecules)

    assert not mad.molecule_needs_annotation("CCOCCOCS")
    assert not mad.molecule_needs_annotation("C[Pd]C")
    assert mad.molecule_needs_annotation("CC[Pd]CC")
    assert mad.molecule_needs_annotation("[Zn+2].[Mg+2]")
    assert mad.molecule_needs_annotation("C[Al](C)C")


def test_from_annotations():
    annotations = [
        MoleculeAnnotation("C[Pd]C", None, "accept", []),
        MoleculeAnnotation("O[Pd]O", None, "reject", []),
        MoleculeAnnotation("N[Pd]N", "[NH3+][Pd-2][NH3+]", "accept", []),
    ]

    # NB: use the default AnnotationCriterion - i.e. extended transition metals
    mad = MissingAnnotationDetector.from_molecule_annotations(annotations)

    # Transition metals with no annotation yet -> requires annotation
    assert mad.molecule_needs_annotation("[Pd+].[Cl-]")
    assert mad.molecule_needs_annotation("Cl[Zn]Cl")

    # No transition metal
    assert not mad.molecule_needs_annotation("COC")

    # Already in annotations
    assert not mad.molecule_needs_annotation("C[Pd]C")
    assert not mad.molecule_needs_annotation("N[Pd]N")

    # Also the updated SMILES of an annotation does not require annotation
    assert not mad.molecule_needs_annotation("[NH3+][Pd-2][NH3+]")
