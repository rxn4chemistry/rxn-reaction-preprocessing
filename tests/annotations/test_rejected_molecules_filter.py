from rxn_reaction_preprocessing.annotations.molecule_annotation import (
    MoleculeAnnotation,
)
from rxn_reaction_preprocessing.annotations.rejected_molecules_filter import (
    RejectedMoleculesFilter,
)


def test_molecule_filter_with_molecule_smiles():
    rejected_molecules = ["O", "C[O-].[Na+]"]
    molecule_filter = RejectedMoleculesFilter(rejected_molecules)

    # Reject when given the exact same SMILES
    assert not molecule_filter.is_valid_molecule_smiles("O")
    assert not molecule_filter.is_valid_molecule_smiles("C[O-].[Na+]")

    # No reject when part of a bigger complex
    assert molecule_filter.is_valid_molecule_smiles("O.O")

    # No reject when tilde is used instead of dot
    assert molecule_filter.is_valid_molecule_smiles("C[O-]~[Na+]")

    # Does not try to re-order the fragments
    assert molecule_filter.is_valid_molecule_smiles("[Na+].C[O-]")


def test_molecule_filter_with_reaction_smiles():
    rejected_molecules = ["O", "C[O-].[Na+]"]
    molecule_filter = RejectedMoleculesFilter(rejected_molecules)

    # Reject when given the exact same SMILES
    assert not molecule_filter.is_valid_reaction_smiles("O.CC>>CCO")
    assert not molecule_filter.is_valid_reaction_smiles("CC>O>CCO")

    # Reject when fragments are matching (but only if correct fragment bond is given)
    assert not molecule_filter.is_valid_reaction_smiles(
        "C.C[O-]~[Na+]>>CCO", fragment_bond="~"
    )
    assert not molecule_filter.is_valid_reaction_smiles(
        "C.C[O-]%[Na+]>>CCO", fragment_bond="%"
    )
    assert molecule_filter.is_valid_reaction_smiles(
        "C.C[O-].[Na+]>>CCO", fragment_bond="~"
    )
    assert molecule_filter.is_valid_reaction_smiles("C.C[O-]~[Na+]>>CCO")

    # No reject when part of a bigger complex
    assert molecule_filter.is_valid_molecule_smiles("O.O")

    # No reject when tilde is used instead of dot
    assert molecule_filter.is_valid_molecule_smiles("C[O-]~[Na+]")


def test_molecule_filter_from_annotations():
    annotations = [
        MoleculeAnnotation("O", None, "accept", []),
        MoleculeAnnotation("C[O-]~[Na+]", None, "reject", []),
    ]

    molecule_filter = RejectedMoleculesFilter.from_molecule_annotations(annotations)

    # Does not reject 'O', which has a decision "accept" above
    assert molecule_filter.is_valid_molecule_smiles("O")

    # Rejects 'C[O-].[Na+]' but only if not using the fragment bond at that stage
    assert not molecule_filter.is_valid_molecule_smiles("C[O-].[Na+]")
    assert molecule_filter.is_valid_molecule_smiles("C[O-]~[Na+]")
