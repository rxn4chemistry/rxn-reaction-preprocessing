import pytest
from rxn.chemutils.conversion import canonicalize_smiles
from rxn.chemutils.exceptions import InvalidSmiles
from rxn.chemutils.reaction_equation import ReactionEquation

from rxn.reaction_preprocessing.annotations.molecule_annotation import (
    MoleculeAnnotation,
)
from rxn.reaction_preprocessing.cleaner import remove_isotope_information
from rxn.reaction_preprocessing.molecule_standardizer import (
    MissingAnnotation,
    MoleculeStandardizer,
    RejectedMolecule,
)


def test_callable() -> None:
    standardizer = MoleculeStandardizer()

    smiles_strings = ["C(C)", "CC", "C1=CC=CC=C1"]

    # calling the object is the same as calling the standardize() function
    for smiles in smiles_strings:
        assert standardizer(smiles) == standardizer.standardize(smiles)


def test_tilde_as_fragment_bond() -> None:
    standardizer = MoleculeStandardizer()

    fragment_ok = "C.O"
    fragment_not_ok = "C~O"

    # Fragment bond with a dot -> all ok
    assert standardizer(fragment_ok) == [fragment_ok]

    # Fragment bond with a tilde -> not permitted
    with pytest.raises(ValueError) as exc_info:
        _ = standardizer(fragment_not_ok)

    # It must not be a StandardizationError
    assert exc_info.type is ValueError


def test_invalid_smiles() -> None:
    standardizer = MoleculeStandardizer()

    smiles = "C(C)(C)(C)"
    invalid_smiles = "Invalid"

    # Default behavior: canonicalization
    assert standardizer.standardize(smiles) != [smiles]
    assert standardizer.standardize(smiles) == [canonicalize_smiles(smiles)]

    # Invalid SMILES
    with pytest.raises(InvalidSmiles):
        _ = standardizer.standardize(invalid_smiles)

    # non-canonicalizing standardizer - still raises for invalid SMILES
    assert MoleculeStandardizer(canonicalize=False).standardize(smiles) == [smiles]
    with pytest.raises(InvalidSmiles):
        _ = MoleculeStandardizer(canonicalize=False).standardize(invalid_smiles)


def test_annotated_as_rejected() -> None:
    rejected_smiles = "CCC"
    annotations = [
        MoleculeAnnotation(
            original_smiles=rejected_smiles,
            updated_smiles=None,
            decision="reject",
            categories=[],
        )
    ]

    standardizer_without_annotations = MoleculeStandardizer()
    standardizer_with_annotations = MoleculeStandardizer(annotations=annotations)

    assert standardizer_without_annotations(rejected_smiles) == [rejected_smiles]
    with pytest.raises(RejectedMolecule):
        _ = standardizer_with_annotations(rejected_smiles)


def test_missing_annotation() -> None:
    # transition metals require an annotation
    tm_smiles = "[C][Pd][C]"
    annotations = [
        MoleculeAnnotation(
            original_smiles=tm_smiles,
            updated_smiles=None,
            decision="accept",
            categories=[],
        )
    ]

    # By default: unannotated are kept
    nondiscarding_standardizer = MoleculeStandardizer()
    assert nondiscarding_standardizer(tm_smiles) == [tm_smiles]

    # with the flag, no annotation: exception raised
    discarding_standardizer = MoleculeStandardizer(discard_missing_annotations=True)
    with pytest.raises(MissingAnnotation):
        _ = discarding_standardizer(tm_smiles)

    # with the flag, molecule annotated: ok
    standardizer_with_annotations = MoleculeStandardizer(
        annotations=annotations, discard_missing_annotations=True
    )
    assert standardizer_with_annotations(tm_smiles) == [tm_smiles]


def test_replacements() -> None:
    annotations = [
        MoleculeAnnotation(
            original_smiles="CC", updated_smiles=None, decision="accept", categories=[]
        ),
        MoleculeAnnotation(
            original_smiles="CCC",
            updated_smiles="CC~C",
            decision="accept",
            categories=[],
        ),
        MoleculeAnnotation(
            original_smiles="CCCC",
            updated_smiles="CC.CC",
            decision="accept",
            categories=[],
        ),
    ]

    standardizer = MoleculeStandardizer(annotations=annotations)

    # Annotated molecule with no replacement
    assert standardizer("CC") == ["CC"]

    # Annotated molecule, replacement has fragment
    assert standardizer("CCC") == ["CC.C"]

    # Annotated molecule, replacement is two separate molecules
    assert standardizer("CCCC") == ["CC", "CC"]


def test_canonicalization_before_rejection_check() -> None:
    non_canonical_smiles = "CC(C)"
    rejected_smiles = "CCC"
    annotations = [
        MoleculeAnnotation(
            original_smiles=rejected_smiles,
            updated_smiles=None,
            decision="reject",
            categories=[],
        )
    ]

    # The rejection check must be done after canonicalization
    with pytest.raises(RejectedMolecule):
        _ = MoleculeStandardizer(annotations=annotations)(non_canonical_smiles)

    # Without canonicalization, no exception: the molecule is not recognized as a rejected one
    _ = MoleculeStandardizer(annotations=annotations, canonicalize=False)(
        non_canonical_smiles
    )


def test_missing_annotation_check_after_rejection_check() -> None:
    # The RejectedMolecule exception has priority over the MissingAnnotation one.
    rejected_smiles = "[C][Pd][C]"
    annotations = [
        MoleculeAnnotation(
            original_smiles=rejected_smiles,
            updated_smiles=None,
            decision="reject",
            categories=[],
        )
    ]

    standardizer = MoleculeStandardizer(
        annotations=annotations, discard_missing_annotations=True
    )
    with pytest.raises(RejectedMolecule):
        _ = standardizer(rejected_smiles)

    # Removing the annotations, it is now a MissingAnnotation being raised
    standardizer = MoleculeStandardizer(
        annotations=[], discard_missing_annotations=True
    )
    with pytest.raises(MissingAnnotation):
        _ = standardizer(rejected_smiles)


def test_discarding_isotope_info() -> None:
    # As a check: doing the isotope info removal after canonicalization does
    # not always lead to a canonical SMILES
    assert remove_isotope_information(
        canonicalize_smiles("[14CH3]COC")
    ) != canonicalize_smiles("CCOC")

    # On the other hand, in the molecule standardizer, the standardized
    # molecules must be canonical even if there was some isotope information.
    standardizer = MoleculeStandardizer()
    assert standardizer.standardize("[14CH3]COC") == [canonicalize_smiles("CCOC")]


def test_standardize_in_equation() -> None:
    annotations = [
        MoleculeAnnotation(
            original_smiles="CC", updated_smiles=None, decision="accept", categories=[]
        ),
        MoleculeAnnotation(
            original_smiles="CCC",
            updated_smiles="CC~C",
            decision="accept",
            categories=[],
        ),
        MoleculeAnnotation(
            original_smiles="CCCC",
            updated_smiles="CC.CC",
            decision="accept",
            categories=[],
        ),
        MoleculeAnnotation(
            original_smiles="C", updated_smiles=None, decision="reject", categories=[]
        ),
    ]

    # Standardizes following the annotations, and additionally canonicalizes OC -> CO
    standardizer = MoleculeStandardizer(annotations=annotations)
    reaction = ReactionEquation.from_string("CC.CCC.O>CCCC>OC")
    assert (
        standardizer.standardize_in_equation(reaction).to_string("~")
        == "CC.CC~C.O>CC.CC>CO"
    )

    # Rejected because 'C' is rejected from the annotations
    rejected_reaction = ReactionEquation.from_string("C.O>>CO")
    with pytest.raises(RejectedMolecule):
        standardizer.standardize_in_equation(rejected_reaction)
