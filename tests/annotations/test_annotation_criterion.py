import pytest
from rxn.chemutils.exceptions import InvalidSmiles

from rxn.reaction_preprocessing.annotations.annotation_criterion import (
    AnnotationCriterion,
)


def test_extended_transition_metals():
    # Test the function to get the extended transition metals
    elements_to_check = set(AnnotationCriterion.extended_transition_metals())

    assert "Al" in elements_to_check
    assert "Pb" in elements_to_check
    assert "Sc" in elements_to_check
    assert "Hg" in elements_to_check
    assert "Pd" in elements_to_check
    assert "La" in elements_to_check
    assert "Er" in elements_to_check
    assert "U" in elements_to_check

    assert "C" not in elements_to_check
    assert "Na" not in elements_to_check
    assert "Cl" not in elements_to_check
    assert "Sn" not in elements_to_check


def test_annotation_criterion():
    # Here we indicate that in addition to the transition metals, phosphorus
    # requires an annotation, and also that Pd does not require annotations.
    annotation_criterion = AnnotationCriterion(
        additional_elements_to_consider=["P"], elements_not_to_consider=["Pd"]
    )

    # Not requiring annotations: anything without transition metals
    assert not annotation_criterion("CCOCC")
    assert not annotation_criterion("CCOCC.[Na+].[Cl-]")
    assert not annotation_criterion("CI")

    # Not requiring annotations: Pd (because specified above)
    assert not annotation_criterion("CC[Pd]OC(=O)C")

    # Not requiring annotations: Sn
    assert not annotation_criterion("CCCC[Sn](CCCC)(CCCC)Cl")

    # Requiring annotation: other transition metals
    assert annotation_criterion("CC[Ni]OC(=O)C")
    assert annotation_criterion("[Co]")

    # Requiring annotation: phosphorus (because specified above)
    assert annotation_criterion("P(C)(C)C")

    # Raises for invalid SMILES
    with pytest.raises(InvalidSmiles):
        _ = annotation_criterion("CC((N")
