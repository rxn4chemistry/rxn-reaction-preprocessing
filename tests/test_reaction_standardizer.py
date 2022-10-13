import pytest
from rxn.chemutils.reaction_equation import ReactionEquation

from rxn.reaction_preprocessing.reaction_standardizer import ReactionStandardizer


@pytest.fixture
def standardizer() -> ReactionStandardizer:
    return ReactionStandardizer()


def test_merges_reactants_and_reagents(standardizer: ReactionStandardizer) -> None:
    reaction = ReactionEquation.from_string("A.B.C>D.E>F.G")
    assert standardizer(reaction).to_string() == "A.B.C.D.E>>F.G"


def test_removes_precursors_from_product(standardizer: ReactionStandardizer) -> None:
    # in the reactants
    reaction = ReactionEquation.from_string("A.B.C>D.E>F.A.G")
    assert standardizer(reaction).to_string() == "A.B.C.D.E>>F.G"

    # in the agents
    reaction = ReactionEquation.from_string("A.B.C>D.E>F.E.G")
    assert standardizer(reaction).to_string() == "A.B.C.D.E>>F.G"

    # do not remove if part of a fragment bond on either side:
    # - F is not removed because it is part of a fragment in the precursors
    # - G~A is not removed because only A is a precursor
    reaction = ReactionEquation.from_string("A.B~F>>F.G~A", "~")
    assert standardizer(reaction).to_string("~") == "A.B~F>>F.G~A"

    # But remove if the full molecule with fragment matches
    reaction = ReactionEquation.from_string("A.B~F>>B~F.G~A", "~")
    assert standardizer(reaction).to_string("~") == "A.B~F>>G~A"


def test_remove_duplicate_compounds(standardizer: ReactionStandardizer) -> None:
    reaction = ReactionEquation.from_string("A.B.C.A.D>B>F.G.F")
    assert standardizer(reaction).to_string() == "A.B.C.D>>F.G"

    # With fragment bonds - removes only if full compound matches
    reaction = ReactionEquation.from_string("A~B.C~D.A.B.C~D>>E", "~")
    assert standardizer(reaction).to_string("~") == "A.A~B.B.C~D>>E"


def test_sorts_the_compounds(standardizer: ReactionStandardizer) -> None:
    reaction = ReactionEquation.from_string("N.A>D.E>M.F")
    assert standardizer(reaction).to_string() == "A.D.E.N>>F.M"


def test_does_not_modify_original_reaction(standardizer: ReactionStandardizer) -> None:
    # example with multiple things that will be changed
    rxn_string = "D.B.A>C>D.E"
    expected = "A.B.C.D>>E"

    original_reaction = ReactionEquation.from_string(rxn_string)
    updated_reaction = standardizer(original_reaction)

    # The updated reaction has the expected string, while the original one did not change
    assert original_reaction.to_string() == rxn_string
    assert updated_reaction.to_string() == expected
