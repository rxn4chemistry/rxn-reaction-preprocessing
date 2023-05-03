# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2022
# ALL RIGHTS RESERVED
import pytest
from rxn.chemutils.reaction_equation import ReactionEquation

from rxn.reaction_preprocessing import MixedReactionFilter
from rxn.reaction_preprocessing.mixed_reaction_filter import ReactionFilterError


@pytest.fixture
def filter() -> MixedReactionFilter:
    return MixedReactionFilter(
        max_reactants=5,
        max_agents=0,
        max_products=1,
        min_reactants=2,
        min_agents=0,
        min_products=1,
        max_reactants_tokens=300,
        max_agents_tokens=0,
        max_products_tokens=200,
        max_absolute_formal_charge=2,
    )


@pytest.fixture
def good_reaction() -> ReactionEquation:
    return ReactionEquation.from_string(
        "O=[N+]([O-])c1cc(-c2nc3ccccc3o2)ccc1F.C~C>>Nc1cc(-c2nc3ccccc3o2)ccc1NCC(=O)N1CCOCC1",
        fragment_bond="~",
    )


@pytest.fixture
def bad_reaction() -> ReactionEquation:
    return ReactionEquation.from_string(
        "[C].C.[O--].[O--].O.O=[N+]([O-])c1cc(-c2nc3ccccc3o2)ccc1F.C.C>O>O.C"
    )


@pytest.fixture
def small_reaction() -> ReactionEquation:
    return ReactionEquation.from_string("C>O>")


@pytest.fixture
def big_reaction() -> ReactionEquation:
    return ReactionEquation.from_string(
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC>CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC>CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
    )


@pytest.fixture
def reaction_with_no_product() -> ReactionEquation:
    return ReactionEquation.from_string("CCC.CCCO>OC>")


@pytest.fixture
def alchemic_reaction() -> ReactionEquation:
    return ReactionEquation.from_string("C>[Hg]>[Au]")


def test_max_reactants_exceeded(
    filter: MixedReactionFilter,
    good_reaction: ReactionEquation,
    bad_reaction: ReactionEquation,
) -> None:
    assert not filter.max_reactants_exceeded(good_reaction)
    assert filter.max_reactants_exceeded(bad_reaction)


def test_max_agents_exceeded(
    filter: MixedReactionFilter,
    good_reaction: ReactionEquation,
    bad_reaction: ReactionEquation,
) -> None:
    assert not filter.max_agents_exceeded(good_reaction)
    assert filter.max_agents_exceeded(bad_reaction)


def test_max_products_exceeded(
    filter: MixedReactionFilter,
    good_reaction: ReactionEquation,
    bad_reaction: ReactionEquation,
) -> None:
    assert not filter.max_products_exceeded(good_reaction)
    assert filter.max_products_exceeded(bad_reaction)


def test_min_reactants_subceeded(
    filter: MixedReactionFilter,
    good_reaction: ReactionEquation,
    small_reaction: ReactionEquation,
) -> None:
    assert not filter.min_reactants_subceeded(good_reaction)
    assert filter.min_reactants_subceeded(small_reaction)


def test_min_agents_subceeded(
    filter: MixedReactionFilter,
    good_reaction: ReactionEquation,
    small_reaction: ReactionEquation,
) -> None:
    assert not filter.min_agents_subceeded(good_reaction)


def test_min_products_subceeded(
    filter: MixedReactionFilter,
    good_reaction: ReactionEquation,
    small_reaction: ReactionEquation,
) -> None:
    assert not filter.min_products_subceeded(good_reaction)
    assert filter.min_products_subceeded(small_reaction)


def test_products_subset_of_reactants(
    filter: MixedReactionFilter,
    good_reaction: ReactionEquation,
    bad_reaction: ReactionEquation,
    reaction_with_no_product: ReactionEquation,
) -> None:
    assert not filter.products_subset_of_reactants(reaction_with_no_product)
    assert not filter.products_subset_of_reactants(good_reaction)
    assert filter.products_subset_of_reactants(bad_reaction)


def test_products_single_atoms(
    filter: MixedReactionFilter,
    good_reaction: ReactionEquation,
    bad_reaction: ReactionEquation,
    reaction_with_no_product: ReactionEquation,
) -> None:
    assert not filter.products_single_atoms(reaction_with_no_product)
    assert not filter.products_single_atoms(good_reaction)
    assert filter.products_single_atoms(bad_reaction)


def test_max_reactant_tokens_exceeded(
    filter: MixedReactionFilter,
    good_reaction: ReactionEquation,
    big_reaction: ReactionEquation,
) -> None:
    assert not filter.max_reactant_tokens_exceeded(good_reaction)
    assert filter.max_reactant_tokens_exceeded(big_reaction)


def test_max_agent_tokens_exceeded(
    filter: MixedReactionFilter,
    good_reaction: ReactionEquation,
    big_reaction: ReactionEquation,
) -> None:
    assert not filter.max_agent_tokens_exceeded(good_reaction)
    assert filter.max_agent_tokens_exceeded(big_reaction)


def test_max_product_tokens_exceeded(
    filter: MixedReactionFilter,
    good_reaction: ReactionEquation,
    big_reaction: ReactionEquation,
) -> None:
    assert not filter.max_product_tokens_exceeded(good_reaction)
    assert filter.max_product_tokens_exceeded(big_reaction)


def test_formal_charge_exceeded(
    filter: MixedReactionFilter,
    good_reaction: ReactionEquation,
    bad_reaction: ReactionEquation,
) -> None:
    assert not filter.formal_charge_exceeded(good_reaction)
    assert filter.formal_charge_exceeded(bad_reaction)


def test_different_atom_types(
    filter: MixedReactionFilter,
    good_reaction: ReactionEquation,
    alchemic_reaction: ReactionEquation,
) -> None:
    assert not filter.different_atom_types(good_reaction)
    assert filter.different_atom_types(alchemic_reaction)

    # hydrogen should not count as a different atom type
    reaction = ReactionEquation.from_string(r"C=CC.CC.C>>[H]/C=C\C")
    assert not filter.different_atom_types(reaction)


def test_invalid_smiles(filter: MixedReactionFilter) -> None:
    # "[J]" is not a valid molecule
    reaction_with_invalid_smiles = ReactionEquation.from_string(
        "CCCCO.CCCCN.[J]>>CCCCOCCCCN"
    )
    assert not filter.is_valid(reaction_with_invalid_smiles)

    assert filter.validate_reasons(reaction_with_invalid_smiles) == (
        False,
        ["rdkit_molfromsmiles_failed"],
    )


def test_smiles_with_asterisks(filter: MixedReactionFilter) -> None:
    reactions_with_asterisks = [
        ReactionEquation.from_string("CCCCO.CCCCN>>CCCCOCCCCN*"),
        ReactionEquation.from_string("CCCCO.CCCCN*>>CCCCOCCCCN"),
        ReactionEquation.from_string("CCCCO.CCCCN*>>CCCCOCCCCN*"),
    ]
    for reaction in reactions_with_asterisks:
        assert not filter.is_valid(reaction)


def test_exception(
    filter: MixedReactionFilter,
    good_reaction: ReactionEquation,
    bad_reaction: ReactionEquation,
) -> None:
    # Nothing raised
    filter.validate(good_reaction)

    # exception raised
    with pytest.raises(ReactionFilterError):
        filter.validate(bad_reaction)
