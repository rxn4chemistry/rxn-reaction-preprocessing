# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
import pytest
from rxn_chemutils.reaction_equation import ReactionEquation

from rxn_reaction_preprocessing.mixed_reaction_filter import (
    MixedReactionFilter, ReactionFilterError
)


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
        'O=[N+]([O-])c1cc(-c2nc3ccccc3o2)ccc1F.C~C>>Nc1cc(-c2nc3ccccc3o2)ccc1NCC(=O)N1CCOCC1',
        fragment_bond='~'
    )


@pytest.fixture
def bad_reaction() -> ReactionEquation:
    return ReactionEquation.from_string(
        '[C].C.[O--].[O--].O.O=[N+]([O-])c1cc(-c2nc3ccccc3o2)ccc1F.C.C>O>O.C'
    )


@pytest.fixture
def small_reaction() -> ReactionEquation:
    return ReactionEquation.from_string('C>O>')


@pytest.fixture
def big_reaction() -> ReactionEquation:
    return ReactionEquation.from_string(
        'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC>CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC>CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
    )


@pytest.fixture
def reaction_with_no_product() -> ReactionEquation:
    return ReactionEquation.from_string('CCC.CCCO>OC>')


@pytest.fixture
def alchemic_reaction() -> ReactionEquation:
    return ReactionEquation.from_string('C>[Hg]>[Au]')


def test_max_reactants_exceeded(filter, good_reaction, bad_reaction):
    assert not filter.max_reactants_exceeded(good_reaction)
    assert filter.max_reactants_exceeded(bad_reaction)


def test_max_agents_exceeded(filter, good_reaction, bad_reaction):
    assert not filter.max_agents_exceeded(good_reaction)
    assert filter.max_agents_exceeded(bad_reaction)


def test_max_products_exceeded(filter, good_reaction, bad_reaction):
    assert not filter.max_products_exceeded(good_reaction)
    assert filter.max_products_exceeded(bad_reaction)


def test_min_reactants_subceeded(filter, good_reaction, small_reaction):
    assert not filter.min_reactants_subceeded(good_reaction)
    assert filter.min_reactants_subceeded(small_reaction)


def test_min_agents_subceeded(filter, good_reaction, small_reaction):
    assert not filter.min_agents_subceeded(good_reaction)


def test_min_products_subceeded(filter, good_reaction, small_reaction):
    assert not filter.min_products_subceeded(good_reaction)
    assert filter.min_products_subceeded(small_reaction)


def test_products_subset_of_reactants(
    filter, good_reaction, bad_reaction, reaction_with_no_product
):
    assert not filter.products_subset_of_reactants(reaction_with_no_product)
    assert not filter.products_subset_of_reactants(good_reaction)
    assert filter.products_subset_of_reactants(bad_reaction)


def test_products_single_atoms(filter, good_reaction, bad_reaction, reaction_with_no_product):
    assert not filter.products_single_atoms(reaction_with_no_product)
    assert not filter.products_single_atoms(good_reaction)
    assert filter.products_single_atoms(bad_reaction)


def test_max_reactant_tokens_exceeded(filter, good_reaction, big_reaction):
    assert not filter.max_reactant_tokens_exceeded(good_reaction)
    assert filter.max_reactant_tokens_exceeded(big_reaction)


def test_max_agent_tokens_exceeded(filter, good_reaction, big_reaction):
    assert not filter.max_agent_tokens_exceeded(good_reaction)
    assert filter.max_agent_tokens_exceeded(big_reaction)


def test_max_product_tokens_exceeded(filter, good_reaction, big_reaction):
    assert not filter.max_product_tokens_exceeded(good_reaction)
    assert filter.max_product_tokens_exceeded(big_reaction)


def test_formal_charge_exceeded(filter, good_reaction, bad_reaction):
    assert not filter.formal_charge_exceeded(good_reaction)
    assert filter.formal_charge_exceeded(bad_reaction)


def test_different_atom_types(filter, good_reaction, alchemic_reaction):
    assert not filter.different_atom_types(good_reaction)
    assert filter.different_atom_types(alchemic_reaction)


def test_invalid_smiles(filter: MixedReactionFilter):
    # "[J]" is not a valid molecule
    reaction_with_invalid_smiles = ReactionEquation.from_string('CCCCO.CCCCN.[J]>>CCCCOCCCCN')
    assert not filter.is_valid(reaction_with_invalid_smiles)

    assert filter.validate_reasons(reaction_with_invalid_smiles
                                   ) == (False, ['rdkit_molfromsmiles_failed'])


def test_exception(filter, good_reaction, bad_reaction):
    # Nothing raised
    _ = filter.validate(good_reaction)

    # exception raised
    with pytest.raises(ReactionFilterError):
        _ = filter.validate(bad_reaction)
