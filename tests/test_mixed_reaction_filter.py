import pytest
from data_preprocessor import MixedReactionFilter, Reaction


@pytest.fixture
def filter():
    return MixedReactionFilter(
        max_reactants=5,
        max_products=1,
        max_agents=0,
        max_reactants_tokens=300,
        max_agents_tokens=0,
        max_products_tokens=200,
        max_absolute_formal_charge=2,
    )


@pytest.fixture
def good_reaction():
    return Reaction(
        "O=[N+]([O-])c1cc(-c2nc3ccccc3o2)ccc1F.C~C>>Nc1cc(-c2nc3ccccc3o2)ccc1NCC(=O)N1CCOCC1"
    )


@pytest.fixture
def bad_reaction():
    return Reaction(
        "[C].C.[O--].[O--].O.O=[N+]([O-])c1cc(-c2nc3ccccc3o2)ccc1F.C.C>O>O.C"
    )


@pytest.fixture
def big_reaction():
    return Reaction(
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC>CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC>CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
    )


@pytest.fixture
def alchemic_reaction():
    return Reaction("C>[Hg]>[Au]")


@pytest.fixture
def alchemic_reaction():
    return Reaction("C>[Hg]>[Au]")


def test_max_reactants_exceeded(filter, good_reaction, bad_reaction):
    assert not filter.max_reactants_exceeded(good_reaction)
    assert filter.max_reactants_exceeded(bad_reaction)


def test_max_agents_exceeded(filter, good_reaction, bad_reaction):
    assert not filter.max_agents_exceeded(good_reaction)
    assert filter.max_agents_exceeded(bad_reaction)


def test_max_products_exceeded(filter, good_reaction, bad_reaction):
    assert not filter.max_products_exceeded(good_reaction)
    assert filter.max_products_exceeded(bad_reaction)


def test_products_subset_of_reactants(filter, good_reaction, bad_reaction):
    assert not filter.products_subset_of_reactants(good_reaction)
    assert filter.products_subset_of_reactants(bad_reaction)


def test_products_single_atoms(filter, good_reaction, bad_reaction):
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