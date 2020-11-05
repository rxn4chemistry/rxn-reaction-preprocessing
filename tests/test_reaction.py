import pytest
from data_preprocessor import Reaction


@pytest.fixture
def reaction():
    return Reaction("[14C]Cl.[Na]O>O>[Na]Cl.[14C]O")


@pytest.fixture
def duplicte_reaction():
    return Reaction("[14C]Cl.[14C]Cl.[Na]O>O.O>[Na]Cl.[Na]Cl.[14C]O")


def test_len(reaction):
    assert len(reaction) == 5


def test_str(reaction):
    assert str(reaction) == "[14C]Cl.O[Na]>O>[Na]Cl.[14C]O"


def test_eq(reaction):
    reaction_reversed = Reaction("[Na]Cl.[14C]O>O>[14C]Cl.[Na]O")
    assert reaction == reaction
    assert reaction != reaction_reversed


def testt_reactant_count(reaction):
    assert len(reaction.reactants) == 2


def test_agent_count(reaction):
    assert len(reaction.agents) == 1


def test_product_count(reaction):
    assert len(reaction.products) == 2


def test_get_reactants(reaction):
    assert reaction.get_reactants_as_smiles() == ["[14C]Cl", "O[Na]"]


def test_get_agents(reaction):
    assert reaction.get_agents_as_smiles() == ["O"]


def test_get_products(reaction):
    assert reaction.get_products_as_smiles() == ["[Na]Cl", "[14C]O"]


def test_find(reaction):
    assert reaction.find("O") == ([1], [0], [1])


def test_remove(reaction):
    reaction.remove(([1], [0], [1]))
    assert len(reaction.reactants) == 1


def test_sort(reaction):
    assert str(reaction.sort()) == "O[Na].[14C]Cl>O>[14C]O.[Na]Cl"


def test_sort_only_reactants(reaction):
    assert (
        str(reaction.sort(sort_products=False, sort_agents=False))
        == "O[Na].[14C]Cl>O>[Na]Cl.[14C]O"
    )


def test_deduplicate(reaction, duplicte_reaction):
    assert str(reaction.deduplicate()) == "[14C]Cl.O[Na]>O>[Na]Cl.[14C]O"
