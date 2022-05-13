# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
import pytest

from rxn_reaction_preprocessing import Reaction, ReactionPart


@pytest.fixture
def reaction():
    return Reaction("[14C]Cl.[Na]O>O>[Na]Cl.[14C]O")


@pytest.fixture
def duplicate_reaction():
    return Reaction(
        "[14C]Cl.[14C]Cl.[Na]O>O.O>[Na]Cl.[Na]Cl.[14C]O", remove_duplicates=True
    )


@pytest.fixture
def dirty_reaction():
    return Reaction("[14C]Cl.[Na]O>O>[Na]Cl.[14C]O.O.[14C]Cl")


def test_remove_duplicates(reaction, duplicate_reaction):
    assert reaction == duplicate_reaction


def test_len(reaction):
    assert len(reaction) == 5


def test_str(reaction):
    assert str(reaction) == "[14C]Cl.O[Na]>O>[Na]Cl.[14C]O"


def test_eq(reaction):
    reaction_reversed = Reaction("[Na]Cl.[14C]O>O>[14C]Cl.[Na]O")
    assert reaction == reaction
    assert reaction != reaction_reversed


def test_reactant_count(reaction):
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


def test_find_in(reaction):
    assert reaction.find_in("O", ReactionPart.reactants) == [1]
    assert reaction.find_in("O", ReactionPart.agents) == [0]
    assert reaction.find_in("O", ReactionPart.products) == [1]


def test_remove(reaction):
    reaction.remove(([1], [0], [1]))
    assert len(reaction.reactants) == 1


def test_filter(reaction):
    reaction.filter(([1], [0], [1]))
    assert len(reaction.reactants) == 1


def test_sort(reaction):
    assert str(reaction.sort()) == "O[Na].[14C]Cl>O>[14C]O.[Na]Cl"


def test_sort_only_reactants(reaction):
    assert (
        str(reaction.sort(sort_products=False, sort_agents=False))
        == "O[Na].[14C]Cl>O>[Na]Cl.[14C]O"
    )


def test_remove_precursors_from_products(dirty_reaction):
    assert (
        str(dirty_reaction.remove_precursors_from_products())
        == "[14C]Cl.O[Na]>O>[Na]Cl.[14C]O"
    )


def test_has_none(reaction):
    assert not reaction.has_none()
    reaction.products.append(None)
    assert reaction.has_none()


def test_remove_none(reaction):
    reaction.products.append(None)
    assert len(reaction.remove_none().products) == 2
