from rxn.chemutils.conversion import canonicalize_smiles
from rxn.chemutils.reaction_equation import ReactionEquation

from rxn.reaction_preprocessing.special_tokens import (
    HEAT_TOKEN,
    LIGHT_TOKEN,
    add_heat_token,
    add_light_token,
    contains_heat_token,
    contains_light_token,
    strip_all_special_tokens,
    strip_heat_token,
    strip_light_token,
)


def test_special_tokens_are_canonical() -> None:
    assert canonicalize_smiles(LIGHT_TOKEN) == LIGHT_TOKEN
    assert canonicalize_smiles(HEAT_TOKEN) == HEAT_TOKEN


def test_add_light_token() -> None:
    # in-place
    reaction = ReactionEquation.from_string("A.B>>C")
    updated = add_light_token(reaction, in_place=True)
    assert updated is reaction
    assert updated.to_string() == f"A.B.{LIGHT_TOKEN}>>C"

    # not in-place
    reaction = ReactionEquation.from_string("A.B>>C")
    updated = add_light_token(reaction, in_place=False)
    assert updated is not reaction
    assert updated.to_string() == f"A.B.{LIGHT_TOKEN}>>C"
    assert reaction.to_string() == "A.B>>C"


def test_add_heat_token() -> None:
    # in-place
    reaction = ReactionEquation.from_string("A.B>>C")
    updated = add_heat_token(reaction, in_place=True)
    assert updated is reaction
    assert updated.to_string() == f"A.B.{HEAT_TOKEN}>>C"

    # not in-place
    reaction = ReactionEquation.from_string("A.B>>C")
    updated = add_heat_token(reaction, in_place=False)
    assert updated is not reaction
    assert updated.to_string() == f"A.B.{HEAT_TOKEN}>>C"
    assert reaction.to_string() == "A.B>>C"


def test_contains_light_token() -> None:
    # no light token
    reaction = ReactionEquation.from_string("A.B>>C")
    assert not contains_light_token(reaction)
    reaction = ReactionEquation.from_string(f"A.B.{HEAT_TOKEN}>>C")
    assert not contains_light_token(reaction)

    # Match independently of the location in the reactants
    reaction = ReactionEquation.from_string(f"A.B.{LIGHT_TOKEN}>>C")
    assert contains_light_token(reaction)
    reaction = ReactionEquation.from_string(f"{LIGHT_TOKEN}.A.B>>C")
    assert contains_light_token(reaction)

    # Match if in the agents or products instead of the reactants
    reaction = ReactionEquation.from_string(f"A.B>{LIGHT_TOKEN}>C")
    assert contains_light_token(reaction)
    reaction = ReactionEquation.from_string(f"A.B>>C.{LIGHT_TOKEN}")
    assert contains_light_token(reaction)

    # For the highly unlikely case that the token is part of a bigger molecule:
    # should not be considered to be present
    reaction = ReactionEquation.from_string(f"A.B.O{LIGHT_TOKEN}>>C")
    assert not contains_light_token(reaction)


def test_contains_heat_token() -> None:
    # No heat token
    reaction = ReactionEquation.from_string("A.B>>C")
    assert not contains_heat_token(reaction)
    reaction = ReactionEquation.from_string(f"A.B.{LIGHT_TOKEN}>>C")
    assert not contains_heat_token(reaction)

    # Match independently of the location in the reactants
    reaction = ReactionEquation.from_string(f"A.B.{HEAT_TOKEN}>>C")
    assert contains_heat_token(reaction)
    reaction = ReactionEquation.from_string(f"{HEAT_TOKEN}.A.B>>C")
    assert contains_heat_token(reaction)

    # Match if in the agents or products instead of the reactants
    reaction = ReactionEquation.from_string(f"A.B>{HEAT_TOKEN}>C")
    assert contains_heat_token(reaction)
    reaction = ReactionEquation.from_string(f"A.B>>C.{HEAT_TOKEN}")
    assert contains_heat_token(reaction)

    # For the highly unlikely case that the token is part of a bigger molecule:
    # should not be considered to be present
    reaction = ReactionEquation.from_string(f"A.B.O{HEAT_TOKEN}>>C")
    assert not contains_heat_token(reaction)


def test_strip_all_special_tokens() -> None:
    # No special token - no change needed
    reaction = ReactionEquation.from_string("A.B>>C")
    assert strip_all_special_tokens(reaction).to_string() == "A.B>>C"

    reaction = ReactionEquation.from_string(f"A.B>{LIGHT_TOKEN}>C")
    assert strip_all_special_tokens(reaction).to_string() == "A.B>>C"

    reaction = ReactionEquation.from_string(f"A.B>{LIGHT_TOKEN}>C.{HEAT_TOKEN}")
    assert strip_all_special_tokens(reaction).to_string() == "A.B>>C"

    reaction = ReactionEquation.from_string(
        f"{LIGHT_TOKEN}.A.B>{LIGHT_TOKEN}>C.{HEAT_TOKEN}"
    )
    assert strip_all_special_tokens(reaction).to_string() == "A.B>>C"

    # For the highly unlikely case that the token is part of a bigger molecule:
    # will not be removed
    reaction = ReactionEquation.from_string(f"A.B.O{HEAT_TOKEN}>>C")
    assert strip_all_special_tokens(reaction).to_string() == f"A.B.O{HEAT_TOKEN}>>C"

    # in-place
    reaction = ReactionEquation.from_string(f"A.B>{LIGHT_TOKEN}>C")
    updated = strip_all_special_tokens(reaction, in_place=True)
    assert updated is reaction

    # not in-place
    reaction = ReactionEquation.from_string(f"A.B>{LIGHT_TOKEN}>C")
    updated = strip_all_special_tokens(reaction, in_place=False)
    assert updated is not reaction
    assert updated != reaction


def test_strip_light_token() -> None:
    # No special token - no change needed
    reaction = ReactionEquation.from_string("A.B>>C")
    assert strip_light_token(reaction).to_string() == "A.B>>C"

    reaction = ReactionEquation.from_string(f"A.B>{LIGHT_TOKEN}>C")
    assert strip_light_token(reaction).to_string() == "A.B>>C"

    reaction = ReactionEquation.from_string(f"A.B>{LIGHT_TOKEN}>C.{LIGHT_TOKEN}")
    assert strip_light_token(reaction).to_string() == "A.B>>C"

    # Does not remove the heat token
    reaction = ReactionEquation.from_string(f"{LIGHT_TOKEN}.A.B.{HEAT_TOKEN}>>C")
    assert strip_light_token(reaction).to_string() == f"A.B.{HEAT_TOKEN}>>C"

    # in-place
    reaction = ReactionEquation.from_string(f"A.B>{LIGHT_TOKEN}>C")
    updated = strip_light_token(reaction, in_place=True)
    assert updated is reaction

    # not in-place
    reaction = ReactionEquation.from_string(f"A.B>{LIGHT_TOKEN}>C")
    updated = strip_light_token(reaction, in_place=False)
    assert updated is not reaction
    assert updated != reaction


def test_strip_heat_token() -> None:
    # No special token - no change needed
    reaction = ReactionEquation.from_string("A.B>>C")
    assert strip_heat_token(reaction).to_string() == "A.B>>C"

    reaction = ReactionEquation.from_string(f"A.B>{HEAT_TOKEN}>C")
    assert strip_heat_token(reaction).to_string() == "A.B>>C"

    reaction = ReactionEquation.from_string(f"A.B>{HEAT_TOKEN}>C.{HEAT_TOKEN}")
    assert strip_heat_token(reaction).to_string() == "A.B>>C"

    # Does not remove the light token
    reaction = ReactionEquation.from_string(f"{LIGHT_TOKEN}.A.B.{HEAT_TOKEN}>>C")
    assert strip_heat_token(reaction).to_string() == f"{LIGHT_TOKEN}.A.B>>C"

    # in-place
    reaction = ReactionEquation.from_string(f"A.B>{HEAT_TOKEN}>C")
    updated = strip_heat_token(reaction, in_place=True)
    assert updated is reaction

    # not in-place
    reaction = ReactionEquation.from_string(f"A.B>{HEAT_TOKEN}>C")
    updated = strip_heat_token(reaction, in_place=False)
    assert updated is not reaction
    assert updated != reaction


def test_add_on_lists() -> None:
    # The add_* functions work not only on ReactionEquation instances,
    # but also on lists of strings

    # in-place
    smiles = ["A", "B", "C"]
    updated = add_light_token(smiles, in_place=True)
    assert updated is smiles
    assert updated == ["A", "B", "C", LIGHT_TOKEN]
    add_light_token(smiles, in_place=True)
    add_heat_token(smiles, in_place=True)
    assert updated == ["A", "B", "C", LIGHT_TOKEN, LIGHT_TOKEN, HEAT_TOKEN]

    # not in-place
    smiles = ["A", "B", "C"]
    updated = add_light_token(smiles, in_place=False)
    assert updated is not smiles
    assert updated == ["A", "B", "C", LIGHT_TOKEN]
    assert smiles == ["A", "B", "C"]


def test_contains_on_iterables() -> None:
    # The contains_* functions work not only on ReactionEquation instances,
    # but also on iterables of strings

    # on lists
    assert not contains_heat_token(["A", "B", "C.D"])
    assert contains_heat_token(["A", "B", HEAT_TOKEN, "C.D"])
    assert not contains_heat_token(["A", f"B.{HEAT_TOKEN}", "C.D"])

    # on sets
    assert contains_heat_token({"A", "B", HEAT_TOKEN, "C.D"})
    assert not contains_light_token({"A", "B", HEAT_TOKEN, "C.D"})

    # on tuples
    assert contains_heat_token(("A", "B", HEAT_TOKEN, "C.D"))
    assert not contains_light_token(("A", "B", HEAT_TOKEN, "C.D"))
    assert contains_light_token(("A", "B", LIGHT_TOKEN, "C.D"))


def test_strip_on_list() -> None:
    # The strip_* functions work not only on ReactionEquation instances,
    # but also on lists of strings

    smiles_list = ["A", "B"]
    assert strip_light_token(smiles_list) == ["A", "B"]

    smiles_list = ["A", "B", LIGHT_TOKEN, "D"]
    assert strip_heat_token(smiles_list) == ["A", "B", LIGHT_TOKEN, "D"]
    assert strip_light_token(smiles_list) == ["A", "B", "D"]

    smiles_list = [HEAT_TOKEN, "A", "B", LIGHT_TOKEN, "D"]
    assert strip_all_special_tokens(smiles_list) == ["A", "B", "D"]

    # in-place
    smiles_list = [HEAT_TOKEN, "A", "B", LIGHT_TOKEN, "D"]
    updated = strip_all_special_tokens(smiles_list, in_place=True)
    assert updated is smiles_list

    # not in-place
    smiles_list = [HEAT_TOKEN, "A", "B", LIGHT_TOKEN, "D"]
    updated = strip_all_special_tokens(smiles_list, in_place=False)
    assert updated is not smiles_list
    assert updated != smiles_list
