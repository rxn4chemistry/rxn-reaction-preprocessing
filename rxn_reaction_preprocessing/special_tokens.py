"""
Module for functionality related to special tokens (such as light and heat) to
include in reaction SMILES.

Important comments:
1. Other modules should always try to use the functions defined here, i.e. try to:
   - never explicitly mention "[Lv]" or "[Ts]"
   - not refer directly to `LIGHT_TOKEN` or `HEAT_TOKEN`: normally the
     functions in this module should be sufficient (adding such tokens,
     removing them, querying whether a reaction contains them).
2. The functions here rely on ReactionEquation instead of SMILES strings. This
   makes the code independent of the reaction SMILES format. Conversion to and
   from strings should be done outside this module.
3. The Python objects starting with an underscore are meant not to be used elsewhere.
"""
from enum import Enum
from typing import Iterable

from rxn_chemutils.reaction_equation import ReactionEquation

LIGHT_TOKEN = '[Lv]'
HEAT_TOKEN = '[Ts]'


class _SpecialToken(Enum):
    """
    Enum class for special reaction SMILES tokens.

    Useful to avoid dealing with the token strings where not actually necessary.
    """
    LIGHT = LIGHT_TOKEN
    HEAT = HEAT_TOKEN


def _add_special_tokens(reaction: ReactionEquation, tokens: Iterable[_SpecialToken]) -> None:
    """Add the required tokens to the reactants of a reaction (in-place)."""
    for token in tokens:
        reaction.reactants.append(token.value)


def add_light_token(reaction: ReactionEquation) -> None:
    """Add the light token to the precursors of a reaction (in-place)."""
    _add_special_tokens(reaction, [_SpecialToken.LIGHT])


def add_heat_token(reaction: ReactionEquation) -> None:
    """Add the heat token to the precursors of a reaction (in-place)."""
    _add_special_tokens(reaction, [_SpecialToken.HEAT])


def _contains_token(reaction: ReactionEquation, token: _SpecialToken) -> bool:
    """Whether a reaction contains the specified token."""
    return any(compound == token.value for compound in reaction.iter_all_smiles())


def contains_light_token(reaction: ReactionEquation) -> bool:
    """Whether a reaction contains the light token."""
    return _contains_token(reaction, _SpecialToken.LIGHT)


def contains_heat_token(reaction: ReactionEquation) -> bool:
    """Whether a reaction contains the heat token."""
    return _contains_token(reaction, _SpecialToken.HEAT)


def _strip_special_tokens(reaction: ReactionEquation, tokens: Iterable[_SpecialToken]) -> None:
    """Strip the specified tokens from a reaction (in-place)."""
    token_strings_to_remove = [token.value for token in tokens]

    for reaction_group in reaction:
        for token in token_strings_to_remove:
            try:
                reaction_group.remove(token)
            except ValueError:
                # NB: remove() raises ValueError if the value is not in the list
                pass


def strip_all_special_tokens(reaction: ReactionEquation) -> None:
    """Strip all the special tokens from a reaction (in-place)."""
    # NB: calling list on an enum class gets all the possible values.
    _strip_special_tokens(reaction, list(_SpecialToken))


def strip_heat_token(reaction: ReactionEquation) -> None:
    """Strip the heat from a reaction (in-place)."""
    _strip_special_tokens(reaction, [_SpecialToken.HEAT])


def strip_light_token(reaction: ReactionEquation) -> None:
    """Strip the light token from a reaction (in-place)."""
    _strip_special_tokens(reaction, [_SpecialToken.LIGHT])
