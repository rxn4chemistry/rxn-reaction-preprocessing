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

import copy
from enum import Enum
from typing import Iterable, List, TypeVar, Union

from rxn.chemutils.reaction_equation import ReactionEquation

LIGHT_TOKEN = "[Lv]"
HEAT_TOKEN = "[Ts]"


class _SpecialToken(Enum):
    """
    Enum class for special reaction SMILES tokens.

    Useful to avoid dealing with the token strings where not actually necessary.
    """

    LIGHT = LIGHT_TOKEN
    HEAT = HEAT_TOKEN


ReactionOrList = TypeVar("ReactionOrList", ReactionEquation, List[str])
ReactionOrIterable = Union[ReactionEquation, Iterable[str]]


def _add_special_tokens_to_list(
    smiles_list: List[str], tokens: Iterable[_SpecialToken], in_place: bool
) -> List[str]:
    """Add the required tokens to a list of SMILES."""
    if not in_place:
        smiles_list = copy.deepcopy(smiles_list)

    for token in tokens:
        smiles_list.append(token.value)

    return smiles_list


def _add_special_tokens(
    reaction_or_list: ReactionOrList, tokens: Iterable[_SpecialToken], in_place: bool
) -> ReactionOrList:
    """Add the required tokens to the reactants of a reaction or list of SMILES."""
    if isinstance(reaction_or_list, ReactionEquation):
        # Create a copy of the ReactionEquation if not in-place - the copy can
        # then be updated in-place.
        if not in_place:
            reaction_or_list = copy.deepcopy(reaction_or_list)

        _add_special_tokens_to_list(reaction_or_list.reactants, tokens, in_place=True)

        return reaction_or_list
    else:
        # i.e., already a list
        return _add_special_tokens_to_list(reaction_or_list, tokens, in_place=in_place)


def add_light_token(
    reaction_or_list: ReactionOrList, in_place: bool = False
) -> ReactionOrList:
    """Add the light token to the precursors of a reaction or list of SMILES."""
    return _add_special_tokens(
        reaction_or_list, [_SpecialToken.LIGHT], in_place=in_place
    )


def add_heat_token(
    reaction_or_list: ReactionOrList, in_place: bool = False
) -> ReactionOrList:
    """Add the heat token to the precursors of a reaction or list of SMILES."""
    return _add_special_tokens(
        reaction_or_list, [_SpecialToken.HEAT], in_place=in_place
    )


def _contains_token(
    reaction_or_iterable: ReactionOrIterable, token: _SpecialToken
) -> bool:
    """Whether a reaction (or set of SMILES strings) contains the specified token."""
    smiles_iterable: Iterable[str]
    if isinstance(reaction_or_iterable, ReactionEquation):
        smiles_iterable = reaction_or_iterable.iter_all_smiles()
    else:
        smiles_iterable = reaction_or_iterable
    return any(compound == token.value for compound in smiles_iterable)


def contains_light_token(reaction_or_iterable: ReactionOrIterable) -> bool:
    """Whether a reaction (or set of SMILES strings) contains the light token."""
    return _contains_token(reaction_or_iterable, _SpecialToken.LIGHT)


def contains_heat_token(reaction_or_iterable: ReactionOrIterable) -> bool:
    """Whether a reaction (or set of SMILES strings) contains the heat token."""
    return _contains_token(reaction_or_iterable, _SpecialToken.HEAT)


def _strip_special_tokens_from_list(
    smiles_list: List[str], token_strings: Iterable[str], in_place: bool
) -> List[str]:
    """Strip the specified tokens from a list of SMILES strings."""

    if not in_place:
        smiles_list = copy.deepcopy(smiles_list)

    for token in token_strings:
        try:
            smiles_list.remove(token)
        except ValueError:
            # NB: remove() raises ValueError if the value is not in the list
            pass

    return smiles_list


def _strip_special_tokens(
    reaction_or_list: ReactionOrList, tokens: Iterable[_SpecialToken], in_place: bool
) -> ReactionOrList:
    """Strip the specified tokens from a reaction or list of SMILES strings."""
    token_strings_to_remove = [token.value for token in tokens]

    if isinstance(reaction_or_list, ReactionEquation):
        # Create a copy of the ReactionEquation if not in-place - the copy can
        # then be updated in-place.
        if not in_place:
            reaction_or_list = copy.deepcopy(reaction_or_list)

        for reaction_group in reaction_or_list:
            _strip_special_tokens_from_list(
                reaction_group, token_strings_to_remove, in_place=True
            )

        return reaction_or_list
    else:
        # i.e., already a list
        return _strip_special_tokens_from_list(
            reaction_or_list, token_strings_to_remove, in_place=in_place
        )


def strip_all_special_tokens(
    reaction_or_list: ReactionOrList, in_place: bool = False
) -> ReactionOrList:
    """Strip all the special tokens from a reaction or list of SMILES strings."""
    # NB: calling list on an enum class gets all the possible values.
    return _strip_special_tokens(
        reaction_or_list, list(_SpecialToken), in_place=in_place
    )


def strip_heat_token(
    reaction_or_list: ReactionOrList, in_place: bool = False
) -> ReactionOrList:
    """Strip the heat from a reaction or list of SMILES strings."""
    return _strip_special_tokens(
        reaction_or_list, [_SpecialToken.HEAT], in_place=in_place
    )


def strip_light_token(
    reaction_or_list: ReactionOrList, in_place: bool = False
) -> ReactionOrList:
    """Strip the light token from a reaction or list of SMILES strings."""
    return _strip_special_tokens(
        reaction_or_list, [_SpecialToken.LIGHT], in_place=in_place
    )
