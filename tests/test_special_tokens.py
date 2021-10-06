from rxn_chemutils.conversion import canonicalize_smiles
from rxn_chemutils.reaction_equation import ReactionEquation

from rxn_reaction_preprocessing.special_tokens import add_heat_token
from rxn_reaction_preprocessing.special_tokens import add_light_token
from rxn_reaction_preprocessing.special_tokens import contains_heat_token
from rxn_reaction_preprocessing.special_tokens import contains_light_token
from rxn_reaction_preprocessing.special_tokens import HEAT_TOKEN
from rxn_reaction_preprocessing.special_tokens import LIGHT_TOKEN
from rxn_reaction_preprocessing.special_tokens import strip_all_special_tokens
from rxn_reaction_preprocessing.special_tokens import strip_heat_token
from rxn_reaction_preprocessing.special_tokens import strip_light_token


def test_special_tokens_are_canonical():
    assert canonicalize_smiles(LIGHT_TOKEN) == LIGHT_TOKEN
    assert canonicalize_smiles(HEAT_TOKEN) == HEAT_TOKEN


def test_add_light_token():
    reaction = ReactionEquation.from_string('A.B>>C')
    add_light_token(reaction)
    assert reaction.to_string() == f'A.B.{LIGHT_TOKEN}>>C'


def test_add_heat_token():
    reaction = ReactionEquation.from_string('A.B>>C')
    add_heat_token(reaction)
    assert reaction.to_string() == f'A.B.{HEAT_TOKEN}>>C'


def test_contains_light_token():
    # no light token
    reaction = ReactionEquation.from_string('A.B>>C')
    assert not contains_light_token(reaction)
    reaction = ReactionEquation.from_string(f'A.B.{HEAT_TOKEN}>>C')
    assert not contains_light_token(reaction)

    # Match independently of the location in the reactants
    reaction = ReactionEquation.from_string(f'A.B.{LIGHT_TOKEN}>>C')
    assert contains_light_token(reaction)
    reaction = ReactionEquation.from_string(f'{LIGHT_TOKEN}.A.B>>C')
    assert contains_light_token(reaction)

    # Match if in the agents or products instead of the reactants
    reaction = ReactionEquation.from_string(f'A.B>{LIGHT_TOKEN}>C')
    assert contains_light_token(reaction)
    reaction = ReactionEquation.from_string(f'A.B>>C.{LIGHT_TOKEN}')
    assert contains_light_token(reaction)

    # For the highly unlikely case that the token is part of a bigger molecule:
    # should not be considered to be present
    reaction = ReactionEquation.from_string(f'A.B.O{LIGHT_TOKEN}>>C')
    assert not contains_light_token(reaction)


def test_contains_heat_token():
    # No heat token
    reaction = ReactionEquation.from_string('A.B>>C')
    assert not contains_heat_token(reaction)
    reaction = ReactionEquation.from_string(f'A.B.{LIGHT_TOKEN}>>C')
    assert not contains_heat_token(reaction)

    # Match independently of the location in the reactants
    reaction = ReactionEquation.from_string(f'A.B.{HEAT_TOKEN}>>C')
    assert contains_heat_token(reaction)
    reaction = ReactionEquation.from_string(f'{HEAT_TOKEN}.A.B>>C')
    assert contains_heat_token(reaction)

    # Match if in the agents or products instead of the reactants
    reaction = ReactionEquation.from_string(f'A.B>{HEAT_TOKEN}>C')
    assert contains_heat_token(reaction)
    reaction = ReactionEquation.from_string(f'A.B>>C.{HEAT_TOKEN}')
    assert contains_heat_token(reaction)

    # For the highly unlikely case that the token is part of a bigger molecule:
    # should not be considered to be present
    reaction = ReactionEquation.from_string(f'A.B.O{HEAT_TOKEN}>>C')
    assert not contains_heat_token(reaction)


def test_strip_all_special_tokens():
    # No special token - no change needed
    reaction = ReactionEquation.from_string('A.B>>C')
    strip_all_special_tokens(reaction)
    assert reaction.to_string() == 'A.B>>C'

    reaction = ReactionEquation.from_string(f'A.B>{LIGHT_TOKEN}>C')
    strip_all_special_tokens(reaction)
    assert reaction.to_string() == 'A.B>>C'

    reaction = ReactionEquation.from_string(f'A.B>{LIGHT_TOKEN}>C.{HEAT_TOKEN}')
    strip_all_special_tokens(reaction)
    assert reaction.to_string() == 'A.B>>C'

    reaction = ReactionEquation.from_string(f'{LIGHT_TOKEN}.A.B>{LIGHT_TOKEN}>C.{HEAT_TOKEN}')
    strip_all_special_tokens(reaction)
    assert reaction.to_string() == 'A.B>>C'

    # For the highly unlikely case that the token is part of a bigger molecule:
    # will not be removed
    reaction = ReactionEquation.from_string(f'A.B.O{HEAT_TOKEN}>>C')
    strip_all_special_tokens(reaction)
    assert reaction.to_string() == f'A.B.O{HEAT_TOKEN}>>C'


def test_strip_light_token():
    # No special token - no change needed
    reaction = ReactionEquation.from_string('A.B>>C')
    strip_light_token(reaction)
    assert reaction.to_string() == 'A.B>>C'

    reaction = ReactionEquation.from_string(f'A.B>{LIGHT_TOKEN}>C')
    strip_light_token(reaction)
    assert reaction.to_string() == 'A.B>>C'

    reaction = ReactionEquation.from_string(f'A.B>{LIGHT_TOKEN}>C.{LIGHT_TOKEN}')
    strip_light_token(reaction)
    assert reaction.to_string() == 'A.B>>C'

    # Does not remove the heat token
    reaction = ReactionEquation.from_string(f'{LIGHT_TOKEN}.A.B.{HEAT_TOKEN}>>C')
    strip_light_token(reaction)
    assert reaction.to_string() == f'A.B.{HEAT_TOKEN}>>C'


def test_strip_heat_token():
    # No special token - no change needed
    reaction = ReactionEquation.from_string('A.B>>C')
    strip_heat_token(reaction)
    assert reaction.to_string() == 'A.B>>C'

    reaction = ReactionEquation.from_string(f'A.B>{HEAT_TOKEN}>C')
    strip_heat_token(reaction)
    assert reaction.to_string() == 'A.B>>C'

    reaction = ReactionEquation.from_string(f'A.B>{HEAT_TOKEN}>C.{HEAT_TOKEN}')
    strip_heat_token(reaction)
    assert reaction.to_string() == 'A.B>>C'

    # Does not remove the light token
    reaction = ReactionEquation.from_string(f'{LIGHT_TOKEN}.A.B.{HEAT_TOKEN}>>C')
    strip_heat_token(reaction)
    assert reaction.to_string() == f'{LIGHT_TOKEN}.A.B>>C'
