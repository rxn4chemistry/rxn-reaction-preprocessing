import pandas as pd
import pytest

import rxn_reaction_preprocessing as rrp
from rxn_reaction_preprocessing import Augmenter
from rxn_reaction_preprocessing.utils import reset_random_seed


@pytest.fixture
def augmenter():
    reset_random_seed()
    df = pd.DataFrame(
        {
            'smiles':
                [
                    '[Na]Cl~[K+].[K+].CCCC.CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
                    'CC(C)(C)O[K]~CCC.O',
                ],
        }
    )
    return Augmenter(df, fragment_bond='~')


def test_restricted(augmenter):
    new_df = augmenter.augment(rrp.RandomType['restricted'], permutations=1, tokenize=False)

    expected = [
        'CCCC.Cl[Na]~[K+].Cn1c2c(c(=O)n(C)c1=O)n(C)cn2.[K+]',
        'CC(C)(C)O[K]~CCC.O',
    ]
    assert all([new_df['restricted'].values[i] == expected[i] for i in range(len(expected))])


def test_unrestricted(augmenter):
    new_df = augmenter.augment(rrp.RandomType['unrestricted'], permutations=1, tokenize=False)

    expected = [
        'C(C)CC.[K+].[Na]Cl~[K+].c12ncn(c1c(n(c(n2C)=O)C)=O)C',
        'C(C)(C)(C)O[K]~CCC.O',
    ]
    assert all([new_df['unrestricted'].values[i] == expected[i] for i in range(len(expected))])


def test_rotated(augmenter):
    new_df = augmenter.augment(rrp.RandomType['rotated'], permutations=1, tokenize=False)

    expected = [
        'C(CC)C.[K+].[Na]Cl~[K+].n1c2c(c(=O)n(C)c(=O)n2C)n(C)c1',
        'C(C)(C)(O[K])C~CCC.O',
    ]
    assert all([new_df['rotated'].values[i] == expected[i] for i in range(len(expected))])


def test_molecules_order(augmenter):
    new_df = augmenter.augment(rrp.RandomType['molecules'], permutations=1, tokenize=False)

    expected = [
        '[Na]Cl~[K+].CN1C=NC2=C1C(=O)N(C(=O)N2C)C.[K+].CCCC',
        'CC(C)(C)O[K]~CCC.O',
    ]
    assert all([new_df['molecules'].values[i] == expected[i] for i in range(len(expected))])


def test_multiple_augmentation_molecules(augmenter):
    new_df = augmenter.augment(rrp.RandomType['molecules'], permutations=3, tokenize=False)

    expected = [
        '[Na]Cl~[K+].CN1C=NC2=C1C(=O)N(C(=O)N2C)C.[K+].CCCC',
        'CN1C=NC2=C1C(=O)N(C(=O)N2C)C.[K+].CCCC.[Na]Cl~[K+]',
        '[Na]Cl~[K+].[K+].CCCC.CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'O.CC(C)(C)O[K]~CCC',
        'CC(C)(C)O[K]~CCC.O'
    ]
    assert all([new_df['molecules'].values[i] == expected[i] for i in range(len(expected))])


def test_multiple_augmentation_rotated(augmenter):
    new_df = augmenter.augment(rrp.RandomType['rotated'], permutations=3, tokenize=False)

    expected = [
        'C(CC)C.[K+].[Na]Cl~[K+].n1c2c(c(=O)n(C)c(=O)n2C)n(C)c1',
        'CCCC.Cn1cnc2c1c(=O)n(C)c(=O)n2C.[K+].[Na]Cl~[K+]',
        'C(CC)C.[K+].[Na]Cl~[K+].n1(C)c(=O)n(C)c2ncn(C)c2c1=O', 'O.O([K])C(C)(C)C~CCC',
        'O.[K]OC(C)(C)C~CCC', 'C(C)(C)(O[K])C~C(C)C.O'
    ]
    assert all([new_df['rotated'].values[i] == expected[i] for i in range(len(expected))])
