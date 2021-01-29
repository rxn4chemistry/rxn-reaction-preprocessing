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
    new_df = augmenter.augment(rrp.RandomType.restricted, permutations=1, tokenize=False)

    expected = [
        'CCCC.CN1C2=C(C(=O)N(C)C1=O)N(C)C=N2.Cl[Na]~[K+].[K+]',
        'CC(C)(C)O[K]~CCC.O',
    ]
    assert new_df['restricted'].tolist() == expected


def test_unrestricted(augmenter):
    new_df = augmenter.augment(rrp.RandomType.unrestricted, permutations=1, tokenize=False)

    expected = [
        'C(C)CC.C12N=CN(C=1C(N(C(N2C)=O)C)=O)C.[K+].[Na]Cl~[K+]',
        'C(C)(C)(C)O[K]~CCC.O',
    ]
    assert new_df['unrestricted'].tolist() == expected


def test_rotated(augmenter):
    new_df = augmenter.augment(rrp.RandomType.rotated, permutations=1, tokenize=False)

    expected = [
        'C1(=O)N(C)C2=C(N(C)C=N2)C(=O)N1C.CCCC.Cl[Na]~[K+].[K+]',
        'O.[K]OC(C)(C)C~CCC',
    ]
    assert new_df['rotated'].tolist() == expected


def test_molecules_order(augmenter):
    new_df = augmenter.augment(rrp.RandomType.molecules, permutations=1, tokenize=False)

    expected = [
        '[Na]Cl~[K+].CN1C=NC2=C1C(=O)N(C(=O)N2C)C.[K+].CCCC',
        'CC(C)(C)O[K]~CCC.O',
    ]
    assert new_df['molecules'].tolist() == expected


def test_multiple_augmentation_molecules(augmenter):
    new_df = augmenter.augment(rrp.RandomType.molecules, permutations=3, tokenize=False)

    expected = [
        '[Na]Cl~[K+].CN1C=NC2=C1C(=O)N(C(=O)N2C)C.[K+].CCCC',
        'CN1C=NC2=C1C(=O)N(C(=O)N2C)C.[K+].CCCC.[Na]Cl~[K+]',
        '[Na]Cl~[K+].[K+].CCCC.CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'O.CC(C)(C)O[K]~CCC',
        'CC(C)(C)O[K]~CCC.O'
    ]
    assert new_df['molecules'].tolist() == expected


def test_multiple_augmentation_rotated(augmenter):
    new_df = augmenter.augment(rrp.RandomType.rotated, permutations=3, tokenize=False)

    expected = [
        'C1(=O)N(C)C2=C(N(C)C=N2)C(=O)N1C.CCCC.Cl[Na]~[K+].[K+]',
        'C(CC)C.Cl[Na]~[K+].N1=CN(C)C2=C1N(C)C(=O)N(C)C2=O.[K+]',
        'CCCC.CN1C2=C(N(C)C=N2)C(=O)N(C)C1=O.[K+].[Na]Cl~[K+]', 'O.[K]OC(C)(C)C~CCC',
        'O.O([K])C(C)(C)C~CCC', 'CC(C)(C)O[K]~CCC.O'
    ]
    assert new_df['rotated'].tolist() == expected
