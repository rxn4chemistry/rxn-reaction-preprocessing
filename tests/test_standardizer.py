import numpy as np
import pandas as pd
import pytest

from rxn_reaction_preprocessing import Patterns
from rxn_reaction_preprocessing import Standardizer


@pytest.fixture
def standardizer():
    df = pd.DataFrame(
        {
            'rxn':
                [
                    '[Li]O>>[Na]Cl',
                    'CC(C)(C)O[K]~CCC>>[Li]O',
                    r'CC(=O)/C=C(\C)O[V](=O)O/C(C)=C/C(C)=O~O=C(O[K])O[K].[Li]O>>O[K]',
                ],
        }
    )
    pt = Patterns(jsonpath='tests/test_patterns', fragment_bond='~')
    return Standardizer(df, pt, 'rxn', fragment_bond=None)


@pytest.fixture
def standardizer_with_same_fragment():
    df = pd.DataFrame(
        {
            'rxn':
                [
                    '[Li]O.O[Na]>>[Na]Cl', 'CC(C)(C)O[K].CCC>>[Li]O',
                    r'CC(=O)/C=C(\C)O[V](=O)O/C(C)=C/C(C)=O.O=C(O[K])O[K].[Li]O>>O[K]',
                    'CO.COC(=O)C(C)(C)c1ccc(C(=O)CCCN2CCC(C(O)(c3ccccc3)c3ccccc3)CC2)cc1.Cl.O[Na]>>CC(C)(C(=O)O)c1ccc(C(=O)CCCN2CCC(C(O)(c3ccccc3)c3ccccc3)CC2)cc1~Cl'
                ],
        }
    )
    pt = Patterns(jsonpath='tests/test_patterns', fragment_bond='~')
    return Standardizer(df, pt, 'rxn', fragment_bond='~')


@pytest.fixture
def standardizer_with_different_fragment():
    df = pd.DataFrame(
        {
            'rxn':
                [
                    '[Li]O>>[Na]Cl',
                    'CC(C)(C)O[K]$CCC>>[Li]O',
                    r'CC(=O)/C=C(\C)O[V](=O)O/C(C)=C/C(C)=O$O=C(O[K])O[K].[Li]O>>O[K]',
                ],
        }
    )
    pt = Patterns(jsonpath='tests/test_patterns', fragment_bond='~')
    return Standardizer(df, pt, 'rxn', fragment_bond='$')


def test_patterns():
    pt = Patterns(jsonpath='tests/test_patterns', fragment_bond='~')

    assert type(pt.patterns) is dict
    assert list(pt.patterns.keys()) == ['exception_patterns', 'potassium_compounds']


def test_patterns_manipulation():
    modified = [
        [
            r'(?:^|(?<=\~|\.|>))O\[K\](?=\~|\.|>|$)',
            r'(?:^|(?<=\~|\.|>))\[Li\]O(?=\~|\.|>|$)',
            r'(?:^|(?<=\~|\.|>))O\[Na\](?=\~|\.|>|$)',
        ],
        [
            r'(?:^|(?<=\~|\.|>))O=C\(O\[K\]\)O\[K\](?=\~|\.|>|$)',
            r'(?:^|(?<=\~|\.|>))CC\(C\)\(C\)O\[K\](?=\~|\.|>|$)',
            '(?:^|(?<=\\~|\\.|>))CC\\(=O\\)/C=C\\(\\\\C\\)O\\[V\\]\\(=O\\)O/C\\(C\\)=C/C\\(C\\)=O(?=\\~|\\.|>|$)'
        ]
    ]

    pt = Patterns(jsonpath='tests/test_patterns', fragment_bond='~')

    for i, elem in enumerate(pt.patterns.values()):
        for j, pat in enumerate(elem):
            assert pat[0] == modified[i][j]


def test_standardization(standardizer):
    pass
    # Not working
    # new_df = standardizer.standardize().df

    # converted_rxns = [
    #     "[Li+]~[OH-]>>[Na]Cl",
    #     "CC(C)(C)[O-]~[K+]~CCC>>[Li+]~[OH-]",
    #     "CC(=O)C=C(C)[O-]~CC(=O)C=C(C)[O-]~O=[V+2]~O=C([O-])[O-]~[K+]~[K+].[Li+]~[OH-]>>[K+]~[OH-]",
    # ]

    # assert all(
    #     [new_df['_rxn_std'].values[i] == converted_rxns[i] for i in range(len(converted_rxns))]
    # )


def test_standardization_with_same_fragment(standardizer_with_same_fragment):
    pass
    # Not working
    # new_df = standardizer_with_same_fragment.standardize().df

    # converted_rxns = [
    #     "[Li+]~[OH-].[Na+]~[OH-]>>[Na]Cl", "CC(C)(C)[O-]~[K+].CCC>>[Li+]~[OH-]",
    #     "CC(=O)C=C(C)[O-]~CC(=O)C=C(C)[O-]~O=[V+2].O=C([O-])[O-]~[K+]~[K+].[Li+]~[OH-]>>[K+]~[OH-]",
    #     "CO.COC(=O)C(C)(C)c1ccc(C(=O)CCCN2CCC(C(O)(c3ccccc3)c3ccccc3)CC2)cc1.Cl.[Na+]~[OH-]>>CC(C)(C(=O)O)c1ccc(C(=O)CCCN2CCC(C(O)(c3ccccc3)c3ccccc3)CC2)cc1~Cl"
    # ]

    # assert all(
    #     [new_df['_rxn_std'].values[i] == converted_rxns[i] for i in range(len(converted_rxns))]
    # )


def test_standardization_with_different_fragment(standardizer_with_different_fragment):
    pass
    # Not working
    # new_df = standardizer_with_different_fragment.standardize().df

    # converted_rxns = [
    #     "[Li+]$[OH-]>>[Na]Cl",
    #     "CC(C)(C)[O-]$[K+]$CCC>>[Li+]$[OH-]",
    #     "CC(=O)C=C(C)[O-]$CC(=O)C=C(C)[O-]$O=[V+2]$O=C([O-])[O-]$[K+]$[K+].[Li+]$[OH-]>>[K+]$[OH-]",
    # ]

    # assert all(
    #     [new_df['_rxn_std'].values[i] == converted_rxns[i] for i in range(len(converted_rxns))]
    # )