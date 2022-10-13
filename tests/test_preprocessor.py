# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2022
# ALL RIGHTS RESERVED
import pandas as pd
import pytest

from rxn.reaction_preprocessing import MixedReactionFilter, Preprocessor


@pytest.fixture
def preprocessor() -> Preprocessor:
    df = pd.DataFrame(
        {
            "rxn": [
                "[14C]Cl.O[Na]>O>[Na]Cl.[14C]O",
                "[14C]Cl.O[Na]>>[Na]Cl",
                "[C].C.[O--].[O--].O.O=[N+]([O-])c1cc(-c2nc3ccccc3o2)ccc1F.C.C>O>O.C",
                "=)(/?",
            ],
            "class": [2, 2, 6, 0],
        }
    )
    return Preprocessor(df, "rxn")


@pytest.fixture
def filter() -> MixedReactionFilter:
    return MixedReactionFilter(
        max_reactants=5,
        max_agents=0,
        max_products=1,
        min_reactants=2,
        min_agents=0,
        min_products=1,
        max_reactants_tokens=300,
        max_agents_tokens=0,
        max_products_tokens=200,
        max_absolute_formal_charge=2,
    )


def test_filter(preprocessor: Preprocessor, filter: MixedReactionFilter) -> None:
    preprocessor.filter(filter, False)
    assert not preprocessor.df._rxn_valid.iloc[0]
    assert preprocessor.df._rxn_valid.iloc[1]
    assert not preprocessor.df._rxn_valid.iloc[2]


def test_filter_verbose(
    preprocessor: Preprocessor, filter: MixedReactionFilter
) -> None:
    preprocessor.filter(filter, True)
    assert len(preprocessor.df._rxn_valid_messages.iloc[0]) == 3
    assert len(preprocessor.df._rxn_valid_messages.iloc[1]) == 0
    assert len(preprocessor.df._rxn_valid_messages.iloc[2]) == 7


def test_preprocess_with_rdkit_bug_for_concatenated_smiles() -> None:
    # In an earlier version, preprocessing failed for rare cases when RDKit cannot
    # convert concatenated SMILES into an RDKit Mol.
    # Was related to https://github.com/rdkit/rdkit/issues/4266

    df = pd.DataFrame({"rxn": ["c1ccccc1.C123C45C16C21C34C561>>CC"]})
    preprocessor = Preprocessor(df, "rxn")
    mrf = MixedReactionFilter()

    # This used to raise an exception
    preprocessor.filter(mrf)
