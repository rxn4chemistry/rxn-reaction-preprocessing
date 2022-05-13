# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
import random
from string import ascii_lowercase
from typing import List, Sequence

import pandas as pd
import pytest

from rxn_reaction_preprocessing import StableDataSplitter
from rxn_reaction_preprocessing.utils import reset_random_seed


def random_strings(n: int) -> List[str]:
    reset_random_seed()
    return [
        "".join([random.choice(ascii_lowercase) for _ in range(random.randint(5, 20))])
        for _ in range(n)
    ]


@pytest.fixture
def data():
    return pd.DataFrame(data={"col_1": random_strings(1000)})


def test_split(data):
    train, validate, test = StableDataSplitter.split(
        data, "rxn", "col_1", split_ratio=0.05
    )
    assert len(train) == 889
    assert len(validate) == 48
    assert len(test) == 63
    assert validate.index.tolist() == [
        18,
        24,
        54,
        60,
        74,
        160,
        173,
        174,
        183,
        220,
        235,
        274,
        307,
        316,
        322,
        336,
        401,
        409,
        426,
        455,
        475,
        518,
        527,
        536,
        570,
        579,
        582,
        612,
        678,
        704,
        707,
        712,
        714,
        742,
        781,
        806,
        809,
        841,
        848,
        856,
        859,
        866,
        918,
        937,
        941,
        948,
        967,
        994,
    ]


def test_split_with_different_hash_seed(data):
    train1, validate1, test1 = StableDataSplitter.split(data, "rxn", "col_1")
    train2, validate2, test2 = StableDataSplitter.split(
        data, "rxn", "col_1", hash_seed=123
    )

    # The generated splits must be different if the seed was different
    assert not train2.equals(train1)
    assert not validate2.equals(validate1)
    assert not test2.equals(test1)

    # Check that the size of the splits is in the accepted
    assert len(train2) == 899
    assert len(validate2) == 58
    assert len(test2) == 43


def all_samples_in_one_split(splits: Sequence[pd.DataFrame]) -> bool:
    assert len(splits) == 3
    number_non_empty_splits = 0
    for split in splits:
        if len(split) != 0:
            number_non_empty_splits += 1
    return number_non_empty_splits == 1


def test_split_on_products():
    # With products splitting, all the reactions end up in the same split
    df = pd.DataFrame(
        data={"col_1": ["C>>CCC", "CC>>CCC", "CCC>>CCC", "CCCC>>CCC", "CCCCC>>CCC"]}
    )
    splits = StableDataSplitter.split(
        df, reaction_column_name="col_1", index_column="products"
    )
    assert all_samples_in_one_split(splits)

    # With normal splitting, does not all end up in the same split
    splits = StableDataSplitter.split(
        df, reaction_column_name="col_1", index_column="col_1"
    )
    assert not all_samples_in_one_split(splits)

    # With precursors splitting, does not all end up in the same split
    splits = StableDataSplitter.split(
        df, reaction_column_name="col_1", index_column="precursors"
    )
    assert not all_samples_in_one_split(splits)


def test_split_on_reactants():
    # With precursors splitting, all the reactions end up in the same split
    df = pd.DataFrame(
        data={
            "col_1": [
                "C.C.C>>C",
                "C.C.C>>CC",
                "C.C.C>>CCC",
                "C.C.C.C>>CCCC",
                "C.C.C>>CCCCC",
            ]
        }
    )
    splits = StableDataSplitter.split(
        df, reaction_column_name="col_1", index_column="precursors"
    )
    assert all_samples_in_one_split(splits)

    # With normal splitting, does not all end up in the same split
    splits = StableDataSplitter.split(
        df, reaction_column_name="col_1", index_column="col_1"
    )
    assert not all_samples_in_one_split(splits)

    # With products splitting, does not all end up in the same split
    splits = StableDataSplitter.split(
        df, reaction_column_name="col_1", index_column="products"
    )
    assert not all_samples_in_one_split(splits)


def test_train_split_is_shuffled():
    df = pd.DataFrame(data={"col_1": random_strings(1000)})

    train, validate, test = StableDataSplitter.split(
        df, "rxn", "col_1", shuffle_seed=123
    )

    # The train indices should not be sorted, while the validation and test
    # indices are sorted (i.e. for validation and test, the order is the same as
    # in the original data).
    train_indices = train.index.tolist()
    validation_indices = validate.index.tolist()
    test_indices = test.index.tolist()
    assert sorted(train_indices) != train_indices
    assert sorted(validation_indices) == validation_indices
    assert sorted(test_indices) == test_indices

    # repeating with the same seed leads to the same order
    train_2, _, _ = StableDataSplitter.split(df, "rxn", "col_1", shuffle_seed=123)
    assert train_2.equals(train)

    # repeating with a different seed leads to a different order, but the same content
    train_3, _, _ = StableDataSplitter.split(df, "rxn", "col_1", shuffle_seed=234)
    assert not train_3.equals(train)
    assert sorted(train_3.index.tolist()) == sorted(train.index.tolist())
