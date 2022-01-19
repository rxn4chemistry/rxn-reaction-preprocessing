# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
import random
from string import ascii_lowercase
from typing import List

import pandas as pd
import pytest

from rxn_reaction_preprocessing import StableDataSplitter
from rxn_reaction_preprocessing.utils import reset_random_seed


def random_strings(n: int) -> List[str]:
    reset_random_seed()
    return [
        ''.join([random.choice(ascii_lowercase) for _ in range(random.randint(5, 20))])
        for _ in range(n)
    ]


@pytest.fixture
def data():
    return pd.DataFrame(data={'col_1': random_strings(1000)})


def test_split(data):
    train, validate, test = StableDataSplitter.split(data, 'rxn', 'col_1', split_ratio=0.05)
    assert train.sum() == 889
    assert validate.sum() == 48
    assert test.sum() == 63
    assert validate.index[validate.tolist()].tolist() == [
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


def test_split_with_different_seed(data):
    train1, validate1, test1 = StableDataSplitter.split(data, 'rxn', 'col_1', split_ratio=0.05)
    train2, validate2, test2 = StableDataSplitter.split(
        data, 'rxn', 'col_1', split_ratio=0.05, seed=123
    )

    # The generated splits must be different if the seed was different
    assert not train2.equals(train1)
    assert not validate2.equals(validate1)
    assert not test2.equals(test1)

    # Check that the size of the splits is in the accepted
    assert train2.sum() == 899
    assert validate2.sum() == 58
    assert test2.sum() == 43


def test_split_with_max_valid_samples():
    df = pd.DataFrame(data={'col_1': random_strings(10000)})

    # Split without restricting the size - leads to roughly 5% in valid and 5% in test
    train, validate, test = StableDataSplitter.split(
        df, 'rxn', 'col_1', split_ratio=0.05, max_in_valid=None
    )
    assert 8900 < train.sum() < 9100
    assert 450 < validate.sum() < 550
    assert 450 < test.sum() < 550

    # Split limiting the size of the validation set to +/- 50 samples
    train, validate, test = StableDataSplitter.split(
        df, 'rxn', 'col_1', split_ratio=0.05, max_in_valid=50
    )
    assert 9400 < train.sum() < 9600
    assert 40 < validate.sum() < 60
    assert 450 < test.sum() < 550
