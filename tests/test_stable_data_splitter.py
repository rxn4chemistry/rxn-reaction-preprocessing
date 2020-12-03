# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
import random
from string import ascii_lowercase

import pandas as pd
import pytest

from data_preprocessor import StableDataSplitter

random.seed(42)


@pytest.fixture
def data():
    random_strings = [
        ''.join([random.choice(ascii_lowercase) for i in range(random.randint(5, 20))])
        for j in range(1000)
    ]

    return pd.DataFrame(data={'col_1': random_strings})


def test_split(data):
    train, validate, test = StableDataSplitter.split(data, 'col_1', split_ratio=0.05)
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
