# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2022
# ALL RIGHTS RESERVED
import random
from pathlib import Path
from string import ascii_lowercase
from typing import Generator, List, Tuple

import pandas as pd
import pytest
from rxn.utilities.csv import iterate_csv_column
from rxn.utilities.files import named_temporary_directory

from rxn.reaction_preprocessing import StableDataSplitter
from rxn.reaction_preprocessing.utils import reset_random_seed

RXN_COLUMN = "col_1"


def random_strings(n: int) -> List[str]:
    reset_random_seed()
    return [
        "".join([random.choice(ascii_lowercase) for _ in range(random.randint(5, 20))])
        for _ in range(n)
    ]


class SplitsDirectory:
    """Class to return as a fixture, determines paths to use in the tests, and
    populates the "input" data file."""

    def __init__(self, tmp_directory: Path, content: List[str]):
        self.directory = tmp_directory
        self.input_csv = self.directory / "input.csv"
        self.train_csv = self.directory / "train.csv"
        self.valid_csv = self.directory / "validation.csv"
        self.test_csv = self.directory / "test.csv"

        # Write random data to the input
        input_df = pd.DataFrame(
            data={RXN_COLUMN: content, "idx": list(range(len(content)))}
        )
        input_df.to_csv(self.input_csv, index=False)

    def train_content(self, column: str) -> List[str]:
        return list(iterate_csv_column(self.train_csv, column))

    def valid_content(self, column: str) -> List[str]:
        return list(iterate_csv_column(self.valid_csv, column))

    def test_content(self, column: str) -> List[str]:
        return list(iterate_csv_column(self.test_csv, column))

    def contents(self, column: str) -> Tuple[List[str], List[str], List[str]]:
        return (
            self.train_content(column),
            self.valid_content(column),
            self.test_content(column),
        )

    def contents_as_ints(self, column: str) -> Tuple[List[int], List[int], List[int]]:
        return (
            [int(v) for v in self.train_content(column)],
            [int(v) for v in self.valid_content(column)],
            [int(v) for v in self.test_content(column)],
        )


@pytest.fixture
def data_directory() -> Generator[SplitsDirectory, None, None]:
    with named_temporary_directory() as path:
        yield SplitsDirectory(path, content=random_strings(1000))


def test_split(data_directory: SplitsDirectory) -> None:
    # for conciseness
    d = data_directory

    splitter = StableDataSplitter("rxn", RXN_COLUMN, split_ratio=0.05)
    splitter.split_file(d.input_csv, d.train_csv, d.valid_csv, d.test_csv)

    train, valid, test = d.contents(RXN_COLUMN)
    assert len(train) == 889
    assert len(valid) == 48
    assert len(test) == 63
    assert [int(v) for v in d.valid_content("idx")] == [
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


def test_split_with_different_hash_seed(data_directory: SplitsDirectory) -> None:
    # for conciseness
    d = data_directory

    # First: use default value for the hash seed
    splitter = StableDataSplitter("rxn", RXN_COLUMN)
    splitter.split_file(d.input_csv, d.train_csv, d.valid_csv, d.test_csv)
    train1, valid1, test1 = d.contents(RXN_COLUMN)

    # Second: Change the hash seed
    splitter.hash_seed = 123
    splitter.split_file(d.input_csv, d.train_csv, d.valid_csv, d.test_csv)
    train2, valid2, test2 = d.contents(RXN_COLUMN)

    # The generated splits must be different if the seed was different
    assert train1 != train2
    assert valid1 != valid2
    assert test1 != test2

    # Check that the size of the splits is in the accepted
    assert len(train2) == 899
    assert len(valid2) == 58
    assert len(test2) == 43


def all_samples_in_one_split(test_directory: SplitsDirectory, column: str) -> bool:
    number_non_empty_splits = 0
    for split in test_directory.contents(column):
        if len(split) != 0:
            number_non_empty_splits += 1
    return number_non_empty_splits == 1


def test_split_on_products() -> None:
    splitter = StableDataSplitter(reaction_column_name=RXN_COLUMN, index_column="tbd")

    with named_temporary_directory() as path:
        d = SplitsDirectory(
            path, content=["C>>CCC", "CC>>CCC", "CCC>>CCC", "CCCC>>CCC", "CCCCC>>CCC"]
        )

        # With products splitting, all the reactions end up in the same split
        splitter.index_column = "products"
        splitter.split_file(d.input_csv, d.train_csv, d.valid_csv, d.test_csv)
        assert all_samples_in_one_split(d, RXN_COLUMN)

        # With normal splitting, does not all end up in the same split
        splitter.index_column = RXN_COLUMN
        splitter.split_file(d.input_csv, d.train_csv, d.valid_csv, d.test_csv)
        assert not all_samples_in_one_split(d, RXN_COLUMN)

        # With precursors splitting, does not all end up in the same split
        splitter.index_column = "precursors"
        splitter.split_file(d.input_csv, d.train_csv, d.valid_csv, d.test_csv)
        assert not all_samples_in_one_split(d, RXN_COLUMN)


def test_split_on_reactants() -> None:
    splitter = StableDataSplitter(reaction_column_name=RXN_COLUMN, index_column="tbd")

    with named_temporary_directory() as path:
        d = SplitsDirectory(
            path,
            content=[
                "C.C.C>>C",
                "C.C.C>>CC",
                "C.C.C>>CCC",
                "C.C.C.C>>CCCC",
                "C.C.C>>CCCCC",
            ],
        )

        # With precursors splitting, all the reactions end up in the same split
        splitter.index_column = "precursors"
        splitter.split_file(d.input_csv, d.train_csv, d.valid_csv, d.test_csv)
        assert all_samples_in_one_split(d, RXN_COLUMN)

        # With normal splitting, does not all end up in the same split
        splitter.index_column = RXN_COLUMN
        splitter.split_file(d.input_csv, d.train_csv, d.valid_csv, d.test_csv)
        assert not all_samples_in_one_split(d, RXN_COLUMN)

        # With products splitting, does not all end up in the same split
        splitter.index_column = "products"
        splitter.split_file(d.input_csv, d.train_csv, d.valid_csv, d.test_csv)
        assert not all_samples_in_one_split(d, RXN_COLUMN)


def test_train_split_is_shuffled(data_directory: SplitsDirectory) -> None:
    # for conciseness
    d = data_directory

    splitter = StableDataSplitter("rxn", RXN_COLUMN, shuffle_seed=123)
    splitter.split_file(d.input_csv, d.train_csv, d.valid_csv, d.test_csv)

    train, valid, test = d.contents_as_ints("idx")

    # The train indices should not be sorted, while the validation and test
    # indices are sorted (i.e. for validation and test, the order is the same as
    # in the original data).
    assert sorted(train) != train
    assert sorted(valid) == valid
    assert sorted(test) == test

    # repeating with the same seed leads to the same order
    splitter.split_file(d.input_csv, d.train_csv, d.valid_csv, d.test_csv)
    train_2, _, _ = d.contents_as_ints("idx")
    assert train_2 == train

    # repeating with a different seed leads to a different order, but the same content
    splitter.shuffle_seed = 234
    splitter.split_file(d.input_csv, d.train_csv, d.valid_csv, d.test_csv)
    train_3, _, _ = d.contents_as_ints("idx")
    assert train_3 != train
    assert sorted(train_3) == sorted(train)
