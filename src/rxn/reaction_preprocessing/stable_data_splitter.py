# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2022
# ALL RIGHTS RESERVED
"""A utility class to split data sets in a stable manner."""
import csv
import functools
from pathlib import Path
from typing import Callable, Hashable, Iterable, List

from rxn.utilities.csv import CsvIterator
from rxn.utilities.files import PathLike, stable_shuffle
from typing_extensions import Protocol
from xxhash import xxh64_intdigest

from rxn.reaction_preprocessing.config import SplitConfig
from rxn.reaction_preprocessing.utils import DataSplit


class _CsvWriter(Protocol):
    """Useful because csv.writer can't be used as a type annotation."""

    def writerow(self, row: List[str]) -> None: ...

    def writerows(self, rows: Iterable[List[str]]) -> None: ...


class StableSplitter:
    """
    Split data in a reproducible manner, based on the hash of values required
    to always be in the same split.

    Useful for instance to ensure that a reaction product with a given SMILES
    will always be in the same split.
    """

    HASH_SIZE = 2**64

    def __init__(
        self,
        split_ratio: float,
        seed: int = 0,
    ):
        """
        Args:
            split_ratio: The approximate split ratio for test and validation set.
            seed: seed to use for hashing. The default of 0 corresponds to the
                default value in the xxhash implementation.
        """
        self.hash_fn = functools.partial(xxh64_intdigest, seed=seed)

        self.test_ratio = split_ratio
        self.valid_ratio = split_ratio

        # Compute these here to avoid repeating the calculations all the time
        # in the get_split function
        self._test_threshold = self.test_ratio * self.HASH_SIZE
        self._validation_threshold = (
            self.test_ratio + self.valid_ratio
        ) * self.HASH_SIZE

    def get_split(self, split_value: Hashable) -> DataSplit:
        value = self.hash_fn(split_value)  # type:ignore
        if value < self._test_threshold:
            return DataSplit.TEST
        if value < self._validation_threshold:
            return DataSplit.VALIDATION
        return DataSplit.TRAIN


class StableDataSplitter:
    def __init__(
        self,
        reaction_column_name: str,
        index_column: str,
        split_ratio: float = 0.05,
        hash_seed: int = 0,
        shuffle_seed: int = 42,
    ):
        """
        Args:
            reaction_column_name: Name of the reaction column for the data file.
            index_column: The name of the column used to generate the hash which ensures
                stable splitting. "products" and "precursors" are also allowed even if
                they do not exist as columns.
            split_ratio: The split ratio. Defaults to 0.05.
            hash_seed: seed to use for hashing. The default of 0 corresponds to
                the default value in the xxhash implementation.
            shuffle_seed: Seed for shuffling the train split.
        """
        self.rxn_column = reaction_column_name
        self.index_column = index_column
        self.split_ratio = split_ratio
        self.hash_seed = hash_seed
        self.shuffle_seed = shuffle_seed

    def split_file(
        self,
        input_csv: PathLike,
        train_csv: PathLike,
        valid_csv: PathLike,
        test_csv: PathLike,
    ) -> None:
        """
        Split an input file into train, validation, and test CSVs.
        """
        with open(input_csv, "rt") as f_input, open(train_csv, "wt") as f_train, open(
            valid_csv, "wt"
        ) as f_valid, open(test_csv, "wt") as f_test:
            input_iterator = CsvIterator.from_stream(f_input)

            # initialize the writers
            writers = []
            for f in [f_train, f_valid, f_test]:
                writer = csv.writer(f)
                writer.writerow(input_iterator.columns)
                writers.append(writer)

            self._split_iterator(input_iterator, *writers)

        # Shuffle the training split
        stable_shuffle(train_csv, train_csv, seed=self.shuffle_seed, is_csv=True)

    def _split_iterator(
        self,
        input_iterator: CsvIterator,
        train_writer: _CsvWriter,
        valid_writer: _CsvWriter,
        test_writer: _CsvWriter,
    ) -> None:
        """
        Split an input CSV iterator into train, validation, and test iterators,
        by writing immediately to the corresponding CSV writers.
        """
        fn = self._callable_for_value_to_hash(input_iterator)
        splitter = StableSplitter(
            split_ratio=self.split_ratio,
            seed=self.hash_seed,
        )

        for row in input_iterator.rows:
            value_to_hash = fn(row)
            split = splitter.get_split(value_to_hash)
            if split is DataSplit.TRAIN:
                train_writer.writerow(row)
            elif split is DataSplit.VALIDATION:
                valid_writer.writerow(row)
            elif split is DataSplit.TEST:
                test_writer.writerow(row)

    def _callable_for_value_to_hash(
        self, csv_iterator: CsvIterator
    ) -> Callable[[List[str]], Hashable]:
        if self.index_column == "products":
            rxn_column = csv_iterator.column_index(self.rxn_column)
            return lambda x: x[rxn_column].split(">>")[1]
        elif self.index_column == "precursors":
            rxn_column = csv_iterator.column_index(self.rxn_column)
            return lambda x: x[rxn_column].split(">>")[0]
        elif self.index_column in csv_iterator.columns:
            column_index = csv_iterator.column_index(self.index_column)
            return lambda x: x[column_index]
        raise RuntimeError(
            f'Can\'t determine what value to hash from index_column "{self.index_column}".'
        )


def split(cfg: SplitConfig) -> None:
    output_directory = Path(cfg.output_directory)
    if not Path(cfg.input_file_path).exists():
        raise ValueError(
            f"Input file for standardization does not exist: {cfg.input_file_path}"
        )

    splitter = StableDataSplitter(
        reaction_column_name=cfg.reaction_column_name,
        index_column=cfg.index_column,
        hash_seed=cfg.hash_seed,
        split_ratio=cfg.split_ratio,
        shuffle_seed=cfg.shuffle_seed,
    )

    # Get the file name without the extension
    stem = Path(cfg.input_file_path).stem

    splitter.split_file(
        input_csv=cfg.input_file_path,
        train_csv=output_directory / (stem + ".train.csv"),
        valid_csv=output_directory / (stem + ".validation.csv"),
        test_csv=output_directory / (stem + ".test.csv"),
    )
