# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
""" A utility class to split data sets in a stable manner. """
import functools
import os
from pathlib import Path
from typing import Hashable, Optional, Tuple

import pandas as pd
from xxhash import xxh64_intdigest

from rxn.reaction_preprocessing.config import SplitConfig
from rxn.reaction_preprocessing.utils import DataSplit


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
        total_size: Optional[int] = None,
        seed: int = 0,
    ):
        """
        Args:
            split_ratio: The approximate split ratio for test and validation set.
            total_size: Total size of the dataset, must be provided only if
                max_if_valid is not None.
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
        value = self.hash_fn(split_value)
        if value < self._test_threshold:
            return DataSplit.TEST
        if value < self._validation_threshold:
            return DataSplit.VALIDATION
        return DataSplit.TRAIN


class StableDataSplitter:
    @staticmethod
    def split(
        df: pd.DataFrame,
        reaction_column_name: str,
        index_column: str,
        split_ratio: float = 0.05,
        hash_seed: int = 0,
        shuffle_seed: int = 42,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Creates a stable split into training, validation, and test sets. Returns a boolean mask
        for each set, to not duplicate the DataFrame while preserving the original.

        Args:
            df: The pandas DataFrame to be split into training, validation, and test sets.
            reaction_column_name: Name of the reaction column for the data file.
            index_column: The name of the column used to generate the hash which ensures
                stable splitting.
            split_ratio: The split ratio. Defaults to 0.05.
            hash_seed: seed to use for hashing. The default of 0 corresponds to
                the default value in the xxhash implementation.
            shuffle_seed: Seed for shuffling the train split.

        Returns:
            A tuple of three pandas DataFrames for the three splits.
        """
        splitter = StableSplitter(
            split_ratio=split_ratio,
            total_size=len(df),
            seed=hash_seed,
        )

        # Get a pandas Series for the splits by first getting a pandas Series
        # for the values to hash, and then do the actual hashing / split attribution.
        pre_hash_series = StableDataSplitter._pre_hashing_series(
            df, index_column=index_column, reaction_column_name=reaction_column_name
        )
        data_split = pre_hash_series.apply(splitter.get_split)

        train_df = df[data_split == DataSplit.TRAIN]
        validation_df = df[data_split == DataSplit.VALIDATION]
        test_df = df[data_split == DataSplit.TEST]

        # Shuffle the training split
        train_df = train_df.sample(frac=1, random_state=shuffle_seed)

        return train_df, validation_df, test_df

    @staticmethod
    def _pre_hashing_series(
        df: pd.DataFrame, index_column: str, reaction_column_name: str
    ) -> pd.Series:
        if index_column == "products":
            return df[reaction_column_name].apply(lambda value: value.split(">>")[1])
        elif index_column == "precursors":
            return df[reaction_column_name].apply(lambda value: value.split(">>")[0])
        elif index_column in df.columns:
            return df[index_column]
        else:
            raise KeyError(index_column)


def split(cfg: SplitConfig) -> None:
    output_directory = Path(cfg.output_directory)
    if not Path(cfg.input_file_path).exists():
        raise ValueError(
            f"Input file for standardization does not exist: {cfg.input_file_path}"
        )

    df = pd.read_csv(cfg.input_file_path, lineterminator="\n")
    # Split into train, validation, and test sets, but do not export yet
    train, validation, test = StableDataSplitter.split(
        df,
        reaction_column_name=cfg.reaction_column_name,
        index_column=cfg.index_column,
        split_ratio=cfg.split_ratio,
        hash_seed=cfg.hash_seed,
        shuffle_seed=cfg.shuffle_seed,
    )

    # Get the file name without the extension
    stem = Path(cfg.input_file_path).stem

    # Example of exporting one of the sets
    train.to_csv(os.path.join(output_directory, stem + ".train.csv"), index=False)
    validation.to_csv(
        os.path.join(output_directory, stem + ".validation.csv"), index=False
    )
    test.to_csv(os.path.join(output_directory, stem + ".test.csv"), index=False)
