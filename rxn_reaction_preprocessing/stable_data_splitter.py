# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
""" A utility class to split data sets in a stable manner. """
import functools
import os
from pathlib import Path
from typing import Tuple, Optional

import pandas as pd
from xxhash import xxh64_intdigest

from rxn_reaction_preprocessing.config import SplitConfig


class StableDataSplitter:

    @staticmethod
    def split(
        df: pd.DataFrame,
        reaction_column_name: str,
        index_column: str,
        split_ratio: float = 0.05,
        max_in_valid: Optional[int] = None,
        seed: int = 0
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Creates a stable split into training, validation, and test sets. Returns a boolean mask
        for each set, to not duplicate the DataFrame while preserving the original.

        Args:
            df (pd.DataFrame): The pandas DataFrame to be split into training, validation, and
                test sets.
            reaction_column_name: Name of the reaction column for the data file.
            index_column: The name of the column used to generate the hash which ensures
                stable splitting.
            split_ratio: The split ratio. Defaults to 0.05.
            max_in_valid: Maximal number of reactions to keep in the validation set.
                This can be useful to avoid unnecessarily prolonging training. It
                is considered as an approximate limit (due to randomness in hashing).
                Defaults to no restriction.
            seed: seed to use for hashing. The default of 0 corresponds to the default value
                in the xxhash implementation.

        Returns:
            Tuple[pd.Series, pd.Series, pd.Series]: A tuple of pandas Series containing a boolean
                mask of the training, validation and testing set rows within the original DataFrame.
        """
        hash_fn = functools.partial(xxh64_intdigest, seed=seed)

        if index_column == 'products':
            df['hash'] = df[reaction_column_name].apply(lambda value: value.split('>>')[1]
                                                        ).apply(lambda value: hash_fn(value))
        elif index_column == 'precursors':
            df['hash'] = df[reaction_column_name].apply(lambda value: value.split('>>')[0]
                                                        ).apply(lambda value: hash_fn(value))
        elif index_column in df.columns:
            df['hash'] = df[index_column].apply(lambda value: hash_fn(value))
        else:
            raise KeyError(index_column)

        test_ratio = split_ratio
        valid_ratio = split_ratio
        if max_in_valid is not None:
            valid_ratio = max_in_valid / len(df)

        return (
            df['hash'].apply(lambda value: value >= (test_ratio + valid_ratio) * 2**64),
            df['hash'].apply(
                lambda value: (test_ratio * 2**64 <= value < (test_ratio + valid_ratio) * 2**64)
            ),
            df['hash'].apply(lambda value: value < test_ratio * 2**64),
        )


def split(cfg: SplitConfig) -> None:
    output_directory = Path(cfg.output_directory)
    if not Path(cfg.input_file_path).exists():
        raise ValueError(f'Input file for standardization does not exist: {cfg.input_file_path}')

    df = pd.read_csv(cfg.input_file_path, lineterminator='\n')
    # Split into train, validation, and test sets, but do not export yet
    train, validation, test = StableDataSplitter.split(
        df,
        reaction_column_name=cfg.reaction_column_name,
        index_column=cfg.index_column,
        split_ratio=cfg.split_ratio,
        seed=cfg.seed
    )

    # Get the file name without the extension
    stem = Path(cfg.input_file_path).stem

    # Example of exporting one of the sets
    df[train].to_csv(os.path.join(output_directory, stem + '.train.csv'), index=False)
    df[validation].to_csv(os.path.join(output_directory, stem + '.validation.csv'), index=False)
    df[test].to_csv(os.path.join(output_directory, stem + '.test.csv'), index=False)
