# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
""" A utility class to split data sets in a stable manner. """
from typing import Tuple

import pandas as pd
from xxhash import xxh64_intdigest


class StableDataSplitter:

    @staticmethod
    def split(df: pd.DataFrame,
              index_column: str,
              split_ratio: float = 0.05) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Creates a stable split into training, validation, and test sets. Returns a boolean mask for each set, to not duplicate the DataFrame while preserving the original.

        Args:
            df (pd.DataFrame): The pandas DataFrame to be split into training, validation, and test sets.
            index_column (str): The name of the column used to generate the hash which ensures stable splitting.
            split_ratio (float, optional): The split ratio. Defaults to 0.05.

        Returns:
            Tuple[pd.Series, pd.Series, pd.Series]: A tuple of pandas Series containing a boolean mask of the training, validation and testing set rows within the original DataFrame.
        """

        df['hash'] = df[index_column].apply(lambda value: xxh64_intdigest(value))

        return (
            df['hash'].apply(lambda value: value >= split_ratio * 2 * 2**64),
            df['hash'].apply(
                lambda value: (split_ratio * 2**64 <= value) and (value < split_ratio * 2 * 2**64)
            ),
            df['hash'].apply(lambda value: value < split_ratio * 2**64),
        )
