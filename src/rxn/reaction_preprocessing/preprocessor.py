# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2022
# ALL RIGHTS RESERVED
""" The preprocessor class abstracts the workflow for preprocessing reaction data sets. """
import logging
import typing
from collections import Counter
from pathlib import Path
from typing import Callable, List, Optional

import numpy as np
import pandas as pd
from rdkit import RDLogger
from rxn.chemutils.reaction_equation import ReactionEquation
from rxn.chemutils.reaction_smiles import parse_any_reaction_smiles
from tabulate import tabulate

from .config import PreprocessConfig
from .mixed_reaction_filter import MixedReactionFilter
from .reaction_standardizer import ReactionStandardizer

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class Preprocessor:
    def __init__(
        self,
        df: pd.DataFrame,
        reaction_column_name: str,
        valid_column: str = "_rxn_valid",
        valid_message_column: str = "_rxn_valid_messages",
        fragment_bond: str = ".",
        clean_data: bool = True,
    ):
        """Creates a new instance of the Preprocessor class.

        Args:
            df: A pandas DataFrame containing the reaction SMARTS.
            reaction_column_name: The name of the DataFrame column containing the reaction SMARTS.
            valid_column: The name of the column to write the validation results to
                (will be created if it doesn't exist). Defaults to "_rxn_valid".
            valid_message_column: The name of the column to write the validation
                messages to in verbose mode (will be created if it doesn't
                exist). Defaults to "_rxn_valid_messages".
            fragment_bond: The token that represents fragment bonds in the reaction SMILES.
            clean_data: Whether to run a simple pre-cleaning of the input data
                (naive check for valid SMARTS reactions based on the number of
                greater-thans in the string). Defaults to True.
        """
        self.df = df
        self.reaction_standardizer = ReactionStandardizer()
        self.__reaction_column_name = reaction_column_name
        self.__valid_column = valid_column
        self.__valid_message_column = valid_message_column
        self.__fragment_bond = fragment_bond

        if clean_data:
            self.__clean_data()

    #
    # Private Methods
    #
    def __clean_data(self) -> None:
        """Drops records from the internal pandas DataFrame that do not contain valid reaction
        SMARTS (WARNING: Very naive, just checks for two greater-thans)"""
        self.df.drop(
            self.df[self.df[self.__reaction_column_name].str.count(">") != 2].index,
            inplace=True,
        )

    def __filter_func(
        self,
        reaction_smiles: str,
        filter: MixedReactionFilter,
    ) -> bool:
        """The default filter function.

        Args:
            reaction_smiles: An input reaction SMILES.
            filter: An instance of MixedReactionFilter to test the reaction against.

        Returns:
            A boolean indicating whether or not the reaction is valid according
            to the supplied parameters.
        """

        reaction = ReactionEquation.from_string(
            reaction_smiles, fragment_bond=self.__fragment_bond
        )
        return filter.is_valid(reaction)

    def __filter_func_verbose(
        self,
        reaction_smiles: str,
        filter: MixedReactionFilter,
    ) -> List[str]:
        """The default verbose filter function.

        Args:
            reaction_smiles: An input reaction SMILES.
            filter: An instance of MixedReactionFilter to test the reaction against.

        Returns:
            A list of reasons for an invalid reaction. An empty list for a valid reaction.
        """
        invalid_reasons = []

        reaction = ReactionEquation.from_string(
            reaction_smiles, fragment_bond=self.__fragment_bond
        )
        _, filter_reasons = filter.validate_reasons(reaction)
        invalid_reasons.extend(filter_reasons)

        return invalid_reasons

    #
    # Public Methods
    #
    def filter(
        self,
        filter: MixedReactionFilter,
        verbose: bool = False,
        filter_func: Optional[Callable[[str, MixedReactionFilter], bool]] = None,
        filter_func_verbose: Optional[
            Callable[[str, MixedReactionFilter], List[str]]
        ] = None,
    ) -> "Preprocessor":
        """Applies filter functions to the reaction. Alternatives for the default filter functions
        can be supplied. The default filters remove reactions containing molecules that could not
        be parse by rdkit's MolFromSmiles.

        Args:
            filter: An instance of MixedReactionFilter to be applied to each reaction.
            verbose: Whether or not to report the reasons for a reaction being
                deemed invalid. Defaults to False.
            filter_func: A custom filter function. Defaults to None.
            filter_func_verbose: A custom verbose filter function. Defaults to None.

        Returns:
            Itself.
        """
        if filter_func is None:
            filter_func = self.__filter_func

        if filter_func_verbose is None:
            filter_func_verbose = self.__filter_func_verbose

        if self.__valid_column not in self.df.columns:
            self.df[self.__valid_column] = True

        if verbose:
            if self.__valid_message_column not in self.df.columns:
                self.df[self.__valid_message_column] = np.empty(
                    (len(self.df), 0)
                ).tolist()

            self.df[self.__valid_message_column] += self.df.apply(
                lambda row: filter_func_verbose(
                    row[self.__reaction_column_name], filter
                ),
                axis=1,
            )

            self.df[self.__valid_column] = np.logical_and(
                np.where(self.df[self.__valid_message_column], False, True),
                self.df[self.__valid_column],
            )
        else:
            self.df[self.__valid_column] = np.logical_and(
                self.df.apply(
                    lambda row: filter_func(row[self.__reaction_column_name], filter),
                    axis=1,
                ),
                self.df[self.__valid_column],
            )

        return self

    def reaction_standardization(self) -> None:
        """
        Reaction standardization, including merging of reactants and reactants,
        removal of duplicates, sorting, etc.

        Note that this standardization relies on the molecules in the reaction
        SMILES to be canonical - which is the case when this function is
        called as part of the full data processing pipeline.
        """

        def standardize_smiles(rxn_smiles: str) -> str:
            """Function standardizing the reaction SMILES directly,
            to pass to pandas.apply()."""
            reaction = parse_any_reaction_smiles(rxn_smiles)
            reaction = self.reaction_standardizer(reaction)
            return reaction.to_string(self.__fragment_bond)

        rxn_column = self.__reaction_column_name
        self.df[rxn_column] = self.df[rxn_column].apply(standardize_smiles)

    def remove_duplicates(self) -> "Preprocessor":
        """A wrapper around pandas' drop_duplicates with the argument subset
        set to the reaction column name.

        Returns:
            Itself.
        """
        self.df.drop_duplicates(subset=[self.__reaction_column_name], inplace=True)
        return self

    def remove_invalids(self) -> "Preprocessor":
        """A wrapper around removing invalid options from pandas DataFrame.

        Returns:
            Itself.
        """
        self.df.drop(
            self.df[self.df[self.__valid_column] == False].index,  # noqa: E712
            inplace=True,
        )
        return self

    def print_stats(self) -> "Preprocessor":
        """Prints statistics of the filtration to stdout.

        Returns:
            Itself.
        """
        logger.info(f"- {len(self.df)} total reactions.")
        if self.__valid_column in self.df.columns:
            counts = self.df[self.__valid_column].value_counts()
            if True in counts:
                logger.info(f"- {counts[True]} valid reactions.")
            if False in counts:
                logger.info(f"- {counts[False]} invalid reactions removed.")

        if self.__valid_message_column in self.df.columns:
            reasons: typing.Counter[str] = Counter()
            for _, value in self.df[self.__valid_message_column].items():
                reasons.update(value)

            if len(reasons) > 0:
                headers = ["Reason", "Number of Reactions"]
                logger.info(
                    f"- The {counts[False]} reactions were removed for the following reasons:\n"
                    f'{tabulate(list(Counter(reasons).items()), headers, tablefmt="fancy_grid")}'
                )

        return self

    #
    # Static Methods
    #
    @staticmethod
    def read_csv(
        filepath: str, reaction_column_name: str, fragment_bond: str = "."
    ) -> "Preprocessor":
        """A helper function to read a list or csv of reactions.

        Args:
            filepath: The path to the text file containing the reactions.
            reaction_column_name: The name of the reaction column (or the name that
                wil be given to the reaction column if the input file has no headers).
            fragment_bond: The token that represents fragment bonds in the raction SMILES.

        Returns:
            Preprocessor: A new preprocessor instance.
        """
        df = pd.read_csv(filepath, lineterminator="\n")
        if len(df.columns) == 1:
            df.rename(columns={df.columns[0]: reaction_column_name}, inplace=True)

        return Preprocessor(df, reaction_column_name, fragment_bond=fragment_bond)


def preprocess(cfg: PreprocessConfig) -> None:
    RDLogger.DisableLog("rdApp.*")

    output_file_path = Path(cfg.output_file_path)
    if not Path(cfg.input_file_path).exists():
        raise ValueError(
            f"Input file for preprocessing does not exist: {cfg.input_file_path}"
        )

    # Create a instance of the mixed reaction filter with default values.
    # Make arguments for all properties in script
    mrf = MixedReactionFilter(
        max_reactants=cfg.max_reactants,
        max_agents=cfg.max_agents,
        max_products=cfg.max_products,
        min_reactants=cfg.min_reactants,
        min_agents=cfg.min_agents,
        min_products=cfg.min_products,
        max_reactants_tokens=cfg.max_reactants_tokens,
        max_agents_tokens=cfg.max_agents_tokens,
        max_products_tokens=cfg.max_products_tokens,
        max_absolute_formal_charge=cfg.max_absolute_formal_charge,
    )

    pp = Preprocessor.read_csv(
        cfg.input_file_path,
        cfg.reaction_column_name,
        fragment_bond=cfg.fragment_bond.value,
    )

    # Remove duplicate reactions (useful for large dataset, this step is repeated later)
    logger.info(f"- {len(pp.df)} initial reactions.")
    pp.remove_duplicates()
    logger.info(f"- {len(pp.df)} reactions after first deduplication.")

    # Remove duplicate molecules, sort, etc.
    # NB: this relies on molecules in the reaction SMILES to be canonical already!
    pp.reaction_standardization()

    # Remove duplicate reactions
    pp.remove_duplicates()
    logger.info(f"- {len(pp.df)} reactions after second deduplication.")

    # Apply the mixed reaction filter instance defined above, enable verbose mode
    pp.filter(mrf, True)

    # Print the detailed stats
    pp.print_stats()

    # Drop the invalid reactions
    pp.remove_invalids()

    # After dropping invalid columns, display stats again (as an example)
    pp.df.to_csv(output_file_path, index=False)
