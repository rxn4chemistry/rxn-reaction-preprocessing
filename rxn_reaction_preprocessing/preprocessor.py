# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
""" The preprocessor class abstracts the workflow for preprocessing reaction data sets. """
import typing
from collections import Counter
from pathlib import Path
from typing import Any
from typing import Callable
from typing import List

import numpy as np
import pandas as pd
from rdkit import RDLogger
from tabulate import tabulate

from .mixed_reaction_filter import MixedReactionFilter
from .reaction import Reaction
from .reaction import ReactionPart
from rxn_reaction_preprocessing.cleaner import remove_isotope_information
from rxn_reaction_preprocessing.config import PreprocessConfig


class Preprocessor:

    def __init__(
        self,
        df: pd.DataFrame,
        reaction_column_name: str,
        valid_column: str = '_rxn_valid',
        valid_message_column: str = '_rxn_valid_messages',
        fragment_bond: str = '.',
        clean_data=True,
    ):
        """Creates a new instance of the Preprocessor class.

        Args:
            df: A pandas DataFrame containing the reaction SMARTS.
            reaction_column_name: The name of the DataFrame column containing the reaction SMARTS.
            valid_column: The name of the column to write the validation results to (will be created if it doesn't exist). Defaults to "_rxn_valid".
            valid_message_column: The name of the column to write the validation messages to in verbose mode (will be created if it doesn't exist). Defaults to "_rxn_valid_messages".
            fragment_bond: The token that represents fragment bonds in the raction SMILES.
            clean_data: Whether to run a simple pre-cleaning of the input data (naive check for valid SMARTS reactions based on the number of greater-thans in the string). Defaults to True.
        """
        self.df = df
        self.__reaction_column_name = reaction_column_name
        self.__valid_column = valid_column
        self.__valid_message_column = valid_message_column
        self.__fragment_bond = fragment_bond

        if clean_data:
            self.__clean_data()

    #
    # Private Methods
    #
    def __clean_data(self):
        """Drops records from the internal pandas DataFrame that do not contain valid reaction SMARTS (WARNING: Very naive, just checks for two greater-thans)"""
        self.df.drop(
            self.df[self.df[self.__reaction_column_name].str.count('>') != 2].index,
            inplace=True,
        )

    def __filter_func(
        self,
        reaction: Reaction,
        filter: MixedReactionFilter,
    ) -> bool:
        """The default filter function.

        Args:
            reaction: An input reaction.
            filter: An instance of MixedReactionFilter to test the reaction against.

        Returns:
            A boolean indicating whether or not the reaction is valid according to the supplied parameters.
        """
        # If reactions contains *any* None values for molecules
        if reaction.has_none():
            return False

        return filter.validate(reaction)

    def __filter_func_verbose(
        self,
        reaction: Reaction,
        filter: MixedReactionFilter,
    ) -> List[str]:
        """The default verbose filter function.

        Args:
            reaction: An input reaction.
            filter: An instance of MixedReactionFilter to test the reaction against.

        Returns:
            A list of reasons for an invalid reaction. An empty list for a valid reaction.
        """
        invalid_reasons = []

        # If reactions contains *any* None values for molecules
        if reaction.has_none():
            invalid_reasons.append('rdkit_molfromsmiles_failed')

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
        filter_func: Callable[[Reaction, MixedReactionFilter], bool] = None,
        filter_func_verbose: Callable[[Reaction, MixedReactionFilter], List[str]] = None,
    ):
        """Applies filter functions to the reaction. Alternatives for the default filter functions can be supplied. The default filters remove reactions containing molecules that could not be parse by rdkit's MolFromSmiles.

        Args:
            filter: An instance of MixedReactionFilter to be applied to each reaction.
            verbose: Whether or not to report the reasons for a reaction being deemed invalid. Defaults to False.
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
                self.df[self.__valid_message_column] = np.empty((len(self.df), 0)).tolist()

            self.df[self.__valid_message_column] += self.df.apply(
                lambda row: filter_func_verbose(
                    Reaction(row[self.__reaction_column_name], fragment_bond=self.__fragment_bond),
                    filter
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
                    lambda row: filter_func(
                        Reaction(
                            row[self.__reaction_column_name], fragment_bond=self.__fragment_bond
                        ), filter
                    ),
                    axis=1,
                ),
                self.df[self.__valid_column],
            )

        return self

    def filter_smarts(
        self,
        pattern: str,
        reaction_part: ReactionPart,
        keep: bool = False,
        verbose: bool = False,
    ):
        """Filter based on a SMARTS pattern. Unless keep is set to True, filters out reaction where the supplied pattern was found.

        Args:
            pattern: A valid SMARTS pattern.
            reaction_part: Which part of reaction to apply the SMARTS pattern to.
            keep: Whether to mark non-matches as invalid. Defaults to False.
            verbose: Whether to report the reason for why a reaction is deemed invalid. Defaults to False.

        Returns:
            Itself.
        """
        if self.__valid_column not in self.df.columns:
            self.df[self.__valid_column] = True

        tmp = np.logical_and(
            self.df.apply(
                lambda row: not (
                    (
                        len(
                            Reaction(
                                row[self.__reaction_column_name],
                                fragment_bond=self.__fragment_bond
                            ).find_in(pattern, reaction_part)
                        ) > 0
                    ) != keep
                ),
                axis=1,
            ),
            self.df[self.__valid_column],
        )

        self.df[self.__valid_column] = np.logical_and(tmp, self.df[self.__valid_column])

        if verbose:
            if self.__valid_message_column not in self.df.columns:
                self.df[self.__valid_message_column] = np.empty((len(self.df), 0)).tolist()

            if keep:
                tmp = ~tmp

            self.df.loc[tmp, self.__valid_message_column] = self.df.loc[
                tmp, self.__valid_message_column].apply(lambda c: c + ['pattern_' + pattern])

        return self

    def apply(
        self,
        func: Callable[[Reaction], Reaction],
        remove_duplicate_molecules: bool = False,
        **kwargs: Any,
    ):
        """Applies the supplied function to each reaction.

        Args:
            func: A function which is applied to each reaction.
            remove_duplicate_molecules: Whether to remove duplicate molecules when a reaction instance is created from a reaction SMARTS. Defaults to False.
            smiles_to_mol_kwargs: Additional parameters for the rdkit method MolFromSmiles.

        Returns:
            Itself.
        """
        kwargs.setdefault('canonical', True)

        self.df[self.__reaction_column_name] = self.df[self.__reaction_column_name].apply(
            lambda rxn: str(
                func(
                    Reaction(
                        rxn,
                        remove_duplicates=remove_duplicate_molecules,
                        fragment_bond=self.__fragment_bond,
                        **kwargs
                    )
                )
            )
        )

        return self

    def remove_duplicates(self):
        """A wrapper around pandas' drop_duplicates with the argument subset set to the reaction column name.

        Returns:
            Itself.
        """
        self.df.drop_duplicates(subset=[self.__reaction_column_name], inplace=True)
        return self

    def remove_invalids(self):
        """A wrapper around removing invalid options from pandas DataFrame.

        Returns:
            Itself.
        """
        self.df.drop(
            self.df[self.df[self.__valid_column] == False].index,  # noqa: E712
            inplace=True
        )
        return self

    def print_stats(self):
        """Prints statistics of the filtration to stdout.

        Returns:
            Itself.
        """
        print(f'- {len(self.df)} total reactions.')
        if self.__valid_column in self.df.columns:
            counts = self.df[self.__valid_column].value_counts()
            if True in counts:
                print(f'\033[92m- {counts[True]} valid reactions.\033[0m')
            if False in counts:
                print(f'\033[93m- {counts[False]} invalid reactions removed.\033[0m')

        if self.__valid_message_column in self.df.columns:
            reasons: typing.Counter[str] = Counter()
            for _, value in self.df[self.__valid_message_column].items():
                reasons.update(value)

            if len(reasons) > 0:
                print(
                    f'\033[93m- The {counts[False]} reactions were removed for the following reasons:'
                )
                headers = ['Reason', 'Number of Reactions']
                print(
                    tabulate(list(Counter(reasons).items()), headers, tablefmt='fancy_grid') +
                    '\033[0m'
                )

        return self

    #
    # Static Methods
    #
    @staticmethod
    def read_csv(filepath: str, reaction_column_name: str, fragment_bond: str = '.'):
        """A helper function to read a list or csv of reactions.

        Args:
            filepath: The path to the text file containing the reactions.
            reaction_column_name: The name of the reaction column (or the name that wil be given to the reaction column if the input file has no headers).
            fragment_bond: The token that represents fragment bonds in the raction SMILES.

        Returns:
            Preprocessor: A new preprocessor instance.
        """
        df = pd.read_csv(filepath, lineterminator='\n')
        if len(df.columns) == 1:
            df.rename(columns={df.columns[0]: reaction_column_name}, inplace=True)

        return Preprocessor(df, reaction_column_name, fragment_bond=fragment_bond)


def preprocess(cfg: PreprocessConfig) -> None:
    RDLogger.DisableLog('rdApp.*')

    output_file_path = Path(cfg.output_file_path)
    if not Path(cfg.input_file_path).exists():
        raise ValueError(f'Input file for preprocessing does not exist: {cfg.input_file_path}')

    # This is the function that is applied to each reaction.
    def merge_reactants_and_reagents(reaction: Reaction) -> Reaction:
        # Move agents to reactants
        reaction.remove_none()
        reaction.reactants.extend(reaction.agents)
        reaction.agents = []
        return reaction

    def remove_precursors_from_products_and_sort(reaction: Reaction) -> Reaction:
        # Remove products that are also reactants
        reaction.remove_precursors_from_products()
        return reaction.sort()

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
        cfg.input_file_path, cfg.reaction_column_name, fragment_bond=cfg.fragment_bond.value
    )

    # Remove duplicate reactions (useful for large dataset, this step is repeated later)
    print(f'\033[92m- {len(pp.df)} initial reactions.\033[0m')
    pp.remove_duplicates()
    print(f'\033[92m- {len(pp.df)} reactions after first deduplication.\033[0m')

    # In a first step the data is cleaned, in this case isotope information is removed
    pp.df[cfg.reaction_column_name
          ] = pp.df[cfg.reaction_column_name].apply(remove_isotope_information)

    # Apply the two functions above to all reactions, the remove_duplicate_molecules argument
    # is set to true to remove duplicate molecules within each reaction part
    pp.apply(merge_reactants_and_reagents)
    pp.apply(remove_precursors_from_products_and_sort, remove_duplicate_molecules=True)

    # Remove duplicate reactions
    pp.remove_duplicates()
    print(f'\033[92m- {len(pp.df)} reactions after second deduplication.\033[0m')

    # Apply the mixed reaction filter instance defined above, enable verbose mode
    pp.filter(mrf, True)

    # Print the detailed stats
    pp.print_stats()

    # Drop the invalid reactions
    pp.remove_invalids()

    # After dropping invalid columns, display stats again (as an example)
    pp.df.to_csv(output_file_path)
