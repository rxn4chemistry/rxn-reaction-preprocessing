# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED

""" A utility class to augment the dataset files """
import math
import random
from enum import Enum
from enum import auto
from typing import List

import pandas as pd
from rxn_chemutils.smiles_randomization import (
    randomize_smiles_rotated, randomize_smiles_restricted, randomize_smiles_unrestricted
)

from rxn_reaction_preprocessing.smiles_tokenizer import SmilesTokenizer


class RandomType(Enum):
    molecules = auto()
    unrestricted = auto()
    restricted = auto()
    rotated = auto()


def molecules_permutation_given_index(molecules_list: List[str],
                                      permutation_index: int) -> List['str']:
    """
    https://stackoverflow.com/questions/5602488/random-picks-from-permutation-generator
    """
    molecules_list = molecules_list[:]
    for i in range(len(molecules_list) - 1):
        permutation_index, j = divmod(permutation_index, len(molecules_list) - i)
        molecules_list[i], molecules_list[i + j] = molecules_list[i + j], molecules_list[i]
    return molecules_list


class Augmenter:

    def __init__(self, df: pd.DataFrame, fragment_bond: str = '.'):
        """Creates a new instance of the Augmenter class.

        Args:
            df (pd.DataFrame): A pandas DataFrame containing the molecules SMILES.
            fragment_bond (str): The fragment bond tolen contained in the SMILES.
        """
        self.df = df
        self.tokenizer = SmilesTokenizer()
        self.fragment_bond = fragment_bond

    #
    # Private Methods
    #

    def __randomize_smiles(self, smiles: str, random_type: RandomType,
                           permutations: int) -> List[str]:
        """
        Randomizes a molecules SMILES string that might contain fragment bonds and returns a number of augmented versions
        of the SMILES equal to permutations.

        Args:
            smiles (str): The molecules SMILES to augment
            random_type (RandomType): The type of randomization to be applied.
            permutations (int): The number of permutations to deliver for the SMILES

        Returns:
            List[str]: The list of randomized SMILES

        """
        # Raise for empty SMILES
        if not smiles:
            raise ValueError

        list_of_smiles: List[str] = []
        for i in range(permutations):
            list_of_smiles.append(
                '.'.join(
                    sorted(
                        [
                            self.fragment_bond.join(
                                [
                                    Augmenter.__randomize_smiles_without_fragment(
                                        fragment, random_type
                                    ) for fragment in group.split(self.fragment_bond)
                                ]
                            ) for group in smiles.split('.')
                        ]
                    )
                )
            )
        return list_of_smiles

    #
    # Private Static Methods
    #

    @staticmethod
    def __randomize_smiles_without_fragment(smiles: str, random_type: RandomType) -> str:
        """
        Generates a random version of a SMILES without a fragment bond

        Args:
            smiles (str): The pandas DataFrame to be split into training, validation, and test sets.
            random_type (RandomType): The type of randomization to be applied.

        Raises:
            InvalidSmiles: for invalid SMILES (raised via rxn_chemutils).
            ValueError: if an invalid randomization type is provided.

        Returns:
            str: the randomized SMILES
        """
        if random_type == RandomType.unrestricted:
            return randomize_smiles_unrestricted(smiles)
        elif random_type == RandomType.restricted:
            return randomize_smiles_restricted(smiles)
        elif random_type == RandomType.rotated:
            return randomize_smiles_rotated(smiles, with_order_reversal=True)
        raise ValueError(f'Invalid random type: {random_type}')

    @staticmethod
    def __randomize_molecules(smiles: str, permutations: int) -> List[str]:
        """
        Randomizes the order of the molecules inside a SMILES string that might contain fragment bonds and returns
        a number of augmented versions of the SMILES equal to permutations.
        For a number of molecules smaller than permutations, returns a number of permutations equal to the number of molecules

        Args:
            smiles (str): The molecules SMILES to augment
            permutations (int): The number of permutations to deliver for the SMILES

        Returns:
            List[str]: The list of randomized SMILES

        """
        # Raise for empty SMILES
        if not smiles:
            raise ValueError

        molecules_list = smiles.split('.')
        total_permutations = range(math.factorial(len(molecules_list)))
        permutation_indices = random.sample(
            total_permutations, min(permutations, len(molecules_list))
        )
        permuted_molecules_smiles = []
        for idx in permutation_indices:
            permuted_precursors = molecules_permutation_given_index(molecules_list, idx)
            permuted_molecules_smiles.append('.'.join(permuted_precursors))

        return permuted_molecules_smiles

    #
    # Public Methods
    #

    def augment(
        self,
        random_type: RandomType = RandomType.unrestricted,
        permutations: int = 1,
        tokenize: bool = True
    ) -> pd.DataFrame:
        """
        Creates samples for the augmentation. Returns a a pandas Series containing the augmented samples.

        Args:
            random_type (RandomType): The string identifying the type of randomization to apply.
                              "molecules" for randomization of the molecules (canonical SMILES kept)
                              "unrestricted" for unrestricted randomization
                              "restricted" for restricted randomization
                              "rotated" for rotated randomization
                              For details on the differences:
                              https://github.com/undeadpixel/reinvent-randomized and https://github.com/GLambard/SMILES-X
            permutations (int): The number of permutations to generate for each SMILES
            tokenize (bool): Whether to return the tokenized version of the SMILES

        Returns:
            pd.DataFrame: A pandas Series containing the augmented samples.
        """

        self.df[f'{random_type.name}'] = self.df.apply(lambda smiles: smiles.replace(' ', ''))

        if random_type != RandomType.molecules:
            self.df[f'{random_type.name}'] = self.df['smiles'].apply(
                lambda smiles: self.__randomize_smiles(smiles, random_type, permutations)
            )
        else:
            self.df[f'{random_type.name}'] = self.df['smiles'].apply(
                lambda smiles: self.__randomize_molecules(smiles, permutations)
            )
        if tokenize:
            self.df[f'{random_type.name}'] = self.df[f'{random_type.name}'].apply(
                lambda smiles: [self.tokenizer.tokenize(smi) for smi in smiles]
            )
        return self.df.explode(f'{random_type.name}').reset_index(drop=True)

    #
    # Public Static Methods
    #

    @staticmethod
    def read_csv(filepath: str, fragment_bond: str = '.'):
        """A helper function to read a list or csv of SMILES.

        Args:
            filepath (str): The path to the text file containing the molecules SMILES.
            fragment_bond (str): The fragment token in the reaction SMILES

        Returns:
            Augmenter: A new augmenter instance.
        """
        df = pd.read_csv(filepath, lineterminator='\n')
        if len(df.columns) == 1:
            df.rename(columns={df.columns[0]: 'smiles'}, inplace=True)

        return Augmenter(df, fragment_bond)
