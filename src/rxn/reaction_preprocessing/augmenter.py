# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED
"""A utility class to augment the dataset files"""
import math
import random
from pathlib import Path
from typing import List, Set

import pandas as pd
from rxn.chemutils.smiles_randomization import (
    randomize_smiles_restricted,
    randomize_smiles_rotated,
    randomize_smiles_unrestricted,
)

from rxn.reaction_preprocessing.config import AugmentConfig
from rxn.reaction_preprocessing.smiles_tokenizer import SmilesTokenizer
from rxn.reaction_preprocessing.utils import RandomType, ReactionSection


def molecules_permutation_given_index(
    molecules_list: List[str], permutation_index: int
) -> List["str"]:
    """
    https://stackoverflow.com/questions/5602488/random-picks-from-permutation-generator
    """
    molecules_list = molecules_list[:]
    for i in range(len(molecules_list) - 1):
        permutation_index, j = divmod(permutation_index, len(molecules_list) - i)
        molecules_list[i], molecules_list[i + j] = (
            molecules_list[i + j],
            molecules_list[i],
        )
    return molecules_list


class Augmenter:
    """Augmenter.

    Note: Unlike the other classes, which are memory-efficient, this one loads the
    whole data in a pandas DataFrame for processing.
    """

    def __init__(
        self, df: pd.DataFrame, reaction_column_name: str, fragment_bond: str = "."
    ):
        """Creates a new instance of the Augmenter class.

        Args:
            df (pd.DataFrame): A pandas DataFrame containing the molecules SMILES.
            reaction_column_name: The name of the DataFrame column containing the reaction SMILES.
            fragment_bond (str): The fragment bond token contained in the SMILES.
        """
        self.df = df
        self.__reaction_column_name = reaction_column_name
        self.tokenizer = SmilesTokenizer()
        self.fragment_bond = fragment_bond
        self.augmented_columns: Set[str] = set()

    #
    # Private Methods
    #

    def __randomize_smiles(
        self, smiles: str, random_type: RandomType, permutations: int
    ) -> List[str]:
        """
        Randomizes a molecules SMILES string that might contain fragment bonds
        and returns a number of augmented versions of the SMILES equal to permutations.

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
                ".".join(
                    [
                        self.fragment_bond.join(
                            [
                                Augmenter.__randomize_smiles_without_fragment(
                                    fragment, random_type
                                )
                                for fragment in group.split(self.fragment_bond)
                            ]
                        )
                        for group in smiles.split(".")
                    ]
                )
            )
        return list_of_smiles

    #
    # Private Static Methods
    #

    @staticmethod
    def __randomize_smiles_without_fragment(
        smiles: str, random_type: RandomType
    ) -> str:
        """
        Generates a random version of a SMILES without a fragment bond

        Args:
            smiles (str): The pandas DataFrame to be split into training, validation, and test sets.
            random_type (RandomType): The type of randomization to be applied.

        Raises:
            InvalidSmiles: for invalid SMILES (raised via rxn.chemutils).
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
        raise ValueError(f"Invalid random type: {random_type}")

    @staticmethod
    def __randomize_molecules(smiles: str, permutations: int) -> List[str]:
        """
        Randomizes the order of the molecules inside a SMILES string that might
        contain fragment bonds and returns a number of augmented versions of the
        SMILES equal to permutations.

        For a number of molecules smaller than permutations, returns a number of
        permutations equal to the number of molecules

        Args:
            smiles (str): The molecules SMILES to augment
            permutations (int): The number of permutations to deliver for the SMILES

        Returns:
            List[str]: The list of randomized SMILES

        """
        # Raise for empty SMILES
        if not smiles:
            raise ValueError

        molecules_list = smiles.split(".")
        total_permutations = range(min(math.factorial(len(molecules_list)), 4000000))
        permutation_indices = random.sample(
            total_permutations, min(permutations, len(molecules_list))
        )
        permuted_molecules_smiles = []
        for idx in permutation_indices:
            permuted_precursors = molecules_permutation_given_index(molecules_list, idx)
            permuted_molecules_smiles.append(".".join(permuted_precursors))

        return permuted_molecules_smiles

    #
    # Public Methods
    #

    def augment(
        self,
        random_type: RandomType = RandomType.unrestricted,
        rxn_section_to_augment: ReactionSection = ReactionSection.precursors,
        permutations: int = 1,
    ) -> pd.DataFrame:
        """
        Creates samples for the augmentation. Returns a a pandas Series containing the
        augmented samples.

        Args:
            random_type (RandomType): The string identifying the type of randomization to apply.
                "molecules" for randomization of the molecules (canonical SMILES kept)
                "unrestricted" for unrestricted randomization
                "restricted" for restricted randomization
                "rotated" for rotated randomization
                For details on the differences:
                https://github.com/undeadpixel/reinvent-randomized and
                https://github.com/GLambard/SMILES-X
            rxn_section_to_augment (ReactionSection): The section of the rxn SMILES to augment.
                "precursors" for augmenting only the precursors
                "products" for augmenting only the products
            permutations (int): The number of permutations to generate for each SMILES

        Returns:
            pd.DataFrame: A pandas Series containing the augmented samples.
        """

        if rxn_section_to_augment is ReactionSection.precursors:
            self.df[f"precursors_{random_type.name}"] = self.df[
                self.__reaction_column_name
            ].apply(lambda smiles: smiles.replace(" ", "").split(">>")[0])
            if "products" not in self.df.keys():
                self.df["products"] = self.df[self.__reaction_column_name].apply(
                    lambda smiles: smiles.replace(" ", "").split(">>")[1]
                )
            columns_to_augment = [f"precursors_{random_type.name}"]
            columns_to_join = [f"precursors_{random_type.name}", "products"]

        elif rxn_section_to_augment is ReactionSection.products:
            self.df[f"products_{random_type.name}"] = self.df[
                self.__reaction_column_name
            ].apply(lambda smiles: smiles.replace(" ", "").split(">>")[1])
            if "precursors" not in self.df.keys():
                self.df["precursors"] = self.df[self.__reaction_column_name].apply(
                    lambda smiles: smiles.replace(" ", "").split(">>")[0]
                )
            columns_to_augment = [f"products_{random_type.name}"]
            columns_to_join = ["precursors", f"products_{random_type.name}"]
        else:
            raise ValueError(
                f"Invalid reaction section to augment: {rxn_section_to_augment.name}"
            )

        for column in columns_to_augment:
            if random_type != RandomType.molecules:
                self.df[column] = self.df[column].apply(
                    lambda smiles: self.__randomize_smiles(
                        smiles, random_type, permutations
                    )
                )
            else:
                self.df[column] = self.df[column].apply(
                    lambda smiles: self.__randomize_molecules(smiles, permutations)
                )

        # Exploding the dataframe columns where I have the list of augmented
        # versions of a SMILES (the list length is the number of permutations)
        self.df = (
            self.df.set_index(
                [col for col in self.df.keys() if col not in columns_to_augment]
            )
            .apply(pd.Series.explode)
            .reset_index()
        )

        augmented_column_name = f"rxn_{random_type.name}"
        self.augmented_columns.add(augmented_column_name)
        self.df[augmented_column_name] = self.df.apply(
            lambda x: ">>".join(x[columns_to_join]), axis=1
        )

        return self.df

    #
    # Public Static Methods
    #

    @staticmethod
    def read_csv(
        filepath: str, reaction_column_name: str, fragment_bond: str = "."
    ) -> "Augmenter":
        """A helper function to read a list or csv of SMILES.

        Args:
            filepath (str): The path to the text file containing the molecules SMILES.
            reaction_column_name: The name of the reaction column (or the name that wil be given
                to the reaction column if the input file has no headers).
            fragment_bond (str): The fragment token in the reaction SMILES

        Returns:
            Augmenter: A new augmenter instance.
        """
        df = pd.read_csv(filepath, lineterminator="\n")
        if len(df.columns) == 1:
            df.rename(columns={df.columns[0]: reaction_column_name}, inplace=True)

        return Augmenter(df, reaction_column_name, fragment_bond)


def augment(cfg: AugmentConfig) -> None:
    output_file_path = Path(cfg.output_file_path)
    if not Path(cfg.input_file_path).exists():
        raise ValueError(
            f"Input file for standardization does not exist: {cfg.input_file_path}"
        )

    # Create a instance of the Augmenter.
    ag = Augmenter.read_csv(
        cfg.input_file_path, cfg.reaction_column_name, cfg.fragment_bond.value
    )
    columns_to_keep = list(ag.df.columns)

    # Perform augmentation
    ag.augment(
        random_type=cfg.random_type,
        rxn_section_to_augment=cfg.rxn_section_to_augment,
        permutations=cfg.permutations,
    )

    columns_to_keep.extend(ag.augmented_columns)
    if not cfg.keep_intermediate_columns:
        ag.df = ag.df[columns_to_keep]

    # Exporting augmented samples
    ag.df.to_csv(output_file_path, index=False)
