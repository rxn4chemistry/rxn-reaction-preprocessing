""" A utility class to augment the dataset files """
import math
import random
from enum import Enum
from typing import List

import pandas as pd
from rdkit import Chem

from rxn_reaction_preprocessing.smiles_tokenizer import SmilesTokenizer
from rxn_reaction_preprocessing.utils import InvalidSmiles


class RandomType(Enum):
    molecules = 1
    unrestricted = 2
    restricted = 3
    rotated = 4


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

    def __init__(self, df: pd.DataFrame, fragment_bond: str = None):
        """Creates a new instance of the Augmenter class.

        Args:
            df (pd.DataFrame): A pandas DataFrame containing the molecules SMILES.
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
            raise InvalidSmiles(smiles)

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

        Returns:
            str: the randomized SMILES
        """
        mol = Chem.MolFromSmiles(smiles)

        if not mol:
            raise InvalidSmiles(smiles)

        if random_type == RandomType.unrestricted:
            return Chem.MolToSmiles(mol, canonical=False, doRandom=True, isomericSmiles=False)
        elif random_type == RandomType.restricted:
            new_atom_order = list(range(mol.GetNumAtoms()))
            random.shuffle(new_atom_order)
            random_mol = Chem.RenumberAtoms(mol, newOrder=new_atom_order)
            return Chem.MolToSmiles(random_mol, canonical=False, isomericSmiles=False)
        elif random_type == RandomType.rotated:
            n_atoms = mol.GetNumAtoms()
            rotation_index = random.randint(0, n_atoms - 1)
            atoms = list(range(n_atoms))
            new_atoms_order = (
                atoms[rotation_index % len(atoms):] + atoms[:rotation_index % len(atoms)]
            )
            rotated_mol = Chem.RenumberAtoms(mol, new_atoms_order)
            return Chem.MolToSmiles(rotated_mol, canonical=False, isomericSmiles=False)
        return ''

    @staticmethod
    def __randomize_molecules(smiles: str, permutations: int) -> List[str]:
        """
        Randomizes the order of the molecules inside a SMILES string that might contain fragment bonds and returns
        a number of augmented versions of the SMILES equal to permutations.

        Args:
            smiles (str): The molecules SMILES to augment
            permutations (int): The number of permutations to deliver for the SMILES

        Returns:
            List[str]: The list of randomized SMILES

        """
        # Raise for empty SMILES
        if not smiles:
            raise InvalidSmiles(smiles)

        molecules_list = smiles.split('.')
        total_permutations = range(math.factorial(len(molecules_list)))
        permutation_indices = random.sample(total_permutations, permutations)
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
            fragment_bond (str): The token representing a fragment bond
            tokenize (str): Whether to return the tokenized version of the SMILES

        Returns:
            pd.DataFrame: A pandas Series containing the augmented samples.
        """

        self.df[f'{random_type.name}'] = self.df.apply(lambda smiles: smiles.replace(' ', ''))

        if random_type != 'molecules':
            self.df[f'{random_type.name}'] = self.df.apply(
                lambda smiles: self.__randomize_smiles(smiles, random_type, permutations)
            )
            self.df[f'{random_type.name}'] = self.df.apply(
                lambda smiles: Augmenter.__randomize_molecules(smiles, permutations)
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
    def read_csv(filepath: str, kwargs={}):
        """A helper function to read a list or csv of SMILES.

        Args:
            filepath (str): The path to the text file containing the molecules SMILES.
            kwargs (dict, optional): Additional arguments to supply to the internal Pandas read_csv method. Defaults to {}.

        Returns:
            Augmenter: A new augmenter instance.
        """
        df = pd.read_csv(filepath, **kwargs)
        if len(df.columns) == 1:
            df.rename(columns={df.columns[0]: 'smiles'}, inplace=True)

        return Augmenter(df)
