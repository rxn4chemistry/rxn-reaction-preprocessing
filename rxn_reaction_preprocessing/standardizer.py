# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED

""" A utility class to apply standardization to the data """
import json
import re
from typing import List, Optional
from typing import Pattern
from typing import Tuple

import pandas as pd
from rdkit import RDLogger
from rxn_chemutils.miscellaneous import is_valid_smiles
from rxn_chemutils.reaction_equation import ReactionEquation

RDLogger.DisableLog('rdApp.*')


class Patterns:

    def __init__(self, jsonpath: str, fragment_bond: str = None):
        """Creates a new instance of the Patterns class.

        Args:
            jsonpath (str): path to json file containing the molecule SMILES strings to be found in the format
            Dict[str,List[List[str]]]
            and the substitution to replace in the reaction SMILES
            fragment_bond (str): the fragment bond used in the patterns.
        """
        with open(jsonpath) as j:
            self.patterns = json.load(j)
        self.fragment_bond = fragment_bond
        self.__escape_special_chars_and_isolate()
        self.compiled_patterns: List[Tuple[Pattern[str], str]] = []

    def __escape_special_chars_and_isolate(self) -> None:
        """
        Modifies the patterns stored in the object by escaping special chars in the molecule smiles strings to find and
        appending/prepending the isolation regex pattern. Isolation = the patterns are applied to entire molecules
        or entire fragments of molecules
        """
        special_chars = ['\\', '.', '+', '*', '?', '^', '$', '(', ')', '[', ']', '}', '{', '|']

        for key, elem in self.patterns.items():
            for i in range(len(elem)):
                for char in special_chars:
                    self.patterns[key][i][0] = self.patterns[key][i][0].replace(
                        char,
                        r'\{}'.format(char)  # rf'\{char}' fails on my system (dpr)
                    )
                self.patterns[key][i][
                    0] = r'(?:^|(?<=\~|\.|>))' + self.patterns[key][i][0] + r'(?=\~|\.|>|$)'

    def __compile_patterns(self) -> List[Tuple[Pattern[str], str]]:
        """
        Generates a list of Tuples with the precompiled patterns and the substitution SMILES string.

        Returns:
            List[Tuple[Pattern[str], str]]: a list of Tuples with the precompiled patterns and the substitution SMILES
             string
        """
        compiled_patterns = []
        for elem in self.patterns.values():
            for pat in elem:
                regex = re.compile(pat[0])
                compiled_patterns.append((regex, pat[1]))
        return compiled_patterns

    def compile_patterns(self):
        """
        public method, see self.__compile_patterns()
        """
        return self.__compile_patterns()


class Standardizer:

    def __init__(
        self,
        df: pd.DataFrame,
        patterns: Patterns,
        reaction_column_name: str,
        fragment_bond: Optional[str] = None
    ):
        """Creates a new instance of the Standardizer class.

        Args:
            df (pd.DataFrame): A pandas DataFrame containing the reaction SMILES.
            patterns (Patterns): An instance of the Patterns class, containing a dictionary of the patterns and the
            substitutions to be found and replaced in the reaction SMILES as well as the fragment used in the patterns.
            reaction_column_name (str): The name of the DataFrame column containing the reaction SMARTS.
            fragment_bond (str): the fragment bond used.
        """
        self.df = df
        self.patterns = patterns
        self.__reaction_column_name = reaction_column_name
        self.fragment_bond = fragment_bond
        self.current_smiles = ''

        # Dealing with the possible mismatch between the fragment-bond token in the patterns and in the provided data
        if self.fragment_bond and self.fragment_bond != '.' and self.patterns.fragment_bond and \
                self.fragment_bond != self.patterns.fragment_bond:
            for key, elem in self.patterns.patterns.items():
                for i in range(len(elem)):
                    self.patterns.patterns[key][i][0] = self.patterns.patterns[key][i][0]\
                        .replace(self.patterns.fragment_bond, self.fragment_bond)
                    self.patterns.patterns[key][i][1] = self.patterns.patterns[key][i][1]\
                        .replace(self.patterns.fragment_bond, self.fragment_bond)
            self.patterns.fragment_bond = self.fragment_bond
            self.patterns.compiled_patterns = self.patterns.compile_patterns()
        else:
            self.patterns.compiled_patterns = self.patterns.compile_patterns()
            self.fragment_bond = self.patterns.fragment_bond

    def __standardize_reaction_smiles(self, smiles: str):
        """
        Standardizes a single reaction SMILES

        Args:
            smiles (str): a reaction SMILES
        """
        self.current_smiles = smiles
        new_smiles = smiles
        for i in range(len(self.patterns.compiled_patterns)):
            regex = self.patterns.compiled_patterns[i][0]
            new_smiles = regex.sub(self.patterns.compiled_patterns[i][1], new_smiles)

        return self.__validate_mild(new_smiles)

    def standardize(self):
        """
         Standardizes the entries of self.df[_valid_column]
        """
        self.df.rename(
            columns={self.__reaction_column_name: f'{self.__reaction_column_name}_before_std'},
            inplace=True
        )
        self.df[f'{self.__reaction_column_name}'] = self.df[f'{self.__reaction_column_name}_before_std'].\
            apply(lambda x: self.__standardize_reaction_smiles(x))
        return self

    def __validate_mild(self, smiles: str) -> str:
        """
        Checks if the input reaction SMILES is valid. Returns the input SMILES if it is.
        If it is not it returns '>>'.

        Args:
            smiles (str): a reaction SMILES
        Returns:
            str: the canonical version of the reaction SMILES
        """
        reaction_equation = ReactionEquation.from_string(smiles, self.fragment_bond)
        if all(is_valid_smiles(molecule) for group in reaction_equation for molecule in group):
            return smiles
        else:
            return '>>'

    @staticmethod
    def read_csv(
        filepath: str, patterns: Patterns, reaction_column_name: str, fragment_bond: str = None
    ):
        """
        A helper function to read a list or csv of VALID reactions (in the sense of RDKIT).

        Args:
            filepath (str): The path to the text file containing the reactions.
            patterns (Patterns): An instance of the Patterns class, containing a dictionary of the patterns and the
             substitutions to be found and replaced in the reaction SMILES as well as the fragment used in the patterns.
            reaction_column_name (str): The name of the reaction column (or the name that wil be given to the reaction
            column if the input file has no headers).
            fragment_bond (str): the fragment bond used.

        Returns:
            : A new preprocessor instance.
        """
        df = pd.read_csv(filepath, lineterminator='\n')
        if len(df.columns) == 1:
            df.rename(columns={df.columns[0]: reaction_column_name}, inplace=True)

        return Standardizer(
            df,
            patterns=patterns,
            reaction_column_name=reaction_column_name,
            fragment_bond=fragment_bond
        )
