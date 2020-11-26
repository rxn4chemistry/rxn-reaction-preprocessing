""" A utility class to apply standardization to the data """
import re
import json
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
from typing import List, Tuple, Pattern, Optional
from data_preprocessor.utils import InvalidSmiles

RDLogger.DisableLog('rdApp.*')

class Patterns:
    def __init__(self, jsonpath: str,
                 fragment_bond: str = None):
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
        self.compiled_patterns: Optional[List[Tuple[Pattern[str], str]],None] = None

    def __escape_special_chars_and_isolate(self) -> None:
        """
        Modifies the patterns stored in the object by escaping special chars in the molecule smiles strings to find and
        appending/prepending the isolation regex pattern. Isolation = the patterns are applied to entire molecules
        or entire fragments of molecules
        """
        special_chars = ["\\" ,".", "+", "*", "?", "^", "$", "(", ")", "[", "]", "}", "{", "|"]

        for key, elem in self.patterns.items():
            for i in range(len(elem)):
                for char in special_chars:
                    self.patterns[key][i][0] = self.patterns[key][i][0].replace(char, f"\{char}")
                self.patterns[key][i][0] = '(?:^|(?<=\~|\.|>))' + self.patterns[key][i][0] + '(?=\~|\.|>|$)'

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
            valid_column: str = '_rxn_valid',
            fragment_bond: str = None
    ):
        """Creates a new instance of the Standardizer class.

        Args:
            df (pd.DataFrame): A pandas DataFrame containing the reaction SMILES.
            patterns (Patterns): An instance of the Patterns class, containing a dictionary of the patterns and the
            substitutions to be found and replaced in the reaction SMILES as well as the fragment used in the patterns.
            valid_column (str): column on which to apply a standardization. Default _rxn_valid
            fragment_bond (str): the fragment bond used.
        """
        self.df = df
        self.patterns = patterns
        self.__valid_column = valid_column
        self.fragment_bond = fragment_bond
        self.current_smiles = ''

        # Check if the input reaction SMILES are (RDKIT) valid
        self.df[self.__valid_column].apply(lambda x: self.__check_correct_chemistry(x))

        if (self.fragment_bond and self.patterns.fragment_bond) and (self.fragment_bond != self.patterns.fragment_bond):
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

        return self.__check_correct_chemistry(new_smiles)

    def standardize(self):
        """
         Standardizes the entries of self.df[_valid_column]
        """
        self.df[f"std_{self.__valid_column}"] = self.df[self.__valid_column].\
            apply(lambda x: self.__standardize_reaction_smiles(x))
        print(len(self.df.loc[self.df[f"std_{self.__valid_column}"] != self.df[f"{self.__valid_column}"]]))
        return self

    def __check_correct_chemistry(self, smiles: str):
        """
        Check on the standardization, as the input reaction SMILES are assumed to be valid, then the stardardization
        should not change this by invalidating. If a smile is invalid an error is raised

        Args:
            smiles (str): a reaction SMILES
        """
        try:
            if self.fragment_bond:
                reactants, products = smiles.replace(self.fragment_bond, ".").split('>>')
                Chem.MolToSmiles(Chem.MolFromSmiles(reactants)), Chem.MolToSmiles(Chem.MolFromSmiles(products))
            else:
                reactants, products = smiles.split('>>')
                Chem.MolToSmiles(Chem.MolFromSmiles(reactants)), Chem.MolToSmiles(Chem.MolFromSmiles(products))
        except InvalidSmiles:
            print(f'Something went wrong with this standardization: {smiles} -> {self.current_smiles}')
            raise
        return smiles

    @staticmethod
    def read_csv(filepath: str, patterns: Patterns, valid_column: str = '_rxn_valid', fragment_bond: str = None,
                 kwargs={}):
        """
        A helper function to read a list or csv of VALID reactions (in the sense of RDKIT).

        Args:
            filepath (str): The path to the text file containing the reactions.
            patterns (Patterns): An instance of the Patterns class, containing a dictionary of the patterns and the
             substitutions to be found and replaced in the reaction SMILES as well as the fragment used in the patterns.
            valid_column (str): column on which to apply a standardization. Default _rxn_valid
            fragment_bond (str): the fragment bond used.
            kwargs (dict, optional): Additional arguments to supply to the internal Pandas read_csv method.
             Defaults to {}.

        Returns:
            : A new preprocessor instance.
        """
        df = pd.read_csv(filepath, **kwargs)
        if len(df.columns) == 1:
            df.rename(columns={df.columns[0]: valid_column}, inplace=True)

        return Standardizer(df, patterns, valid_column, fragment_bond)
