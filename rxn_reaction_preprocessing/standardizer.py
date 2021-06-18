# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED
""" A utility class to apply standardization to the data """
from pathlib import Path
from typing import List
from typing import Optional

import pandas as pd
from rdkit import RDLogger
from rxn_chemutils.miscellaneous import is_valid_smiles
from rxn_chemutils.reaction_equation import canonicalize_compounds
from rxn_chemutils.reaction_equation import ReactionEquation

from rxn_reaction_preprocessing.annotations.missing_annotation_detector import MissingAnnotationDetector
from rxn_reaction_preprocessing.annotations.molecule_annotation import load_annotations_multiple
from rxn_reaction_preprocessing.annotations.molecule_annotation import MoleculeAnnotation
from rxn_reaction_preprocessing.annotations.molecule_replacer import MoleculeReplacer
from rxn_reaction_preprocessing.annotations.rejected_molecules_filter import RejectedMoleculesFilter
from rxn_reaction_preprocessing.config import StandardizeConfig

RDLogger.DisableLog('rdApp.*')


class Standardizer:

    def __init__(
        self,
        df: pd.DataFrame,
        annotations: List[MoleculeAnnotation],
        discard_unannotated_metals: bool,
        reaction_column_name: str,
        fragment_bond: Optional[str] = None,
    ):
        """Creates a new instance of the Standardizer class.

        Args:
            df: A pandas DataFrame containing the reaction SMILES.
            annotations: A list of MoleculeAnnotation objects used to perform the substitutions/rejections
            discard_unannotated_metals: whether reactions containing unannotated
                molecules with transition metals must be rejected.
            reaction_column_name: The name of the DataFrame column containing the reaction SMILES.
            fragment_bond: the fragment bond used in the dataframe.
        """
        self.df = df
        self.annotations = annotations
        self.discard_unannotated_metals = discard_unannotated_metals
        self.missing_annotation_detector = MissingAnnotationDetector.from_molecule_annotations(
            self.annotations
        )
        self.molecule_filter = RejectedMoleculesFilter.from_molecule_annotations(self.annotations)
        self.molecule_replacer = MoleculeReplacer.from_molecule_annotations(self.annotations)
        self.__reaction_column_name = reaction_column_name
        self.fragment_bond = fragment_bond

    def __detect_missing_annotations(self) -> None:
        """
        Runs over the df and checks if some reactions still need some annotations.
        Replaces the reaction with '>>' if the annotation is needed and adds a column 'rxn_needed_annotations'
        """
        self.df['rxn_needed_annotations'] = self.df[self.__reaction_column_name].apply(
            lambda x: list(
                self.missing_annotation_detector.
                missing_in_reaction_smiles(x, fragment_bond=self.fragment_bond)
            )
        )
        if self.discard_unannotated_metals:
            self.df[self.__reaction_column_name] = self.df.apply(
                lambda x: x[self.__reaction_column_name]
                if not x['rxn_needed_annotations'] else '>>',
                axis=1
            )

    def __filter_reactions_with_rejected_molecules(self) -> None:
        """
        Runs over the df and checks if some reactions contain rejected molecules.
        Replaces the reaction with '>>' if it needs to be discarded
        """
        self.df['rxn_contains_rejected_molecules'] = self.df[self.__reaction_column_name].apply(
            lambda x: not self.molecule_filter.
            is_valid_reaction_smiles(x, fragment_bond=self.fragment_bond)
        )
        self.df[self.__reaction_column_name] = self.df.apply(
            lambda x: x[self.__reaction_column_name]
            if not x['rxn_contains_rejected_molecules'] else '>>',
            axis=1
        )

    def __replace_molecules_in_reactions(self) -> None:
        """
        Runs over the df and for each reaction replaces the molecules that have an annotation.
        """
        self.df[self.__reaction_column_name] = self.df[self.__reaction_column_name].apply(
            lambda x: self.molecule_replacer.
            replace_in_reaction_smiles(x, fragment_bond=self.fragment_bond)
        )

    def standardize(self, canonicalize: bool = True):
        """
         Standardizes the entries of self.df[self.__reaction_column_name]
        """
        self.df[f'{self.__reaction_column_name}_before_std'] = self.df[self.__reaction_column_name]
        self.df[self.__reaction_column_name
                ] = self.df[self.__reaction_column_name
                            ].apply(lambda x: self.__validate_mild(x, canonicalize=canonicalize))
        self.__detect_missing_annotations()
        self.__filter_reactions_with_rejected_molecules()
        self.__replace_molecules_in_reactions()
        return self

    def __validate_mild(self, smiles: str, canonicalize: bool) -> str:
        """
        Checks if the input reaction SMILES is valid. Returns the canonical SMILES if it is.
        If it is not it returns '>>'.

        Args:
            smiles (str): a reaction SMILES
        Returns:
            str: the canonical version of the reaction SMILES
        """
        reaction_equation = ReactionEquation.from_string(smiles, self.fragment_bond)
        if not all(is_valid_smiles(molecule) for group in reaction_equation for molecule in group):
            return '>>'

        if not canonicalize:
            return smiles

        return canonicalize_compounds(reaction_equation).to_string(
            fragment_bond=self.fragment_bond
        )

    @staticmethod
    def read_csv(
        filepath: str,
        annotations: List[MoleculeAnnotation],
        discard_unannotated_metals: bool,
        reaction_column_name: str,
        fragment_bond: str = None
    ):
        """
        A helper function to read a list or csv of VALID reactions (in the sense of RDKIT).

        Args:
            filepath (str): The path to the text file containing the reactions.
            annotations: A list of MoleculeAnnotation objects used to perform the substitutions/rejections
            discard_unannotated_metals: whether reactions containing unannotated
                molecules with transition metals must be rejected.
            reaction_column_name: The name of the reaction column (or the name that wil be given to the reaction
                column if the input file has no headers)
            fragment_bond: the fragment bond used.

        Returns:
            : A new standardizer instance.
        """
        df = pd.read_csv(filepath, lineterminator='\n')
        if len(df.columns) == 1:
            df.rename(columns={df.columns[0]: reaction_column_name}, inplace=True)

        return Standardizer(
            df,
            annotations=annotations,
            discard_unannotated_metals=discard_unannotated_metals,
            reaction_column_name=reaction_column_name,
            fragment_bond=fragment_bond,
        )


def standardize(cfg: StandardizeConfig) -> None:
    output_file_path = Path(cfg.output_file_path)
    if not Path(cfg.input_file_path).exists():
        raise ValueError(f'Input file for standardization does not exist: {cfg.input_file_path}')

    # Create a list of MoleculeAnnotations from the json files provided.
    annotations = load_annotations_multiple(cfg.annotation_file_paths)

    # Create an instance of the Standardizer
    std = Standardizer.read_csv(
        cfg.input_file_path,
        annotations,
        discard_unannotated_metals=cfg.discard_unannotated_metals,
        reaction_column_name=cfg.reaction_column_name,
        fragment_bond=cfg.fragment_bond.value
    )

    # Perform standardization
    std.standardize(canonicalize=True)

    # Exporting standardized samples
    std.df.to_csv(output_file_path)
