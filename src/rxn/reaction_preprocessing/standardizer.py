# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED
""" A utility class to apply standardization to the data """
from pathlib import Path
from typing import List, Optional

import pandas as pd
from rdkit import RDLogger
from rxn.chemutils.miscellaneous import remove_chiral_centers
from rxn.chemutils.reaction_smiles import parse_any_reaction_smiles

from rxn.reaction_preprocessing.annotations.molecule_annotation import (
    MoleculeAnnotation,
    load_annotations_multiple,
)
from rxn.reaction_preprocessing.config import StandardizeConfig
from rxn.reaction_preprocessing.molecule_standardizer import MoleculeStandardizer

RDLogger.DisableLog("rdApp.*")


class Standardizer:
    def __init__(
        self,
        df: pd.DataFrame,
        annotations: List[MoleculeAnnotation],
        discard_unannotated_metals: bool,
        reaction_column_name: str,
        fragment_bond: Optional[str] = None,
        remove_stereo_if_not_defined_in_precursors: bool = False,
    ):
        """Creates a new instance of the Standardizer class.

        Args:
            df: A pandas DataFrame containing the reaction SMILES.
            annotations: A list of MoleculeAnnotation objects used to perform the substitutions/rejections
            discard_unannotated_metals: whether reactions containing unannotated
                molecules with transition metals must be rejected.
            reaction_column_name: The name of the DataFrame column containing the reaction SMILES.
            fragment_bond: the fragment bond used in the dataframe.
            remove_stereo_if_not_defined_in_precursors: Remove chiral centers from products.
        """
        self.df = df
        self.molecule_standardizer = MoleculeStandardizer(
            annotations=annotations,
            discard_missing_annotations=discard_unannotated_metals,
            canonicalize=True,
        )
        self.fragment_bond = fragment_bond
        self.remove_stereo_if_not_defined_in_precursors = (
            remove_stereo_if_not_defined_in_precursors
        )

        self.rxn_column = reaction_column_name
        self.rxn_before_std_column = f"{self.rxn_column}_before_std"
        self.invalid_smiles_column = f"{self.rxn_column}_invalid_smiles"
        self.rejected_smiles_column = f"{self.rxn_column}_rejected_smiles"
        self.missing_annotations_column = f"{self.rxn_column}_missing_annotations"

    def __remove_stereo_if_not_defined_in_precursors(self, rxn_smiles: str) -> str:
        """
        Remove stereocenters from products if not explainable by precursors.
        """
        if not self.remove_stereo_if_not_defined_in_precursors:
            return rxn_smiles

        reactants, reagents, products = rxn_smiles.split(">")
        if "@" in products and not ("@" in reactants or "@" in reagents):
            rxn_smiles = remove_chiral_centers(rxn_smiles)  # replaces with the group
        return rxn_smiles

    def standardize(self, canonicalize: bool = True) -> "Standardizer":
        """
        Standardizes the entries of self.df[self.__reaction_column_name]
        """
        self.molecule_standardizer.canonicalize = canonicalize

        # Make a copy of the non-standardized reaction SMILES. Achieved by
        # renaming to enable the "join" operation below without conflict.
        self.df.rename(
            columns={self.rxn_column: self.rxn_before_std_column}, inplace=True
        )

        new_columns: pd.DataFrame = self.df.apply(self.process_row, axis=1)
        new_columns.columns = [
            self.rxn_column,
            self.invalid_smiles_column,
            self.rejected_smiles_column,
            self.missing_annotations_column,
        ]
        # Merge the new columns
        self.df = self.df.join(new_columns)

        return self

    def process_row(self, x: pd.Series) -> pd.Series:
        """
        Function applied to every row of the dataframe to get the new columns.

        Returns:
            Pandas Series with 1) the standardized reaction SMILES, 2) the list
                of invalid molecules, 3) the list of rejected molecules (from
                the annotations), 4) the list of missing annotations.
        """
        # Get RXN SMILES from the column
        rxn_smiles = x[self.rxn_before_std_column]

        # Remove stereo information from products, if needed
        rxn_smiles = self.__remove_stereo_if_not_defined_in_precursors(rxn_smiles)

        # Read the reaction SMILES while allowing for different formats (with
        # fragment bond, extended reaction SMILES, etc.).
        reaction_equation = parse_any_reaction_smiles(rxn_smiles)

        (
            standardized_reaction,
            invalid_smiles,
            rejected_smiles,
            missing_annotations,
        ) = self.molecule_standardizer.standardize_in_equation_with_errors(
            reaction_equation, propagate_exceptions=False
        )

        standardized_smiles = standardized_reaction.to_string(self.fragment_bond)
        return pd.Series(
            [standardized_smiles, invalid_smiles, rejected_smiles, missing_annotations]
        )

    @staticmethod
    def read_csv(
        filepath: str,
        annotations: List[MoleculeAnnotation],
        discard_unannotated_metals: bool,
        reaction_column_name: str,
        fragment_bond: Optional[str] = None,
        remove_stereo_if_not_defined_in_precursors: bool = False,
    ) -> "Standardizer":
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
            remove_stereo_if_not_defined_in_precursors: Remove chiral centers from products.

        Returns:
            : A new standardizer instance.
        """
        df = pd.read_csv(filepath, lineterminator="\n")
        if len(df.columns) == 1:
            df.rename(columns={df.columns[0]: reaction_column_name}, inplace=True)

        return Standardizer(
            df,
            annotations=annotations,
            discard_unannotated_metals=discard_unannotated_metals,
            reaction_column_name=reaction_column_name,
            fragment_bond=fragment_bond,
            remove_stereo_if_not_defined_in_precursors=remove_stereo_if_not_defined_in_precursors,
        )


def standardize(cfg: StandardizeConfig) -> None:
    output_file_path = Path(cfg.output_file_path)
    if not Path(cfg.input_file_path).exists():
        raise ValueError(
            f"Input file for standardization does not exist: {cfg.input_file_path}"
        )

    # Create a list of MoleculeAnnotations from the json files provided.
    annotations = load_annotations_multiple(cfg.annotation_file_paths)

    # Create an instance of the Standardizer
    std = Standardizer.read_csv(
        cfg.input_file_path,
        annotations,
        discard_unannotated_metals=cfg.discard_unannotated_metals,
        reaction_column_name=cfg.reaction_column_name,
        fragment_bond=cfg.fragment_bond.value,
        remove_stereo_if_not_defined_in_precursors=cfg.remove_stereo_if_not_defined_in_precursors,
    )

    columns_to_keep = list(std.df.columns)

    # Perform standardization
    std.standardize(canonicalize=True)

    if not cfg.keep_intermediate_columns:
        std.df = std.df[columns_to_keep]

    # Exporting standardized samples
    std.df.to_csv(output_file_path, index=False)
