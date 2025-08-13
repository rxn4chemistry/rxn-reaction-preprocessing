# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED
"""A utility class to apply standardization to the data"""

from pathlib import Path
from typing import List, Optional

from attr import define
from rxn.chemutils.miscellaneous import remove_chiral_centers
from rxn.chemutils.reaction_smiles import parse_any_reaction_smiles
from rxn.utilities.csv import CsvIterator, StreamingCsvEditor

from rxn.reaction_preprocessing.annotations.molecule_annotation import (
    MoleculeAnnotation,
    load_annotations_multiple,
)
from rxn.reaction_preprocessing.config import StandardizeConfig
from rxn.reaction_preprocessing.molecule_standardizer import MoleculeStandardizer


class Standardizer:
    def __init__(
        self,
        annotations: List[MoleculeAnnotation],
        discard_unannotated_metals: bool,
        reaction_column_name: str,
        fragment_bond: Optional[str] = None,
        remove_stereo_if_not_defined_in_precursors: bool = False,
        keep_intermediate_columns: bool = False,
    ):
        """Creates a new instance of the Standardizer class.

        Args:
            annotations: A list of MoleculeAnnotation objects used to perform the substitutions/rejections
            discard_unannotated_metals: whether reactions containing unannotated
                molecules with transition metals must be rejected.
            reaction_column_name: The name of the DataFrame column containing the reaction SMILES.
            fragment_bond: the fragment bond used in the dataframe.
            remove_stereo_if_not_defined_in_precursors: Remove chiral centers from products.
            keep_intermediate_columns: Whether the columns generated during preprocessing should be kept.
        """
        self.molecule_standardizer = MoleculeStandardizer(
            annotations=annotations,
            discard_missing_annotations=discard_unannotated_metals,
            canonicalize=True,
        )

        self.remove_stereo_if_not_defined_in_precursors = (
            remove_stereo_if_not_defined_in_precursors
        )
        self.fragment_bond = fragment_bond
        self.keep_intermediate_columns = keep_intermediate_columns

        self.rxn_column = reaction_column_name

    def standardize_file(self, input_csv: Path, output_csv: Path) -> None:
        """Standardize the reactions in a CSV file.

        Args:
            input_csv: CSV with the reactions to standardize.
            output_csv: CSV where to save the standardized reactions.
        """
        with open(input_csv, "rt") as f_in, open(output_csv, "wt") as f_out:
            csv_iterator = CsvIterator.from_stream(f_in)
            csv_iterator = self.standardize_iterator(csv_iterator)
            csv_iterator.to_stream(f_out)

    def standardize_iterator(self, csv_iterator: CsvIterator) -> CsvIterator:
        """Standardize the reactions in a CSV iterator.

        Same as ``standardize_file``, except that it acts directly on the iterator.

        Args:
            csv_iterator: input CSV iterator for the reactions to standardize.

        Returns:
            CsvIterator with reactions after the standardization step.
        """
        editor = self._instantiate_csv_editor()
        return editor.process(csv_iterator)

    def _instantiate_csv_editor(self) -> StreamingCsvEditor:
        if self.keep_intermediate_columns:
            return StreamingCsvEditor(
                columns_in=[self.rxn_column],
                columns_out=StandardizationOutput.keys(self.rxn_column),
                transformation=self._standardize_big,
            )
        else:
            return StreamingCsvEditor(
                columns_in=[self.rxn_column],
                columns_out=[self.rxn_column],
                transformation=self._standardize_small,
            )

    def standardize_one(self, rxn_smiles: str) -> "StandardizationOutput":
        original_rxn_smiles = rxn_smiles

        # Remove stereo information from products, if needed
        rxn_smiles = self._remove_stereo_if_not_defined_in_precursors(rxn_smiles)

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
        return StandardizationOutput(
            standardized_rxn_smiles=standardized_smiles,
            original_rxn_smiles=original_rxn_smiles,
            invalid_smiles=invalid_smiles,
            rejected_smiles=rejected_smiles,
            missing_annotations=missing_annotations,
        )

    def _standardize_small(self, rxn_smiles: str) -> str:
        return self.standardize_one(rxn_smiles).standardized_rxn_smiles

    def _standardize_big(self, rxn_smiles: str) -> List[str]:
        return self.standardize_one(rxn_smiles).values()

    def _remove_stereo_if_not_defined_in_precursors(self, rxn_smiles: str) -> str:
        """
        Remove stereocenters from products if not explainable by precursors.
        """
        if not self.remove_stereo_if_not_defined_in_precursors:
            return rxn_smiles

        reactants, reagents, products = rxn_smiles.split(">")
        if "@" in products and not ("@" in reactants or "@" in reagents):
            rxn_smiles = remove_chiral_centers(rxn_smiles)  # replaces with the group
        return rxn_smiles


@define
class StandardizationOutput:
    """Contains the results and additional information for the standardization
    of one reaction SMILES."""

    standardized_rxn_smiles: str
    original_rxn_smiles: str
    invalid_smiles: List[str]
    rejected_smiles: List[str]
    missing_annotations: List[str]

    @staticmethod
    def keys(rxn_column_name: str) -> List[str]:
        """Return the "keys" in the right order (same order as values)."""
        return [
            rxn_column_name,
            f"{rxn_column_name}_before_std",
            f"{rxn_column_name}_invalid_smiles",
            f"{rxn_column_name}_rejected_smiles",
            f"{rxn_column_name}_missing_annotations",
        ]

    def values(self) -> List[str]:
        """Return the "values" in the right order (same order as keys)."""
        return [
            self.standardized_rxn_smiles,
            self.original_rxn_smiles,
            str(self.invalid_smiles),
            str(self.rejected_smiles),
            str(self.missing_annotations),
        ]


def standardize(cfg: StandardizeConfig) -> None:
    output_path = Path(cfg.output_file_path)
    input_path = Path(cfg.input_file_path)
    if not input_path.exists():
        raise ValueError(
            f"Input file for standardization does not exist: {cfg.input_file_path}"
        )

    # Create a list of MoleculeAnnotations from the json files provided.
    annotations = load_annotations_multiple(cfg.annotation_file_paths)

    # Create an instance of the Standardizer
    std = Standardizer(
        annotations,
        discard_unannotated_metals=cfg.discard_unannotated_metals,
        reaction_column_name=cfg.reaction_column_name,
        fragment_bond=cfg.fragment_bond.value,
        remove_stereo_if_not_defined_in_precursors=cfg.remove_stereo_if_not_defined_in_precursors,
        keep_intermediate_columns=cfg.keep_intermediate_columns,
    )

    std.standardize_file(input_path, output_path)
