import logging
from typing import Callable, Optional, TextIO

from rxn.chemutils.reaction_equation import ReactionEquation
from rxn.chemutils.reaction_smiles import parse_any_reaction_smiles
from rxn.chemutils.utils import remove_atom_mapping
from rxn.utilities.csv import CsvIterator, StreamingCsvEditor
from rxn.utilities.files import PathLike

from rxn.reaction_preprocessing.config import InitialDataFormat, RxnImportConfig
from rxn.reaction_preprocessing.special_tokens import add_heat_token, add_light_token

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class RxnImportError(ValueError):
    """Exception for errors in the initial data import."""

    def __init__(self, msg: str):
        super().__init__(msg)


class InvalidColumn(RxnImportError):
    """Exception when a required column does not exist."""

    def __init__(self, column_name: str):
        super().__init__(f'No column named "{column_name}" in the input file.')


def _str2bool(v: str) -> bool:
    return v.lower() in ("yes", "true", "t", "1")


class RxnImporter:
    def __init__(
        self,
        data_format: InitialDataFormat,
        input_csv_column_name: str,
        reaction_column_name: str,
        fragment_bond: str,
        remove_atom_maps: bool,
        column_for_light: Optional[str],
        column_for_heat: Optional[str],
        keep_original_rxn_column: bool,
    ):
        self.data_format = data_format
        self.input_csv_column_name = input_csv_column_name
        self.rxn_column = reaction_column_name
        self.fragment_bond = fragment_bond
        self.remove_atom_maps = remove_atom_maps
        self.column_for_light = column_for_light
        self.column_for_heat = column_for_heat
        self.keep_original_rxn_column = keep_original_rxn_column

    def import_from_file(self, input_file: PathLike, output_csv: PathLike) -> None:
        """Function to import the reactions from one file and write them
        to another file.

        Goes through all the steps involved in the import: parsing of the
        reaction SMILES, removal of the atom mapping, addition of special tokens.
        """
        with open(input_file, "rt") as f_in, open(output_csv, "wt") as f_out:
            csv_iterator = self._load_into_iterator(f_in)
            csv_iterator = self.import_from_iterator(csv_iterator)
            csv_iterator.to_stream(f_out)

    def import_from_iterator(self, csv_iterator: CsvIterator) -> CsvIterator:
        """
        Same as ``import_rxns``, except that it acts on ``CsvIterator``
        instances instead of files.

        Args:
            csv_iterator: input CSV iterator for the reactions to import.

        Returns:
            CsvIterator with reactions after the import step.
        """
        csv_iterator = self._handle_original_rxn_column(csv_iterator)

        csv_iterator = self._parse_reaction_smiles(csv_iterator)

        # Add special tokens when necessary
        csv_iterator = self._maybe_add_light_token(csv_iterator)
        csv_iterator = self._maybe_add_heat_token(csv_iterator)

        csv_iterator = self._maybe_remove_atom_mapping(csv_iterator)

        return csv_iterator

    def _load_into_iterator(self, input_stream: TextIO) -> CsvIterator:
        """
        Load the initial dataframe.

        This functions (and the ones it calls) make sure to rename the column
        containing the data if it would be overwritten when transforming. The
        column with the data ends up with the name `original_column_name`.
        """

        if self.data_format is InitialDataFormat.TXT:
            return self._load_from_txt(input_stream)
        if self.data_format is InitialDataFormat.CSV:
            return self._load_from_csv(input_stream, separator=",")
        if self.data_format is InitialDataFormat.TSV:
            return self._load_from_csv(input_stream, separator="\t")
        raise ValueError(f"Unsupported data type: {self.data_format}")

    def _load_from_txt(self, input_stream: TextIO) -> CsvIterator:
        # We just get the iterator, and associate the final column name to it
        return CsvIterator(
            columns=[self.rxn_column],
            rows=([line.rstrip("\r\n")] for line in input_stream),
        )

    def _load_from_csv(self, input_stream: TextIO, separator: str) -> CsvIterator:
        csv_iterator = CsvIterator.from_stream(input_stream, delimiter=separator)

        if self.input_csv_column_name not in csv_iterator.columns:
            raise InvalidColumn(self.input_csv_column_name)

        # Rename the input column, if needed
        csv_iterator.columns = [
            self.rxn_column if c == self.input_csv_column_name else c
            for c in csv_iterator.columns
        ]

        return csv_iterator

    def _column_name_to_store_original_rxn(self) -> str:
        """Name of the column where the original data will end up."""

        # For txt: always with "_original" postfix
        if self.data_format is InitialDataFormat.TXT:
            return f"{self.rxn_column}_original"

        # For csv: only add "_original" postfix if the input and output columns are identical
        if self.input_csv_column_name == self.rxn_column:
            return f"{self.rxn_column}_original"

        return self.input_csv_column_name

    def _handle_original_rxn_column(self, csv_iterator: CsvIterator) -> CsvIterator:
        """Make sure that the original reaction SMILES is kept (if required)."""
        if not self.keep_original_rxn_column:
            return csv_iterator

        # To keep the original column: we apply a fake transformation that
        # just adds the column
        def identity(v: str) -> str:
            return v

        editor = StreamingCsvEditor(
            [self.rxn_column], [self._column_name_to_store_original_rxn()], identity
        )
        return editor.process(csv_iterator)

    def _parse_reaction_smiles(self, csv_iterator: CsvIterator) -> CsvIterator:
        editor = StreamingCsvEditor(
            [self.rxn_column],
            [self.rxn_column],
            self._reformat_smiles,
        )
        csv_iterator = editor.process(csv_iterator)
        return self._remove_invalid(csv_iterator)

    def _reformat_smiles(self, reaction_smiles: str) -> str:
        """Import a reaction SMILES in any format and convert it to an "IBM" RXN
        SMILES with the specified fragment bond.

        An empty string is returned if the reaction SMILES is not valid.
        """

        try:
            reaction = parse_any_reaction_smiles(reaction_smiles)
            return reaction.to_string(self.fragment_bond)
        except Exception:
            logger.info(f"Invalid reaction: {reaction_smiles}")
            return ""

    def _remove_invalid(self, csv_iterator: CsvIterator) -> CsvIterator:
        rxn_idx = csv_iterator.column_index(self.rxn_column)

        # Filter out the ones that have an empty SMILES string (note: is empty
        # if the import was unsuccessful).
        return CsvIterator(
            csv_iterator.columns, (row for row in csv_iterator.rows if row[rxn_idx])
        )

    def _maybe_add_special_token(
        self,
        csv_iterator: CsvIterator,
        add_special_token_fn: Callable[[ReactionEquation], ReactionEquation],
        special_token_column: Optional[str],
    ) -> CsvIterator:
        """
        Common function for adding light or heat tokens.

        Args:
            csv_iterator: iterator where to update the reaction SMILES.
            add_special_token_fn: function to call on a ReactionEquation to add the required token.
            special_token_column: column to look at in order to know when to add the special tokens.
        """

        # Do nothing if no column was specified for the special token
        if special_token_column is None:
            return csv_iterator

        if special_token_column not in csv_iterator.columns:
            raise InvalidColumn(special_token_column)

        # "partial" in order to make the function compatible with pandas.apply().
        def fn(reaction_smiles: str, special_token: str) -> str:
            return self._add_token(reaction_smiles, special_token, add_special_token_fn)

        editor = StreamingCsvEditor(
            [self.rxn_column, special_token_column],
            [self.rxn_column],
            fn,
        )
        return editor.process(csv_iterator)

    def _maybe_add_heat_token(self, csv_iterator: CsvIterator) -> CsvIterator:
        """Add the heat token to the precursors of the reaction SMILES when necessary."""

        return self._maybe_add_special_token(
            csv_iterator,
            add_special_token_fn=add_heat_token,
            special_token_column=self.column_for_heat,
        )

    def _maybe_add_light_token(self, csv_iterator: CsvIterator) -> CsvIterator:
        """Add the light token to the precursors of the reaction SMILES when necessary."""

        return self._maybe_add_special_token(
            csv_iterator,
            add_special_token_fn=add_light_token,
            special_token_column=self.column_for_light,
        )

    def _maybe_remove_atom_mapping(self, csv_iterator: CsvIterator) -> CsvIterator:
        """Remove the atom mapping if required by the config.

        NB: This will not clean up the SMILES, i.e. "[CH3:1][CH3:2]" is converted
        to "[CH3][CH3]" and not to "CC". The standardization step will do that.
        """

        if not self.remove_atom_maps:
            return csv_iterator

        editor = StreamingCsvEditor(
            [self.rxn_column],
            [self.rxn_column],
            remove_atom_mapping,
        )
        return editor.process(csv_iterator)

    def _add_token(
        self,
        reaction_smiles: str,
        special_token: str,
        add_special_token_fn: Callable[[ReactionEquation], ReactionEquation],
    ) -> str:
        """Function to use with StreamingCsvEditor to update the reaction SMILES."""

        special_flag_active = _str2bool(special_token)

        if not special_flag_active:
            # do nothing if the reaction is not run under light / heat / etc.
            return reaction_smiles

        reaction = parse_any_reaction_smiles(reaction_smiles)
        reaction = add_special_token_fn(reaction)

        return reaction.to_string(self.fragment_bond)


def rxn_import(cfg: RxnImportConfig) -> None:
    """
    Initial import of reaction data, as a first step of the reaction preprocessing.
    """
    importer = RxnImporter(
        data_format=cfg.data_format,
        input_csv_column_name=cfg.input_csv_column_name,
        reaction_column_name=cfg.reaction_column_name,
        fragment_bond=cfg.fragment_bond.value,
        remove_atom_maps=cfg.remove_atom_mapping,
        column_for_light=cfg.column_for_light,
        column_for_heat=cfg.column_for_heat,
        keep_original_rxn_column=cfg.keep_original_rxn_column,
    )

    importer.import_from_file(input_file=cfg.input_file, output_csv=cfg.output_csv)
