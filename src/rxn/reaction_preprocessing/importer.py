import logging
from functools import partial
from typing import Callable, Optional

import pandas as pd
from rxn.chemutils.reaction_equation import ReactionEquation
from rxn.chemutils.reaction_smiles import parse_any_reaction_smiles
from rxn.chemutils.utils import remove_atom_mapping

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

    def __init__(self, column_name: str, file_name: str):
        super().__init__(f'No column named "{column_name}" in the file "{file_name}".')


class InvalidType(RxnImportError):
    """Exception when a column contains data of the incorrect type."""

    def __init__(self, column_name: str, expected_type: str, actual_type: str):
        super().__init__(
            f'The column "{column_name}" is expected to contain '
            f"values of type {expected_type} (actual: {actual_type})."
        )


def _final_column_name_for_original(cfg: RxnImportConfig) -> str:
    """Name of the column where the original data will end up."""

    # For txt: always with "_original" postfix
    if cfg.data_format is InitialDataFormat.TXT:
        return f"{cfg.reaction_column_name}_original"

    # For csv: only add "_original" postfix if the input and output columns are identical
    if cfg.input_csv_column_name == cfg.reaction_column_name:
        return f"{cfg.reaction_column_name}_original"

    return cfg.input_csv_column_name


def _load_from_txt(cfg: RxnImportConfig) -> pd.DataFrame:
    column_original = _final_column_name_for_original(cfg)
    df: pd.DataFrame = pd.read_csv(cfg.input_file, names=[column_original])
    return df


def _load_from_csv(cfg: RxnImportConfig, separator: str) -> pd.DataFrame:
    df: pd.DataFrame = pd.read_csv(cfg.input_file, sep=separator)

    if cfg.input_csv_column_name not in df.columns:
        raise InvalidColumn(cfg.input_csv_column_name, cfg.input_file)

    # if the name of the import column and of the output column collide: rename it.
    if cfg.input_csv_column_name == cfg.reaction_column_name:
        df.rename(
            columns={cfg.input_csv_column_name: _final_column_name_for_original(cfg)},
            inplace=True,
        )

    return df


def _load_initial(cfg: RxnImportConfig) -> pd.DataFrame:
    """
    Load the initial dataframe.

    This functions (and the ones it calls) make sure to rename the column
    containing the data if it would be overwritten when transforming. The
    column with the data ends up with the name `_final_column_name_for_original`.
    """

    if cfg.data_format is InitialDataFormat.TXT:
        return _load_from_txt(cfg)
    if cfg.data_format is InitialDataFormat.CSV:
        return _load_from_csv(cfg, separator=",")
    if cfg.data_format is InitialDataFormat.TSV:
        return _load_from_csv(cfg, separator="\t")
    raise ValueError(f"Unsupported data type: {cfg.data_format}")


def reformat_smiles(
    reaction_smiles: str, fragment_bond: Optional[str]
) -> Optional[str]:
    """Import a reaction SMILES in any format and convert it to an "IBM" RXN
    SMILES with the specified fragment bond.

    `None` is returned if the reaction SMILES is not valid.
    """

    try:
        reaction = parse_any_reaction_smiles(reaction_smiles)
        return reaction.to_string(fragment_bond)
    except Exception:
        logger.info(f"Invalid reaction: {reaction_smiles}")
        return None


def _add_token(
    row: pd.Series,
    add_special_token_fn: Callable[[ReactionEquation], ReactionEquation],
    column_to_check: str,
    cfg: RxnImportConfig,
) -> str:
    """Function to use in pandas.apply to update the reaction SMILES."""

    reaction_smiles = row[cfg.reaction_column_name]
    special_flag_active = row[column_to_check]

    if not isinstance(special_flag_active, bool):
        raise InvalidType(
            column_name=column_to_check,
            expected_type="bool",
            actual_type=str(type(special_flag_active).__name__),
        )

    if not special_flag_active:
        # do nothing if the reaction is not run under light / heat / etc.
        return reaction_smiles

    reaction = parse_any_reaction_smiles(reaction_smiles)
    reaction = add_special_token_fn(reaction)

    return reaction.to_string(cfg.fragment_bond.value)


def _maybe_add_special_token(
    df: pd.DataFrame,
    cfg: RxnImportConfig,
    add_special_token_fn: Callable[[ReactionEquation], ReactionEquation],
    special_token_column: Optional[str],
) -> None:
    """
    Common function for adding light or heat tokens.

    Args:
        df: DataFrame where to update the reaction SMILES.
        cfg: config.
        add_special_token_fn: function to call on a ReactionEquation to add the required token.
        special_token_column: column to look at in order to know when to add the special tokens.
    """

    # Do nothing if no column was specified for the special token
    if special_token_column is None:
        return

    if special_token_column not in df.columns:
        raise InvalidColumn(special_token_column, cfg.input_file)

    # "partial" in order to make the function compatible with pandas.apply().
    fn = partial(
        _add_token,
        add_special_token_fn=add_special_token_fn,
        column_to_check=special_token_column,
        cfg=cfg,
    )
    df[cfg.reaction_column_name] = df.apply(fn, axis=1)


def _maybe_add_heat_token(df: pd.DataFrame, cfg: RxnImportConfig) -> None:
    """Add the heat token to the precursors of the reaction SMILES when necessary."""

    _maybe_add_special_token(
        df=df,
        cfg=cfg,
        add_special_token_fn=add_heat_token,
        special_token_column=cfg.column_for_heat,
    )


def _maybe_add_light_token(df: pd.DataFrame, cfg: RxnImportConfig) -> None:
    """Add the light token to the precursors of the reaction SMILES when necessary."""

    _maybe_add_special_token(
        df=df,
        cfg=cfg,
        add_special_token_fn=add_light_token,
        special_token_column=cfg.column_for_light,
    )


def _maybe_remove_atom_mapping(df: pd.DataFrame, cfg: RxnImportConfig) -> None:
    """Remove the atom mapping if required by the config.

    NB: This will not clean up the SMILES, i.e. "[CH3:1][CH3:2]" is converted
    to "[CH3][CH3]" and not to "CC". The standardization step will do that.
    """

    if not cfg.remove_atom_mapping:
        return

    df[cfg.reaction_column_name] = df[cfg.reaction_column_name].apply(
        remove_atom_mapping
    )


def rxn_import(cfg: RxnImportConfig) -> None:
    """
    Initial import of reaction data, as a first step of the reaction preprocessing.
    """

    df = _load_initial(cfg)

    # Reformat the reaction SMILES
    fn = partial(reformat_smiles, fragment_bond=cfg.fragment_bond.value)
    df[cfg.reaction_column_name] = df[_final_column_name_for_original(cfg)].apply(fn)

    # Discard lines containing None (for unsuccessful reaction SMILES import)
    df = df.dropna()

    # Add special tokens when necessary
    _maybe_add_light_token(df, cfg)
    _maybe_add_heat_token(df, cfg)

    # Remove atom mapping if necessary
    _maybe_remove_atom_mapping(df, cfg)

    df.to_csv(cfg.output_csv, index=False)
