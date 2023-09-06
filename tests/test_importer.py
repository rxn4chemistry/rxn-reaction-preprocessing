import tempfile
from pathlib import Path
from typing import Generator

import pandas as pd
import pytest
from rxn.utilities.files import dump_list_to_file

from rxn.reaction_preprocessing.config import (
    FragmentBond,
    InitialDataFormat,
    RxnImportConfig,
)
from rxn.reaction_preprocessing.importer import InvalidColumn, rxn_import
from rxn.reaction_preprocessing.special_tokens import HEAT_TOKEN, LIGHT_TOKEN


@pytest.fixture
def input_file() -> Generator[str, None, None]:
    """Fixture to create a temporary directory and get the name of a
    (yet nonexistent) input file inside it."""
    with tempfile.TemporaryDirectory() as temporary_dir:
        temporary_path = Path(temporary_dir)
        input_file = str(temporary_path / "input.file")
        yield input_file


@pytest.fixture
def output_file() -> Generator[str, None, None]:
    """Fixture to create a temporary directory and get the name of a
    (yet nonexistent) output file inside it."""
    with tempfile.TemporaryDirectory() as temporary_dir:
        temporary_path = Path(temporary_dir)
        output_file = str(temporary_path / "output.file")
        yield output_file


def test_import_from_txt(input_file: str, output_file: str) -> None:
    # Write some reactions to a TXT file
    reactions = ["CC>>CC", "OO>>OO", "C.C.O>>CCO", "C.[Na+].[Cl-]>>C"]
    dump_list_to_file(reactions, input_file)

    # Do the initial import
    cfg = RxnImportConfig(
        input_file=input_file,
        output_csv=output_file,
        data_format=InitialDataFormat.TXT,
        reaction_column_name="rxn",
        fragment_bond=FragmentBond.TILDE,
        keep_original_rxn_column=True,
    )
    rxn_import(cfg)

    # Verify the content.
    df = pd.read_csv(output_file)
    assert df["rxn"].tolist() == reactions
    assert df["rxn_original"].tolist() == reactions


def test_import_from_csv(input_file: str, output_file: str) -> None:
    # Write some reactions to a CSV file
    reactions = ["CC>>CC", "OO>>OO", "C.C.O>>CCO", "C.[Na+].[Cl-]>>C"]
    dummy_data = ["A", "B", "C", "D"]
    pd.DataFrame({"smiles": reactions, "dummy": dummy_data}).to_csv(
        input_file, index=False
    )

    # Do the initial import
    cfg = RxnImportConfig(
        input_file=input_file,
        output_csv=output_file,
        data_format=InitialDataFormat.CSV,
        input_csv_column_name="smiles",
        reaction_column_name="rxn",
        fragment_bond=FragmentBond.TILDE,
        keep_original_rxn_column=True,
    )
    rxn_import(cfg)

    # Verify the content
    df = pd.read_csv(output_file)
    assert df["rxn"].tolist() == reactions
    assert df["smiles"].tolist() == reactions


def test_import_from_tsv(input_file: str, output_file: str) -> None:
    # Write some reactions to a TSV file
    reactions = ["CC>>CC", "OO>>OO", "C.C.O>>CCO", "C.[Na+].[Cl-]>>C"]
    dummy_data = ["A", "B", "C", "D"]
    pd.DataFrame({"smiles": reactions, "dummy": dummy_data}).to_csv(
        input_file, index=False, sep="\t"
    )

    # Do the initial import
    cfg = RxnImportConfig(
        input_file=input_file,
        output_csv=output_file,
        data_format=InitialDataFormat.TSV,
        input_csv_column_name="smiles",
        reaction_column_name="rxn",
        fragment_bond=FragmentBond.TILDE,
        keep_original_rxn_column=True,
    )
    rxn_import(cfg)

    # Verify the content
    df = pd.read_csv(output_file)
    assert df["rxn"].tolist() == reactions
    assert df["smiles"].tolist() == reactions


def test_smiles_format_is_updated(input_file: str, output_file: str) -> None:
    # The import must support both extended reaction SMILES (with "|f:" notation),
    # as well as reaction SMILES with the fragment bond "~", even in the same file.
    reactions = [
        "CC>>CC",
        "OO>>OO",
        "C.C.O>>CCO",
        "C.[Na+]~[Cl-]>>C",
        "C.[Na+].[Cl-]>>C |f:1.2|",
    ]
    expected = [
        "CC>>CC",
        "OO>>OO",
        "C.C.O>>CCO",
        "C.[Na+]~[Cl-]>>C",
        "C.[Na+]~[Cl-]>>C",
    ]

    # Write some reactions to a CSV file
    pd.DataFrame({"smiles": reactions}).to_csv(input_file, index=False)

    # Do the initial import
    cfg = RxnImportConfig(
        input_file=input_file,
        output_csv=output_file,
        data_format=InitialDataFormat.CSV,
        input_csv_column_name="smiles",
        reaction_column_name="rxn",
        fragment_bond=FragmentBond.TILDE,
        keep_original_rxn_column=True,
    )
    rxn_import(cfg)

    # Verify the content
    df = pd.read_csv(output_file)
    assert df["rxn"].tolist() == expected
    assert df["smiles"].tolist() == reactions


def test_dot_as_fragment_bond(input_file: str, output_file: str) -> None:
    # When the dot is given as a fragment bond, other fragment information is allowed in
    # the input, but will not be present anymore in the output.
    reactions = [
        "CC>>CC",
        "OO>>OO",
        "C.C.O>>CCO",
        "C.[Na+]~[Cl-]>>C",
        "C.[Na+].[Cl-]>>C |f:1.2|",
    ]
    expected = [
        "CC>>CC",
        "OO>>OO",
        "C.C.O>>CCO",
        "C.[Na+].[Cl-]>>C",
        "C.[Na+].[Cl-]>>C",
    ]

    # Write the reactions to a CSV file
    pd.DataFrame({"smiles": reactions}).to_csv(input_file, index=False)

    # Do the initial import
    cfg = RxnImportConfig(
        input_file=input_file,
        output_csv=output_file,
        data_format=InitialDataFormat.CSV,
        input_csv_column_name="smiles",
        reaction_column_name="rxn",
        fragment_bond=FragmentBond.DOT,
        keep_original_rxn_column=True,
    )
    rxn_import(cfg)

    # Verify the content
    df = pd.read_csv(output_file)
    assert df["rxn"].tolist() == expected
    assert df["smiles"].tolist() == reactions


def test_rename_column_to_avoid_overwrite(input_file: str, output_file: str) -> None:
    reactions = ["CC>>CC", "OO>>OO", "C.C.O>>CCO", "C.[Na+].[Cl-]>>C |f:1.2|"]
    expected = ["CC>>CC", "OO>>OO", "C.C.O>>CCO", "C.[Na+]~[Cl-]>>C"]

    # Same test as above, except that the original column name is the same
    # as the one desired for output. In this case, the original column is renamed.
    pd.DataFrame({"rxn": reactions}).to_csv(input_file, index=False)

    # Do the initial import
    cfg = RxnImportConfig(
        input_file=input_file,
        output_csv=output_file,
        data_format=InitialDataFormat.CSV,
        input_csv_column_name="rxn",
        reaction_column_name="rxn",
        fragment_bond=FragmentBond.TILDE,
        keep_original_rxn_column=True,
    )
    rxn_import(cfg)

    # Verify the content
    df = pd.read_csv(output_file)
    assert df["rxn"].tolist() == expected
    assert df["rxn_original"].tolist() == reactions


def test_raises_when_column_name_is_incorrect(
    input_file: str, output_file: str
) -> None:
    # Write some reactions to a CSV file
    reactions = ["CC>>CC", "OO>>OO", "C.C.O>>CCO", "C.[Na+].[Cl-]>>C"]
    pd.DataFrame({"smiles": reactions}).to_csv(input_file, index=False)

    # Note that the input column name does not exist in the created CSV
    cfg = RxnImportConfig(
        input_file=input_file,
        output_csv=output_file,
        data_format=InitialDataFormat.CSV,
        input_csv_column_name="incorrect_name",
        reaction_column_name="rxn",
        fragment_bond=FragmentBond.TILDE,
    )

    with pytest.raises(InvalidColumn):
        rxn_import(cfg)


def test_light_and_heat(input_file: str, output_file: str) -> None:
    reactions = ["CC>>CC", "OO>>OO", "C.C.O>>CCO", "C=C>>CC"]
    has_light = [False, True, "yes", "no"]
    has_heat = [False, True, "no", "yes"]
    # Add light token to second and third reactions, heat token to second and fourth
    expected = [
        "CC>>CC",
        f"OO.{LIGHT_TOKEN}.{HEAT_TOKEN}>>OO",
        f"C.C.O.{LIGHT_TOKEN}>>CCO",
        f"C=C.{HEAT_TOKEN}>>CC",
    ]

    # Write the reactions to a CSV file
    pd.DataFrame({"smiles": reactions, "heat": has_heat, "light": has_light}).to_csv(
        input_file, index=False
    )

    # Note that the input column name does not exist in the created CSV
    cfg = RxnImportConfig(
        input_file=input_file,
        output_csv=output_file,
        data_format=InitialDataFormat.CSV,
        input_csv_column_name="smiles",
        reaction_column_name="rxn",
        fragment_bond=FragmentBond.TILDE,
        column_for_light="light",
        column_for_heat="heat",
        keep_original_rxn_column=True,
    )

    rxn_import(cfg)

    # Verify the content
    df = pd.read_csv(output_file)
    assert df["rxn"].tolist() == expected
    assert df["smiles"].tolist() == reactions


def test_light_and_heat_with_inexisting_column(
    input_file: str, output_file: str
) -> None:
    reactions = ["CC>>CC", "OO>>OO", "C.C.O>>CCO", "C=C>>CC"]
    pd.DataFrame({"smiles": reactions}).to_csv(input_file, index=False)

    # Note that the input column name does not exist in the created CSV
    cfg = RxnImportConfig(
        input_file=input_file,
        output_csv=output_file,
        data_format=InitialDataFormat.CSV,
        input_csv_column_name="smiles",
        reaction_column_name="rxn",
        fragment_bond=FragmentBond.TILDE,
        column_for_light="non_exising_column",
        keep_original_rxn_column=True,
    )

    with pytest.raises(InvalidColumn):
        rxn_import(cfg)


def test_atom_mapping_removal(input_file: str, output_file: str) -> None:
    reactions = ["[CH3:1][CH3:2]>>[CH2:1]=[CH2:2]", "OO>>OO", "[CH4:1].C.O>>[CH3:1]CO"]
    expected = ["[CH3][CH3]>>[CH2]=[CH2]", "OO>>OO", "[CH4].C.O>>[CH3]CO"]
    pd.DataFrame({"smiles": reactions}).to_csv(input_file, index=False)

    # 1) atom mapping removal active
    cfg = RxnImportConfig(
        input_file=input_file,
        output_csv=output_file,
        data_format=InitialDataFormat.CSV,
        input_csv_column_name="smiles",
        reaction_column_name="rxn",
        remove_atom_mapping=True,
        fragment_bond=FragmentBond.TILDE,
        keep_original_rxn_column=True,
    )
    rxn_import(cfg)
    assert pd.read_csv(output_file)["rxn"].tolist() == expected

    # 2) no atom mapping removal -> gives back the original reactions
    cfg.remove_atom_mapping = False
    rxn_import(cfg)
    assert pd.read_csv(output_file)["rxn"].tolist() == reactions


def test_invalid_reactions_are_ignored(input_file: str, output_file: str) -> None:
    # Write some reactions to a CSV file. Two of them (indices 1 and 4) are invalid.
    reactions = [
        "CC>>CC",
        "invalid",
        "OO>>OO",
        "C.C.O>>CCO",
        "CC>>CO>>CN",
        "C.[Na+].[Cl-]>>C",
    ]
    valid_reactions = [reactions[i] for i in [0, 2, 3, 5]]
    pd.DataFrame({"smiles": reactions}).to_csv(input_file, index=False)

    # Do the initial import
    cfg = RxnImportConfig(
        input_file=input_file,
        output_csv=output_file,
        data_format=InitialDataFormat.CSV,
        input_csv_column_name="smiles",
        reaction_column_name="rxn",
        fragment_bond=FragmentBond.TILDE,
        keep_original_rxn_column=True,
    )
    rxn_import(cfg)

    # Verify the content
    df = pd.read_csv(output_file)
    assert df["rxn"].tolist() == valid_reactions
    assert df["smiles"].tolist() == valid_reactions


def test_remove_original_column_txt(input_file: str, output_file: str) -> None:
    # Write some reactions to a TXT file
    reactions = ["CC>>CC", "OO>>OO", "C.C.O>>CCO", "C.[Na+].[Cl-]>>C"]
    dump_list_to_file(reactions, input_file)

    # Do the initial import
    cfg = RxnImportConfig(
        input_file=input_file,
        output_csv=output_file,
        data_format=InitialDataFormat.TXT,
        reaction_column_name="rxn",
        fragment_bond=FragmentBond.TILDE,
        keep_original_rxn_column=False,
    )
    rxn_import(cfg)

    # Verify that there is only one column
    df = pd.read_csv(output_file)
    assert list(df.columns) == ["rxn"]


def test_remove_original_column_csv(input_file: str, output_file: str) -> None:
    # Write some reactions to a CSV file
    reactions = ["CC>>CC", "OO>>OO", "C.C.O>>CCO", "C.[Na+].[Cl-]>>C"]
    dummy_data = ["A", "B", "C", "D"]
    pd.DataFrame({"smiles": reactions, "dummy": dummy_data}).to_csv(
        input_file, index=False
    )

    # Do the initial import
    cfg = RxnImportConfig(
        input_file=input_file,
        output_csv=output_file,
        data_format=InitialDataFormat.CSV,
        input_csv_column_name="smiles",
        reaction_column_name="rxn",
        fragment_bond=FragmentBond.TILDE,
        keep_original_rxn_column=False,
    )
    rxn_import(cfg)

    # Verify that there is only the two required columns - no "smiles" anymore
    df = pd.read_csv(output_file)
    assert set(df.columns) == {"rxn", "dummy"}
