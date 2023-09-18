# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2022
# ALL RIGHTS RESERVED
from pathlib import Path
from typing import Generator

import pytest
from rxn.utilities.files import (
    dump_list_to_file,
    load_list_from_file,
    named_temporary_directory,
)

from rxn.reaction_preprocessing import MixedReactionFilter, Preprocessor


@pytest.fixture
def tmp_dir() -> Generator[Path, None, None]:
    with named_temporary_directory() as dir_path:
        yield dir_path


def test_preprocessor(tmp_dir: Path) -> None:
    # Create fake input file
    input_path = tmp_dir / "input.csv"
    dump_list_to_file(
        [
            "rxn,class",
            "C.O.O>C>CO,11",
            "C~N.N>>CNN,0",
            "[14C]Cl.O[Na]>O>[Na]Cl.[14C]O,2",
            "[14C]Cl.O[Na]>>[Na]Cl,2",
            "[C].C.[O--].[O--].O.O=[N+]([O-])c1cc(-c2nc3ccccc3o2)ccc1F.C.C>O>O.C,6",
            # In an earlier version, preprocessing failed for rare cases when RDKit cannot
            # convert concatenated SMILES into an RDKit Mol, as for the SMILES below.
            # Was related to https://github.com/rdkit/rdkit/issues/4266
            "c1ccccc1.C123C45C16C21C34C561>>CC,1",
            "=)(/?,0",
        ],
        input_path,
    )

    mrf = MixedReactionFilter(
        max_reactants=5,
        max_agents=0,
        max_products=1,
        min_reactants=2,
        min_agents=0,
        min_products=1,
        max_reactants_tokens=300,
        max_agents_tokens=0,
        max_products_tokens=200,
        max_absolute_formal_charge=2,
    )

    preprocessor = Preprocessor(
        mixed_reaction_filter=mrf, reaction_column_name="rxn", fragment_bond="~"
    )

    output_path = tmp_dir / "output.csv"
    preprocessor.process_file(input_path, output_path)

    assert load_list_from_file(output_path) == [
        "rxn,class",
        "C.O>>CO,11",  # removed the duplicated O and C
        "C~N.N>>CNN,0",  # kept it as such
        "O[Na].[14C]Cl>>[Na]Cl,2",  # sorted the compounds
        "C123C45C16C21C34C561.c1ccccc1>>CC,1",
    ]
