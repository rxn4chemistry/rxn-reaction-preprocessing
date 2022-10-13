# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED
from pathlib import Path

import pandas as pd
import pytest

from rxn.reaction_preprocessing import Standardizer
from rxn.reaction_preprocessing.annotations.molecule_annotation import load_annotations

annotations_file = str(
    Path(__file__).parent / "annotations/test_molecule_annotations.json"
)


@pytest.fixture
def standardizer() -> Standardizer:
    df = pd.DataFrame(
        {
            "rxn": [
                "O=C1CCC2(CC1)OCCO2>>O=C1CCC2(C[C@H]1C)OCCO2",  # stereo only in product
                "O=C1CCC2(CC1)OCCO2>>O=C1CCC2(C[C@@H]1C)OCCO2",  # stereo only in product
                "[Na]Cl.CC[Zn]CC~Cc1ccccc1>>[Na]Cl",  # substitution needed
                "[Na]Cl.Cc1ccccc1~CC[Zn]CC>>[Na]Cl",  # substitution needed but performed only if canonicalization
                "CC.CCC>>CCO",  # no substitution
                "CC.[Na+5].CC>>[Na+]~[OH-]",  # invalid smiles -> '>>'
                "CC(C)(C)O[K].CCO~CCO>>[Li]O",  # rejected reaction -> '>>'
                r"CC(=O)/C=C(\C)O[V](=O)O/C(C)=C/C(C)=O.CCCC[N+](CCCC)(CCCC)CCCC~CCCC[N+](CCCC)(CCCC)CCCC~C1(=C(SC(=S)S1)[S-])[S-]~C1(=C(SC(=S)S1)[S-])[S-]~[Pd+2]>>O[K]",
                # requires annotation -> '>>'
            ],
        }
    )
    annotations = load_annotations(annotations_file)
    return Standardizer(df, annotations, True, "rxn", fragment_bond="~")


@pytest.fixture
def standardizer_without_fragment() -> Standardizer:
    df = pd.DataFrame(
        {
            "rxn": [
                "O=C1CCC2(CC1)OCCO2>>O=C1CCC2(C[C@H]1C)OCCO2",  # stereo only in product
                "[Na]Cl.CC[Zn]CC.Cc1ccccc1>>[Na]Cl",  # requires annotation -> '>>'
                "[Na]Cl.Cc1ccccc1.CC[Zn]CC>>[Na]Cl",  # requires annotation -> '>>'
                "CC.CCC>>CCO",  # no substitution
                "CC[Na5+].CC>>[Na+].[OH-]",  # invalid smiles -> '>>'
                "CC(C)(C)O[K].CCO.CCO>>[Li]O",  # not rejected
                r"CC(=O)/C=C(\C)O[V](=O)O/C(C)=C/C(C)=O.CCCC[N+](CCCC)(CCCC)CCCC.CCCC[N+](CCCC)(CCCC)CCCC.C1(=C(SC(=S)S1)[S-])[S-].C1(=C(SC(=S)S1)[S-])[S-].[Pd+2]>>O[K]",
                # requires annotation -> '>>'
            ],
        }
    )
    annotations = load_annotations(annotations_file)
    return Standardizer(df, annotations, True, "rxn", fragment_bond=None)


def test_standardization(standardizer: Standardizer) -> None:
    new_df = standardizer.standardize().df
    converted_rxns = [
        "O=C1CCC2(CC1)OCCO2>>C[C@@H]1CC2(CCC1=O)OCCO2",
        "O=C1CCC2(CC1)OCCO2>>C[C@H]1CC2(CCC1=O)OCCO2",
        "[Na]Cl.CC[Zn]CC.Cc1ccccc1>>[Na]Cl",
        "[Na]Cl.CC[Zn]CC.Cc1ccccc1>>[Na]Cl",
        "CC.CCC>>CCO",
        ">>",  # invalid smiles
        ">>",  # rejected reaction
        ">>",  # annotation needed
    ]

    assert all(
        [
            new_df["rxn"].values[i] == converted_rxns[i]
            for i in range(len(converted_rxns))
        ]
    )


def test_standardization_non_canonical(standardizer: Standardizer) -> None:
    new_df = standardizer.standardize(canonicalize=False).df
    converted_rxns = [
        "O=C1CCC2(CC1)OCCO2>>O=C1CCC2(C[C@H]1C)OCCO2",
        "O=C1CCC2(CC1)OCCO2>>O=C1CCC2(C[C@@H]1C)OCCO2",
        "[Na]Cl.CC[Zn]CC.Cc1ccccc1>>[Na]Cl",
        ">>",  # does not reorder fragments, so it will need annotation
        "CC.CCC>>CCO",
        ">>",  # invalid smiles
        ">>",  # rejected reaction
        ">>",  # annotation needed
    ]

    assert all(
        [
            new_df["rxn"].values[i] == converted_rxns[i]
            for i in range(len(converted_rxns))
        ]
    )


def test_standardization_without_fragment(
    standardizer_without_fragment: Standardizer,
) -> None:
    new_df = standardizer_without_fragment.standardize().df

    converted_rxns = [
        "O=C1CCC2(CC1)OCCO2>>C[C@@H]1CC2(CCC1=O)OCCO2",
        ">>",  # does not find 'CC[Zn]CC' alone. Needs annotation
        ">>",
        "CC.CCC>>CCO",
        ">>",  # invalid smiles
        "CC(C)(C)O[K].CCO.CCO>>[Li]O",
        # not rejected: the check is done at the molecule level and not at the reaction level
        ">>",
    ]
    assert all(
        [
            new_df["rxn"].values[i] == converted_rxns[i]
            for i in range(len(converted_rxns))
        ]
    )


def test_standardization_without_discarding_unannotated(
    standardizer_without_fragment: Standardizer,
) -> None:
    standardizer_without_fragment.molecule_standardizer.discard_unannotated_metals = (
        False
    )

    # To compare with the previous test: here the
    new_df = standardizer_without_fragment.standardize().df

    converted_rxns = [
        "O=C1CCC2(CC1)OCCO2>>C[C@@H]1CC2(CCC1=O)OCCO2",
        "[Na]Cl.CC[Zn]CC.Cc1ccccc1>>[Na]Cl",
        "[Na]Cl.Cc1ccccc1.CC[Zn]CC>>[Na]Cl",
        "CC.CCC>>CCO",
        ">>",  # invalid smiles
        "CC(C)(C)O[K].CCO.CCO>>[Li]O",
        r"CC(=O)/C=C(\C)O[V](=O)O/C(C)=C/C(C)=O.CCCC[N+](CCCC)(CCCC)CCCC.CCCC[N+](CCCC)(CCCC)CCCC.S=c1sc([S-])c([S-])s1.S=c1sc([S-])c([S-])s1.[Pd+2]>>O[K]",
    ]
    assert all(
        [
            new_df["rxn"].values[i] == converted_rxns[i]
            for i in range(len(converted_rxns))
        ]
    )


def test_standardization_remove_stereo_when_only_in_product(
    standardizer: Standardizer,
) -> None:
    standardizer.remove_stereo_if_not_defined_in_precursors = True
    new_df = standardizer.standardize().df
    converted_rxns = [
        "O=C1CCC2(CC1)OCCO2>>CC1CC2(CCC1=O)OCCO2",
        "O=C1CCC2(CC1)OCCO2>>CC1CC2(CCC1=O)OCCO2",
        "[Na]Cl.CC[Zn]CC.Cc1ccccc1>>[Na]Cl",
        "[Na]Cl.CC[Zn]CC.Cc1ccccc1>>[Na]Cl",
        "CC.CCC>>CCO",
        ">>",  # invalid smiles
        ">>",  # rejected reaction
        ">>",  # annotation needed
    ]

    assert all(
        [
            new_df["rxn"].values[i] == converted_rxns[i]
            for i in range(len(converted_rxns))
        ]
    )
