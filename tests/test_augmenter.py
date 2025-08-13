# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2022
# ALL RIGHTS RESERVED
import pandas as pd
import pytest

from rxn import reaction_preprocessing as rrp
from rxn.reaction_preprocessing import Augmenter
from rxn.reaction_preprocessing.utils import ReactionSection, reset_random_seed


@pytest.fixture
def augmenter() -> Augmenter:
    reset_random_seed()
    df = pd.DataFrame(
        {
            "rxn": [
                "[Na]Cl~[K+].[K+].CCCC.CN1C=NC2=C1C(=O)N(C(=O)N2C)C>>CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                "CC(C)(C)O[K]~CCC.O>>CC",
                "CC>>CC",  # this is to check what happens when doing augmentation does not make much sense
            ],
        }
    )
    return Augmenter(df, "rxn", fragment_bond="~")


def test_restricted(augmenter: Augmenter) -> None:
    new_df = augmenter.augment(rrp.RandomType.restricted, permutations=1)

    expected = [
        "[Cl][Na]~[K+].[K+].CCCC.CN1C2=C(C(=O)N(C)C1=O)N(C)C=N2>>CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "CC(C)(C)[O][K]~CCC.O>>CC",
        "CC>>CC",
    ]
    assert new_df["rxn_restricted"].tolist() == expected


def test_unrestricted(augmenter: Augmenter) -> None:
    new_df = augmenter.augment(rrp.RandomType.unrestricted, permutations=1)

    expected = [
        "[Na][Cl]~[K+].[K+].C(C)CC.O=C1N(C(=O)N(C)C2=C1N(C)C=N2)C>>CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "CC(C)([O][K])C~C(C)C.O>>CC",
        "CC>>CC",
    ]
    assert new_df["rxn_unrestricted"].tolist() == expected


def test_rotated(augmenter: Augmenter) -> None:
    new_df = augmenter.augment(rrp.RandomType.rotated, permutations=1)

    expected = [
        "[Cl][Na]~[K+].[K+].CCCC.C1(=O)N(C)C2=C(N(C)C=N2)C(=O)N1C>>CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "[K][O]C(C)(C)C~CCC.O>>CC",
        "CC>>CC",
    ]
    assert new_df["rxn_rotated"].tolist() == expected


def test_molecules_order(augmenter: Augmenter) -> None:
    new_df = augmenter.augment(rrp.RandomType.molecules, permutations=1)

    expected = [
        "[Na]Cl~[K+].CN1C=NC2=C1C(=O)N(C(=O)N2C)C.[K+].CCCC>>CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "CC(C)(C)O[K]~CCC.O>>CC",
        "CC>>CC",
    ]
    assert new_df["rxn_molecules"].tolist() == expected


def test_multiple_augmentation_molecules(augmenter: Augmenter) -> None:
    new_df = augmenter.augment(rrp.RandomType.molecules, permutations=3)

    expected = [
        "[Na]Cl~[K+].CN1C=NC2=C1C(=O)N(C(=O)N2C)C.[K+].CCCC>>CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C.[K+].CCCC.[Na]Cl~[K+]>>CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "[Na]Cl~[K+].[K+].CCCC.CN1C=NC2=C1C(=O)N(C(=O)N2C)C>>CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "O.CC(C)(C)O[K]~CCC>>CC",
        "CC(C)(C)O[K]~CCC.O>>CC",
        "CC>>CC",
    ]
    assert new_df["rxn_molecules"].tolist() == expected


def test_multiple_augmentation_rotated(augmenter: Augmenter) -> None:
    new_df = augmenter.augment(rrp.RandomType.rotated, permutations=3)

    expected = [
        "[Cl][Na]~[K+].[K+].CCCC.C1(=O)N(C)C2=C(N(C)C=N2)C(=O)N1C>>CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "[Cl][Na]~[K+].[K+].C(CC)C.N1=CN(C)C2=C1N(C)C(=O)N(C)C2=O>>CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "[Na][Cl]~[K+].[K+].CCCC.CN1C2=C(N(C)C=N2)C(=O)N(C)C1=O>>CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "[K][O]C(C)(C)C~CCC.O>>CC",
        "[O]([K])C(C)(C)C~CCC.O>>CC",
        "CC(C)(C)[O][K]~CCC.O>>CC",
        "CC>>CC",
        "CC>>CC",
        "CC>>CC",
    ]
    assert new_df["rxn_rotated"].tolist() == expected


def test_multiple_augmentation_rotated_on_products(augmenter: Augmenter) -> None:
    new_df = augmenter.augment(
        rrp.RandomType.rotated,
        rxn_section_to_augment=ReactionSection.products,
        permutations=3,
    )

    expected = [
        "[Na]Cl~[K+].[K+].CCCC.CN1C=NC2=C1C(=O)N(C(=O)N2C)C>>C1(=O)N(C)C(=O)C2=C(N=CN2C)N1C",
        "[Na]Cl~[K+].[K+].CCCC.CN1C=NC2=C1C(=O)N(C(=O)N2C)C>>CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
        "[Na]Cl~[K+].[K+].CCCC.CN1C=NC2=C1C(=O)N(C(=O)N2C)C>>C1=NC2=C(N1C)C(=O)N(C)C(=O)N2C",
        "CC(C)(C)O[K]~CCC.O>>CC",
        "CC(C)(C)O[K]~CCC.O>>CC",
        "CC(C)(C)O[K]~CCC.O>>CC",
        "CC>>CC",
        "CC>>CC",
        "CC>>CC",
    ]
    assert new_df["rxn_rotated"].tolist() == expected
