# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED
from pathlib import Path

from rxn.reaction_preprocessing import Standardizer
from rxn.reaction_preprocessing.annotations.molecule_annotation import load_annotations

annotations_file = str(
    Path(__file__).parent / "annotations/test_molecule_annotations.json"
)


def test_standardization() -> None:
    annotations = load_annotations(annotations_file)
    standardizer = Standardizer(annotations, True, "rxn", fragment_bond="~")
    input_reactions = [
        "O=C1CCC2(CC1)OCCO2>>O=C1CCC2(C[C@H]1C)OCCO2",  # stereo only in product
        "O=C1CCC2(CC1)OCCO2>>O=C1CCC2(C[C@@H]1C)OCCO2",  # stereo only in product
        "[Na]Cl.CC[Zn]CC~Cc1ccccc1>>[Na]Cl",  # substitution needed
        "[Na]Cl.Cc1ccccc1~CC[Zn]CC>>[Na]Cl",  # substitution needed but performed only if canonicalization
        "CC.CCC>>CCO",
        "CC.[NaK].CC>>[Na+]~[OH-]",
        "CC(C)(C)O[K].CCO~CCO>>[Li]O",
        r"CC(=O)/C=C(\C)O[V](=O)O/C(C)=C/C(C)=O.CCCC[N+](CCCC)(CCCC)CCCC~CCCC[N+](CCCC)(CCCC)CCCC~C1(=C(SC(=S)S1)[S-])[S-]~C1(=C(SC(=S)S1)[S-])[S-]~[Pd+2]>>O[K]",
    ]
    expected_rxns = [
        "O=C1CCC2(CC1)OCCO2>>C[C@@H]1CC2(CCC1=O)OCCO2",
        "O=C1CCC2(CC1)OCCO2>>C[C@H]1CC2(CCC1=O)OCCO2",
        "[Na]Cl.CC[Zn]CC.Cc1ccccc1>>[Na]Cl",
        "[Na]Cl.CC[Zn]CC.Cc1ccccc1>>[Na]Cl",
        "CC.CCC>>CCO",
        ">>",  # invalid smiles
        ">>",  # rejected reaction
        ">>",  # annotation needed
    ]

    output_reactions = [standardizer.standardize_small(s) for s in input_reactions]
    assert output_reactions == expected_rxns


def test_standardization_without_fragment() -> None:
    annotations = load_annotations(annotations_file)
    standardizer = Standardizer(annotations, True, "rxn", fragment_bond=None)
    input_reactions = [
        "O=C1CCC2(CC1)OCCO2>>O=C1CCC2(C[C@H]1C)OCCO2",  # stereo only in product
        "[Na]Cl.CC[Zn]CC.Cc1ccccc1>>[Na]Cl",  # requires annotation -> '>>'
        "[Na]Cl.Cc1ccccc1.CC[Zn]CC>>[Na]Cl",  # requires annotation -> '>>'
        "CC.CCC>>CCO",  # no substitution
        "CC[Na5+].CC>>[Na+].[OH-]",  # invalid smiles -> '>>'
        "CC(C)(C)O[K].CCO.CCO>>[Li]O",  # not rejected
        r"CC(=O)/C=C(\C)O[V](=O)O/C(C)=C/C(C)=O.CCCC[N+](CCCC)(CCCC)CCCC.CCCC[N+](CCCC)(CCCC)CCCC.C1(=C(SC(=S)S1)[S-])[S-].C1(=C(SC(=S)S1)[S-])[S-].[Pd+2]>>O[K]",
    ]

    expected_rxns = [
        "O=C1CCC2(CC1)OCCO2>>C[C@@H]1CC2(CCC1=O)OCCO2",
        ">>",  # does not find 'CC[Zn]CC' alone. Needs annotation
        ">>",
        "CC.CCC>>CCO",
        ">>",  # invalid smiles
        "CC(C)(C)O[K].CCO.CCO>>[Li]O",  # not rejected: the check is done at the molecule level and not at the reaction level
        ">>",  # requires annotation -> '>>'
    ]
    output_reactions = [standardizer.standardize_small(s) for s in input_reactions]
    assert output_reactions == expected_rxns


def test_standardization_without_discarding_unannotated() -> None:
    standardizer = Standardizer(
        annotations=[],
        discard_unannotated_metals=False,
        reaction_column_name="rxn",
        fragment_bond=None,
    )

    input_reactions = [
        "O=C1CCC2(CC1)OCCO2>>O=C1CCC2(C[C@H]1C)OCCO2",  # stereo only in product
        "[Na]Cl.CC[Zn]CC.Cc1ccccc1>>[Na]Cl",  # requires annotation -> '>>'
        "[Na]Cl.Cc1ccccc1.CC[Zn]CC>>[Na]Cl",  # requires annotation -> '>>'
        "CC.CCC>>CCO",  # no substitution
        "CC[Na5+].CC>>[Na+].[OH-]",  # invalid smiles -> '>>'
        "CC(C)(C)O[K].CCO.CCO>>[Li]O",  # not rejected
        r"CC(=O)/C=C(\C)O[V](=O)O/C(C)=C/C(C)=O.CCCC[N+](CCCC)(CCCC)CCCC.CCCC[N+](CCCC)(CCCC)CCCC.C1(=C(SC(=S)S1)[S-])[S-].C1(=C(SC(=S)S1)[S-])[S-].[Pd+2]>>O[K]",
    ]

    expected_rxns = [
        "O=C1CCC2(CC1)OCCO2>>C[C@@H]1CC2(CCC1=O)OCCO2",
        "[Na]Cl.CC[Zn]CC.Cc1ccccc1>>[Na]Cl",
        "[Na]Cl.Cc1ccccc1.CC[Zn]CC>>[Na]Cl",
        "CC.CCC>>CCO",
        ">>",  # invalid smiles
        "CC(C)(C)O[K].CCO.CCO>>[Li]O",
        r"CC(=O)/C=C(\C)O[V](=O)O/C(C)=C/C(C)=O.CCCC[N+](CCCC)(CCCC)CCCC.CCCC[N+](CCCC)(CCCC)CCCC.S=c1sc([S-])c([S-])s1.S=c1sc([S-])c([S-])s1.[Pd+2]>>O[K]",
    ]
    output_reactions = [standardizer.standardize_small(s) for s in input_reactions]
    assert output_reactions == expected_rxns


def test_standardization_remove_stereo_when_only_in_product() -> None:
    annotations = load_annotations(annotations_file)
    standardizer = Standardizer(
        annotations=annotations,
        discard_unannotated_metals=True,
        reaction_column_name="rxn",
        fragment_bond="~",
        remove_stereo_if_not_defined_in_precursors=True,
    )
    input_reactions = [
        "O=C1CCC2(CC1)OCCO2>>O=C1CCC2(C[C@H]1C)OCCO2",  # stereo only in product
        "O=C1CCC2(CC1)OCCO2>>O=C1CCC2(C[C@@H]1C)OCCO2",  # stereo only in product
        "[Na]Cl.CC[Zn]CC~Cc1ccccc1>>[Na]Cl",  # substitution needed
        "[Na]Cl.Cc1ccccc1~CC[Zn]CC>>[Na]Cl",  # substitution needed but performed only if canonicalization
        "CC.CCC>>CCO",
        "CC.[NaK].CC>>[Na+]~[OH-]",
        "CC(C)(C)O[K].CCO~CCO>>[Li]O",
        r"CC(=O)/C=C(\C)O[V](=O)O/C(C)=C/C(C)=O.CCCC[N+](CCCC)(CCCC)CCCC~CCCC[N+](CCCC)(CCCC)CCCC~C1(=C(SC(=S)S1)[S-])[S-]~C1(=C(SC(=S)S1)[S-])[S-]~[Pd+2]>>O[K]",
    ]
    expected_rxns = [
        "O=C1CCC2(CC1)OCCO2>>CC1CC2(CCC1=O)OCCO2",
        "O=C1CCC2(CC1)OCCO2>>CC1CC2(CCC1=O)OCCO2",
        "[Na]Cl.CC[Zn]CC.Cc1ccccc1>>[Na]Cl",
        "[Na]Cl.CC[Zn]CC.Cc1ccccc1>>[Na]Cl",
        "CC.CCC>>CCO",
        ">>",  # invalid smiles
        ">>",  # rejected reaction
        ">>",  # annotation needed
    ]
    output_reactions = [standardizer.standardize_small(s) for s in input_reactions]
    assert output_reactions == expected_rxns
