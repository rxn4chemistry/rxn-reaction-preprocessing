from typing import Dict, List, Union

from rxn.chemutils.reaction_equation import ReactionEquation

from rxn.reaction_preprocessing.annotations.molecule_annotation import (
    MoleculeAnnotation,
)
from rxn.reaction_preprocessing.annotations.molecule_replacer import MoleculeReplacer

replacements: Dict[str, Union[str, List[str]]] = {
    "(C)(C)": "CC",
    "C[O-].[Na+]": "CO[Na]",
    "[C].[Pd].CCCC": ["[C].[Pd]", "CCCC"],
}


def test_replace_molecule_smiles():
    molecule_replacer = MoleculeReplacer(replacements)

    # Replaces the molecules only when they exactly match the full SMILES
    assert molecule_replacer.replace_molecule_smiles("(C)(C)") == ["CC"]
    assert molecule_replacer.replace_molecule_smiles("(C)(C)O") == ["(C)(C)O"]
    assert molecule_replacer.replace_molecule_smiles("(C)(C).O") == ["(C)(C).O"]

    # For molecules with more than one fragment, replace only when in the
    # correct order and only with '.' as a fragment separator.
    assert molecule_replacer.replace_molecule_smiles("C[O-].[Na+]") == ["CO[Na]"]
    assert molecule_replacer.replace_molecule_smiles("[Na+].C[O-]") == ["[Na+].C[O-]"]
    assert molecule_replacer.replace_molecule_smiles("C[O-]~[Na+]") == ["C[O-]~[Na+]"]

    # Example of replacing a SMILES by more than one molecule
    assert molecule_replacer.replace_molecule_smiles("[C]~[Pd]~CCCC") == [
        "[C]~[Pd]~CCCC"
    ]
    assert molecule_replacer.replace_molecule_smiles("[C].[Pd].CCCC") == [
        "[C].[Pd]",
        "CCCC",
    ]


def test_replace_in_reaction_smiles():
    molecule_replacer = MoleculeReplacer(replacements)

    # 1) No change if nothing matches
    smiles = "[Na+]~[Cl-].CC.O>>CCO"
    assert molecule_replacer.replace_in_reaction_smiles(smiles, "~") == smiles

    # 2) Simple replacement
    assert molecule_replacer.replace_in_reaction_smiles("(C)(C).O>>CCO") == "CC.O>>CCO"

    # 3) Replacement with fragment bond: no change if not given as fragments
    smiles = "CC.C[O-].[Na+].O>>CCO"
    assert molecule_replacer.replace_in_reaction_smiles(smiles) == smiles
    assert molecule_replacer.replace_in_reaction_smiles(smiles, "~") == smiles

    # 4) Replacement with fragment bond
    smiles = "CC.C[O-]~[Na+].O>>CCO"
    expected = "CC.CO[Na].O>>CCO"
    assert molecule_replacer.replace_in_reaction_smiles(smiles) == smiles
    assert molecule_replacer.replace_in_reaction_smiles(smiles, "~") == expected

    # 5) One fragment being replaced by two molecules
    smiles = "CC.[C]~[Pd]~CCCC.O>>CCO"
    expected = "CC.[C]~[Pd].CCCC.O>>CCO"
    assert molecule_replacer.replace_in_reaction_smiles(smiles) == smiles
    assert molecule_replacer.replace_in_reaction_smiles(smiles, "~") == expected

    # 5) No replacement if fragments are given in another order
    smiles = "CC.[Pd]~[C]~CCCC.O>>CCO"
    assert molecule_replacer.replace_in_reaction_smiles(smiles, "~") == smiles


def test_replace_in_reaction_equation():
    molecule_replacer = MoleculeReplacer(replacements)

    # 1) '(C)(C)' is replaced by 'CC'
    assert molecule_replacer.replace_in_reaction_equation(
        ReactionEquation(["(C)(C)", "O"], [], ["CCO"])
    ) == ReactionEquation(["CC", "O"], [], ["CCO"])

    # 2) 'C[O-].[Na+]' is replaced by 'CO[Na]', two times
    assert molecule_replacer.replace_in_reaction_equation(
        ReactionEquation(["CC", "C[O-].[Na+]", "O"], ["C[O-].[Na+]"], ["CCO"])
    ) == ReactionEquation(["CC", "CO[Na]", "O"], ["CO[Na]"], ["CCO"])

    # 3) '[C]~[Pd]~CCCC' is expanded into two molecules '[C].[Pd]' and 'CCCC'.
    assert molecule_replacer.replace_in_reaction_equation(
        ReactionEquation(["CC", "[C].[Pd].CCCC", "O"], [], ["CCO"])
    ) == ReactionEquation(["CC", "[C].[Pd]", "CCCC", "O"], [], ["CCO"])


def test_molecule_replacer_from_annotations():
    annotations = [
        MoleculeAnnotation("O", None, "accept", []),
        MoleculeAnnotation("(C)(C)", "CC", "reject", []),
        MoleculeAnnotation("[C]~CCCC~[Pd]", "[C]~[Pd].CCCC", "accept", []),
    ]

    molecule_replacer = MoleculeReplacer.from_molecule_annotations(annotations)

    # 'O' has no updated SMILES so it shouldn't be considered
    assert "O" not in molecule_replacer.replacements

    # '(C)(C)' is rejected so it shouldn't be considered
    assert "(C)(C)" not in molecule_replacer.replacements
    assert molecule_replacer.replace_molecule_smiles("(C)(C)") == ["(C)(C)"]

    # When replaced, the third annotation delivers two molecules
    assert molecule_replacer.replace_molecule_smiles("[C].CCCC.[Pd]") == [
        "[C].[Pd]",
        "CCCC",
    ]
