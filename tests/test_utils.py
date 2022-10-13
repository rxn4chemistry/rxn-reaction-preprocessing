import pytest
from rxn.chemutils.exceptions import InvalidSmiles
from rxn.chemutils.reaction_equation import ReactionEquation

from rxn.reaction_preprocessing.utils import MolEquation


def test_mol_equation() -> None:
    # Dummy example to check if the MolEquation contains the correct number of atoms
    reaction_equation = ReactionEquation.from_string("CC.O>Cl>CCO")
    mol_equation = MolEquation.from_reaction_equation(reaction_equation)
    assert len(mol_equation.reactants) == 2
    assert len(mol_equation.agents) == 1
    assert len(mol_equation.products) == 1
    assert mol_equation.reactants[0].GetNumAtoms() == 2
    assert mol_equation.reactants[1].GetNumAtoms() == 1
    assert mol_equation.agents[0].GetNumAtoms() == 1
    assert mol_equation.products[0].GetNumAtoms() == 3

    # Raises InvalidSmiles if any compound is invalid
    invalid_smiles_equation = ReactionEquation.from_string("JCC.C>>CCO")
    with pytest.raises(InvalidSmiles):
        _ = MolEquation.from_reaction_equation(invalid_smiles_equation)
