# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED
import random
from enum import auto
from enum import Enum
from pathlib import Path
from typing import Iterable
from typing import List
from typing import Set

import attr
import numpy
from rdkit.Chem import GetFormalCharge
from rdkit.Chem import Mol
from rxn_chemutils.conversion import smiles_to_mol
from rxn_chemutils.reaction_equation import ReactionEquation


class RandomType(Enum):
    molecules = auto()
    unrestricted = auto()
    restricted = auto()
    rotated = auto()


class ReactionSection(Enum):
    precursors = auto()
    products = auto()


def root_directory() -> Path:
    """
    Returns the path to the root directory of the repository
    """
    return Path(__file__).parent.parent.resolve()


def data_directory() -> Path:
    """
    Returns the path to the data directory at the root of the repository
    """
    return Path(__file__).parent.resolve() / 'data'


def standardization_files_directory() -> Path:
    """
    Returns the path to the data directory at the root of the repository
    """
    return data_directory() / 'standardization-files'


def reset_random_seed() -> None:
    random.seed(42)
    numpy.random.seed(42)


@attr.s(auto_attribs=True)
class MolEquation:
    """
    Same as a ReactionEquation, except that RDKit Mol objects are stored
    instead of the SMILES.
    """
    reactants: List[Mol]
    agents: List[Mol]
    products: List[Mol]

    @classmethod
    def from_reaction_equation(cls, reaction: ReactionEquation) -> 'MolEquation':
        return cls(
            reactants=[smiles_to_mol(s) for s in reaction.reactants],
            agents=[smiles_to_mol(s) for s in reaction.agents],
            products=[smiles_to_mol(s) for s in reaction.products],
        )


def get_formal_charge_for_mols(mols: Iterable[Mol]) -> int:
    """Get the formal charge for a group of RDKit Mols."""
    return sum(GetFormalCharge(mol) for mol in mols)


def get_atoms_for_mols(mols: Iterable[Mol]) -> Set[str]:
    """Get the set of atoms for a list of RDKit Mols."""
    return {atom.GetSymbol() for mol in mols for atom in mol.GetAtoms()}
