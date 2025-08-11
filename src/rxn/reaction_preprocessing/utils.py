# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED
import logging
import random
from enum import Enum, auto
from pathlib import Path
from typing import Iterable, List, Set

import attr
import numpy
from rdkit.Chem import GetFormalCharge, Mol
from rxn.chemutils.conversion import smiles_to_mol
from rxn.chemutils.reaction_equation import ReactionEquation
from rxn.utilities.files import PathLike
from rxn.utilities.logging import LoggingFormat


class DataSplit(Enum):
    TRAIN = auto()
    VALIDATION = auto()
    TEST = auto()


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
    return Path(__file__).parent.resolve() / "data"


def standardization_files_directory() -> Path:
    """
    Returns the path to the data directory at the root of the repository
    """
    return data_directory() / "standardization-files"


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
    def from_reaction_equation(cls, reaction: ReactionEquation) -> "MolEquation":
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
    return {atom.GetSymbol() for mol in mols for atom in mol.GetAtoms()}  # type:ignore


def add_custom_logger_to_file(log_file: PathLike) -> None:
    """
    Set up logging to a file.

    This is a bit more complex than usual because hydra sets up the logger
    automattically, and it is not possible to disable it.

    Args:
        log_file: file where to save the logs.
    """
    root_logger = logging.getLogger()
    fh = logging.FileHandler(log_file, mode="w")
    fh.setLevel(logging.INFO)
    root_logger.addHandler(fh)


def overwrite_logging_format() -> None:
    """
    Reset the log format to our default, instead of using the hydra default.
    """
    root_logger = logging.getLogger()
    for handler in root_logger.handlers:
        formatter = logging.Formatter(LoggingFormat.BASIC.value)
        handler.setFormatter(formatter)
