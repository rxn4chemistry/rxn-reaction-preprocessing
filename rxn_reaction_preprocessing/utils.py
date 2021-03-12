# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED
import random
from enum import auto
from enum import Enum
from pathlib import Path

import numpy


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
