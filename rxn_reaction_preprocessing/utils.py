import random
from pathlib import Path

import numpy


def root_directory() -> Path:
    """
    Returns the path to the root directory of the repository
    """
    return Path(__file__).parent.parent.resolve()


def data_directory() -> Path:
    """
    Returns the path to the data directory at the root of the repository
    """
    return root_directory() / 'data'


def reset_random_seed() -> None:
    random.seed(42)
    numpy.random.seed(42)
