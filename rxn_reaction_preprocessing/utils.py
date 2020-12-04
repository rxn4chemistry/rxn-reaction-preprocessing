from pathlib import Path


class InvalidSmiles(ValueError):

    def __init__(self, smiles: str):
        super().__init__(f'"{smiles}" is not a valid SMILES string')


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
