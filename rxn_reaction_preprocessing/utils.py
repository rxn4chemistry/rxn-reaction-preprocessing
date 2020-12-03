from pathlib import Path


class InvalidSmiles(ValueError):

    def __init__(self, smiles: str):
        super().__init__(f'"{smiles}" is not a valid SMILES string')


def root_directory() -> Path:
    """
    Returns the path to the root directory of the repository
    """
    this_dir = Path(__file__).resolve().parent
    data_dir = this_dir / '..'
    return data_dir.resolve()


def data_directory() -> Path:
    """
    Returns the path to the data directory at the root of the repository
    """
    data_dir = root_directory() / 'data'
    return data_dir.resolve()
