
class InvalidSmiles(ValueError):

    def __init__(self, smiles: str):
        super().__init__(f'"{smiles}" is not a valid SMILES string')