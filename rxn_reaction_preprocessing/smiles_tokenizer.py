# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED

from rxn_chemutils.tokenization import tokenize_smiles


class SmilesTokenizer:

    #
    # Public Methods
    #

    def tokenize(self, smiles: str) -> str:
        """
        Tokenize a SMILES molecule or reaction, and join the tokens with spaces.
        Args:
            smiles: SMILES string to tokenize, for instance 'CC(CO)=N>>CC(C=O)N'.
        Returns:
            SMILES string after tokenization, for instance 'C C ( C O ) = N >> C C ( C = O ) N'.
        """

        return tokenize_smiles(smiles)
