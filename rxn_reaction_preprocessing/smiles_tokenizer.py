# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
import re


class SmilesTokenizer:

    def __init__(self):
        # Moved this into a class, as compiling a regex multiple times defeats the purpose
        # of compiling a regex
        self.regex = re.compile(
            r'(\%\([0-9]{3}\)|\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|\%[0-9]{2}|[0-9])'
        )

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

        return ' '.join([token for token in self.regex.findall(smiles)])
