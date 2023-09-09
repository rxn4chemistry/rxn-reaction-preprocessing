# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2022
# ALL RIGHTS RESERVED
from rxn.chemutils.tokenization import tokenize_smiles
from rxn.utilities.csv import iterate_csv_column

from rxn.reaction_preprocessing.config import TokenizeConfig


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


def tokenize(cfg: TokenizeConfig) -> None:
    # Tokenize the reactions
    tokenizer = SmilesTokenizer()
    tokenize_fn = tokenizer.tokenize

    for pair in cfg.input_output_pairs:
        precursors = f"{pair.out}.precursors_tokens"
        products = f"{pair.out}.products_tokens"
        with open(precursors, "wt") as f_precursors, open(products, "wt") as f_products:
            for rxn in iterate_csv_column(pair.inp, pair.reaction_column_name):
                precursors_smiles, products_smiles = rxn.split(">>")
                f_precursors.write(f"{tokenize_fn(precursors_smiles)}\n")
                f_products.write(f"{tokenize_fn(products_smiles)}\n")
