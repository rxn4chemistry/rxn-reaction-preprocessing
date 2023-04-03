# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2022
# ALL RIGHTS RESERVED
import pandas as pd
from rxn.chemutils.tokenization import tokenize_smiles

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

    for pair in cfg.input_output_pairs:
        df = pd.read_csv(pair.inp, lineterminator="\n")
        if pair.reaction_column_name not in df.columns:
            raise SystemExit(
                f"The following file does not contain an rxn column:\n{pair.inp}"
            )

        df["rxn_precursors"] = df[pair.reaction_column_name].str.split(">>").str[0]
        df["rxn_products"] = df[pair.reaction_column_name].str.split(">>").str[1]

        df.rxn_precursors = df.rxn_precursors.apply(tokenizer.tokenize)
        df.rxn_products = df.rxn_products.apply(tokenizer.tokenize)

        df[["rxn_precursors"]].to_csv(
            f"{pair.out}.precursors_tokens", header=False, index=False
        )
        df[["rxn_products"]].to_csv(
            f"{pair.out}.products_tokens", header=False, index=False
        )
