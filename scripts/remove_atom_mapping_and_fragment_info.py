import click
import pandas as pd
from rxn_chemutils.extended_reaction_smiles import parse_extended_reaction_smiles
from rxn_chemutils.miscellaneous import remove_atom_mapping


@click.command()
@click.option("--input_file_path", "-i", required=True)
@click.option("--output_file_path", "-o", required=True)
@click.option("--rxn_smiles_column", default="rxn")
@click.option("--remove_atom_maps", default=True)
@click.option("--fragment_bond", default="~")
def main(
    input_file_path: str,
    output_file_path: str,
    rxn_smiles_column: str,
    remove_atom_maps: bool,
    fragment_bond: str,
):

    df = pd.read_csv(input_file_path)
    df["original_atom_mapped_rxn"] = df[rxn_smiles_column]
    if remove_atom_maps:
        df[rxn_smiles_column] = df[rxn_smiles_column].apply(
            lambda x: remove_atom_mapping(x)
        )
    df[rxn_smiles_column] = df[rxn_smiles_column].apply(
        lambda x: parse_extended_reaction_smiles(x, remove_atom_maps=False).to_string(
            fragment_bond
        )
    )

    df.to_csv(output_file_path, index=False)


if __name__ == "__main__":
    main()
