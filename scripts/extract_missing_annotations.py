import ast
from collections import Counter

import click
import pandas as pd
from rdkit import RDLogger
from rxn_chemutils.conversion import canonicalize_smiles

RDLogger.DisableLog("rdApp.*")


@click.command()
@click.argument("input_csv", type=str, required=True)
@click.argument("output_csv", type=str, required=True)
@click.option("--threshold-counts", "-tc", type=int, default=5)
def main(input_csv: str, output_csv: str, threshold_counts: int):
    """
    Script to count and list the number of molecules that still need an annotation.
    The input csv file requires 2 columns:
        'original_atom_mapped_rxn': contains the reaction SMILES and it is used to remove duplicates
        'rxn_missing_annotations': for each line of the csv file a list of molecules requiring annotations
    The two columns are automatically generated with the standardization script
    """
    df = pd.read_csv(input_csv)
    print(df.head())

    print(f"Length of dataset before duplicates removal: {len(df)}")
    df.drop_duplicates(["original_atom_mapped_rxn"], inplace=True)
    print(f"Length of dataset after duplicates removal: {len(df)}")
    # converting a "list string" to a list
    df.rxn_missing_annotations = df.rxn_missing_annotations.apply(
        lambda x: ast.literal_eval(x)
    )

    # save a list of canonical missing annotations
    missing_annotations = []
    for elem in df.rxn_missing_annotations.values:
        if elem != []:
            missing_annotations.extend([canonicalize_smiles(smi) for smi in elem])

    c = Counter(missing_annotations)
    print(c.most_common(10))
    print(f"Number of missing annotations: {len(c)}")

    catalysts = []
    counts = []

    for key, count in c.items():
        if count >= threshold_counts:
            catalysts.append(key)
            counts.append(count)
    print(
        f"Number of missing annotations with count >={threshold_counts}: {len(catalysts)}"
    )

    with open(output_csv, "w") as f:
        f.write("catalysts,counts\n")
        for x in zip(catalysts, counts):
            f.write(f"{x[0]},{x[1]}\n")


if __name__ == "__main__":
    main()
