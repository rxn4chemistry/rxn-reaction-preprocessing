import click
from tqdm import tqdm
from action_sequences.chemistry.utils import smi_tokenizer, InvalidSmiles
from rdkit.Chem import MolFromSmiles, Mol, MolToSmiles

def randomize_smiles(smiles: str) -> str:
    """
    Return a random SMILES string for a molecule
    """
    # Raise for empty SMILES
    if not smiles:
        raise InvalidSmiles(smiles)

    try:
        return MolToSmiles(MolFromSmiles(smiles), doRandom=True)
    except Exception:
        raise InvalidSmiles(smiles)

def randomize_smiles_with_fragment_bonds(smiles: str, fragment_bond='~') -> str:
    """
    Randomizes a SMILES string that contains fragment bonds
    """
    # Raise for empty SMILES
    if not smiles:
        raise InvalidSmiles(smiles)
    try:
        return '.'.join(
                    sorted([fragment_bond.join([randomize_smiles(fragment) for fragment in group.split(fragment_bond)]) for group in smiles.split('.')]))
    except Exception:
        raise InvalidSmiles(smiles)

@click.command()
@click.option('--input_file', '-i', required=True, help='Path to smiles file (can be tokenized or not)')
@click.option('--output_file', '-o', required=True, help='File to save')
@click.option('--tokenize', '-t', default=True, help='Do tokenization (Default: True')
def main(input_file: str, output_file: str, tokenize: bool) -> None:
    with open (input_file, 'r') as f:
        smiles = [line.strip().replace(' ', '') for line in f.readlines()]
    random_smiles = [randomize_smiles_with_fragment_bonds(smi) for smi in tqdm(smiles)]

    if tokenize:
        random_smiles = [smi_tokenizer(smi) for smi in random_smiles]

    with open(output_file, 'w') as f:

        f.write('\n'.join(random_smiles))

if __name__ =='__main__':
    main()