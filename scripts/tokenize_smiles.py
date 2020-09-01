from action_sequences.chemistry.utils import *
import click

@click.command()
@click.option('--input_file', '-i', required=True, help='Path to smiles file (can be tokenized or not)')
@click.option('--output_file', '-o', required=True, help='File to save')
def main(input_file: str, output_file: str) -> None:
    with open (input_file, 'r') as f:
        smiles = [line.strip().replace(' ', '') for line in f.readlines()]

    tok_smile = [smi_tokenizer(smi) for smi in smiles]

    with open(output_file, 'w') as f:
        f.write('\n'.join(tok_smile))

if __name__ =='__main__':
    main()
