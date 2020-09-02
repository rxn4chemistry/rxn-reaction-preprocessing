import click
from data_preprocessor.PistFilter import *
from typing import List

@click.command()
@click.option('--datapath', '-i', required=True, help='Path pistachio csv file')
@click.option('--destpath', '-o', required=True, help='Destination path')
@click.option('--filenames', '-f', required=True, default=['pistachio.cansmi.reduced','pistachio.smi.reduced'],
              help='List of filenames to add to the filtering')
@click.option('--atom_mapping','-am', is_flag=True, required=True, default=False, help='Boolean to set the presence of atom mapped reactions')
def main(datapath: str, destpath: str, filenames: List[str], atom_mapping : bool) -> None:

    print("Data path: ", datapath)
    print("Destination path: ", destpath)
    print("filenames: ", filenames)
    print("Atom mapping? ", type(atom_mapping))

    DE = DataExtractionCleaning(datapath=datapath, destpath=destpath, filenames=filenames, df=None, atom_mapping=atom_mapping)
    DE.read_data()
    DE.remove_duplicates(['reactions'])
    DE.df.to_csv(DE.destpath + 'dump_step_0.csv')
    DE.split_class_info()

    DE.get_fragmented_rxn_smiles()
    DE.remove_duplicates(['reactions_mixed'])

    DE.df.to_csv(DE.destpath + 'dump_step_1.csv')
    DE.split_precursors_products()

    DE.do_set_on_precursors_and_products()
    DE.remove_duplicates(['precursors', 'products'])
    DE.filter_dataset(filtername=filter_products_number)
    DE.filter_dataset(filtername=filter_precursors_number)
    DE.filter_dataset(filtername=filter_single_atom)
    DE.filter_dataset(filtername=filter_all_products_in_precursors)
    DE.filter_dataset(filtername=filter_tokens_number)
    DE.filter_dataset(filtername=filter_formal_charge)
    DE.filter_dataset(filtername=filter_atom_types)
    DE.remove_residual_reagents()
    DE.df.to_csv(DE.destpath + 'dump_step_3.csv')
    DE.remove_single_precursor_rxns()
    DE.df.to_csv(DE.destpath + 'dump_step_4.csv')
    DE.tokenize()
    DE.mark_row_with_a_fragment()
    DE.df.to_csv(DE.destpath + 'dump_step_final.csv')

    # df_csv_filename = '/dump_step_final.csv'
    # DE = DataExtractionCleaning.from_df(datapath, destpath, df_csv_filename)
    # DE.df['reactions_mixed'] = DE.df.precursors + ">>" + DE.df.products
    # DE.df['tokenized_reactions_mixed'] = DE.df.tokenized_precursors + " >> " + DE.df.tokenized_products
    DE.generate_training_files(tipo='stable', split_ratio=0.05)

if __name__ == "__main__":
    main()
