"""
Utilities to extract the data from the original file of the pistachio dataset
and apply some coarse filtering
"""

import os
import pandas as pd
from tqdm import tqdm
import numpy as np
import random
import joblib as jl
from reaction_smiles_extractor import ReactionSmilesExtractor
from data_preprocessor.reaction_filter import MixedReactionFilter
from .tokenization import tokenize_smiles
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import Chem

from crc64iso.crc64iso import crc64

"""
In the final version of the dataframe saved as a pkl extension as 'dump_name_final.pkl' the following informations are available

'US_classes': patent name
'classes_ids': full id of the class
'classes_names': name of the reaction class
'reactions_mixed': reaction smile, fragmented
'precursors': smile of the precursors set with just one appearance of ech molecule
'precursors_multi': smile of the precursors set with possible multiple apperances of mulecules
'products': smile of the products set with possible multiple appearance of mulecules
'has_fragments': is True or False depending if rxns have the fragment bond ~
'tokenized_precursors': tokenized version of 'precursors_multi'
'tokenized_products': tokenized version of 'products'
'largest_fragment_products': column with only the largest fragment within prroducts
'tokenized_largest_fragment_products': tokenized version

"""


def hash_key(identifier):
    return int(crc64(identifier), 16)


def split_smile(smi):
    """
    Checks the validity of the smiles string, after having removed the atom mapping
    :param smi: the smile string
    :return: the splitted smile string
    """
    try:
        RF = ReactionSmilesExtractor(smi)
        splittedsmi = RF.get_fragmented_reaction_smiles(
            mixed=True, order_molecules=True, fragment_bond="~"
        )
    except:
        splittedsmi = " "
    return splittedsmi


def passed_filter(x):
    """
    :param x:
    :return: bool, is False if at least one filter has not passed
    """
    rxn_filter = MixedReactionFilter(x, maximum_number_of_products=4, fragment_bond="~")
    result, _ = rxn_filter.apply_all_filters()
    return result


def filter_products_number(x):
    rxn_filter = MixedReactionFilter(x, maximum_number_of_products=4, fragment_bond="~")
    return rxn_filter.apply_max_number_of_products()


def filter_precursors_number(x):
    rxn_filter = MixedReactionFilter(x, maximum_number_of_products=4, fragment_bond="~")
    return rxn_filter.apply_max_number_of_precursors()


def filter_single_atom(x):
    rxn_filter = MixedReactionFilter(x, maximum_number_of_products=4, fragment_bond="~")
    return rxn_filter.apply_product_is_single_atom()


def filter_all_products_in_precursors(x):
    rxn_filter = MixedReactionFilter(x, maximum_number_of_products=4, fragment_bond="~")
    return rxn_filter.apply_all_products_in_precursors()


def filter_tokens_number(x):
    rxn_filter = MixedReactionFilter(x, maximum_number_of_products=4, fragment_bond="~")
    return rxn_filter.apply_max_number_of_tokens()


def filter_formal_charge(x):
    rxn_filter = MixedReactionFilter(x, maximum_number_of_products=4, fragment_bond="~")
    return rxn_filter.apply_formal_charge_precursors_and_products()


def filter_atom_types(x):
    rxn_filter = MixedReactionFilter(x, maximum_number_of_products=4, fragment_bond="~")
    return rxn_filter.apply_same_atom_types()


# Check that it is not a purification step and in that case remove the reaction
def is_in_precursors(row):
    reactants = row["tokenized_precursors"].replace(" ", "").split(".")
    return row["tokenized_largest_fragment_products"] in reactants


def remove_res_reagents(x):
    products = x[1].split(".")
    precursors = x[0].split(".")
    for j, elem in enumerate(products):
        if elem in precursors:
            products.pop(j)
    return ".".join(products)


class DataExtractionCleaning:
    def __init__(self, datapath, destpath, filenames=[], df=None, atom_mapping=False):
        self.datapath = datapath
        self.destpath = destpath
        if not isinstance(filenames, list):
            self.filenames = [filenames]
        else:
            self.filenames = filenames
        self.df = df
        self.atom_mapping = atom_mapping

    @classmethod
    def from_df(cls, datapath, destpath, df_csv_filename, atom_mapping):
        """
        Alternative constructor from dataset stored in pkl or csv file
        :param datapath: the path of the file to load
        :param destpath: the destination path of the dumping
        :param df_csv_filename: the name of the file to load
        :return:
        """
        df = pd.read_csv(destpath + df_csv_filename)
        return cls(datapath, destpath, df=df, atom_mapping=atom_mapping)

    def read_data(self):
        print("Loading dataset ...")

        for f in self.filenames:
            new_df = pd.read_csv(self.datapath + f, sep="\t")
            print(type(self.df))
            self.df = pd.concat([self.df, new_df], axis=1, sort=False)

        if self.atom_mapping:
            print(type(self.atom_mapping))
            self.df.columns = [
                "reactions",
                "classes_long",
                "atom_mapping",
                "classes_long_cp",
            ]
            assert list(self.df.classes_long.values) == list(
                self.df.classes_long_cp.values
            )
            self.df = self.df.drop(columns=["classes_long_cp"])
        else:
            self.df.columns = ["reactions", "classes_long"]

        print(self.df.loc[0].values)
        print("Finished loading dataset. Dimension: {}".format(len(self.df)))

        return

    def remove_duplicates(self, cols):
        print("Removing rxn duplicates...")
        self.df.drop_duplicates(subset=cols, keep="first", inplace=True)
        print("Removed duplicates, length of dataset {}".format(len(self.df)))
        return

    def split_class_info(self):
        print("Splitting class information ...")
        self.df["classes_long"] = self.df["classes_long"].apply(
            lambda x: [l.strip() for l in x.split()]
        )
        self.df["US_classes"] = self.df["classes_long"].apply(lambda x: x[0])
        self.df["classes_ids"] = self.df["classes_long"].apply(
            lambda x: x[1].replace("[", "").replace("]", "")
        )
        self.df["classes_names"] = self.df["classes_long"].apply(
            lambda x: " ".join(x[2:])
        )
        print("Dropping column with whole class information ...")
        self.df = self.df.drop(columns=["classes_long"])
        return

    def get_fragmented_rxn_smiles(self):
        print("Getting fragmented rxn smiles ...")
        self.df["reactions_mixed"] = [
            split_smile(rxn) for rxn in tqdm(self.df.reactions)
        ]
        print("Removing 'reactions' column and generating 'reactions_mixed'")
        self.df = self.df.drop(columns=["reactions"])
        self.df = self.df.loc[self.df["reactions_mixed"] != " "]
        print("Removed invalid smiles, length of dataset {}".format(len(self.df)))
        return

    def filter_dataset(self, filtername):
        fname = str(filtername).split(" ")[1]
        print(f"Applying filter '{fname}'...")
        self.df["filtering"] = [filtername(x) for x in tqdm(self.df.reactions_mixed)]
        print(f"Dumping rxns which did not pass the filter...")
        self.df.loc[self.df["filtering"] == False].to_csv(
            self.destpath + f"/dump_step_2_{fname}_removed.csv"
        )
        print("Removing rxns which did not pass the filter ...")
        self.df = self.df.loc[self.df["filtering"] == True]
        self.df.to_csv(self.destpath + f"/dump_step_2_{fname}.csv")
        self.df = self.df.drop(columns=["filtering"])
        print("Size of database: {}".format(len(self.df)))
        return

    def split_precursors_products(self):
        print("Splitting the rxns into the column 'splitted_reactions_mixed' ...")
        self.df["splitted_reactions_mixed"] = [
            x.split(">>") for x in self.df.reactions_mixed
        ]
        print(
            "Putting the precursors with multiple appearances in 'precursors_multi' ..."
        )
        self.df["precursors_multi"] = [x[0] for x in self.df.splitted_reactions_mixed]
        print(
            "Putting the products with multiple appearances in 'precursors_multi' ..."
        )
        self.df["products_multi"] = [x[1] for x in self.df.splitted_reactions_mixed]
        return

    def mark_row_with_a_fragment(self):
        print("Marking the rows with a fragment in the products ...")
        self.df["has_fragments"] = self.df.products.apply(lambda x: "~" in x)
        return

    def tokenize(self):
        print("Tokenizing the SMILEs ...")
        self.df["tokenized_precursors"] = [
            tokenize_smiles(x) for x in self.df.precursors
        ]
        # self.df.drop(column = ['precursors'])
        self.df["tokenized_products"] = [tokenize_smiles(x) for x in self.df.products]
        self.df["tokenized_reactions_mixed"] = (
            self.df.tokenized_precursors + " >> " + self.df.tokenized_products
        )
        return

    def keep_largest_fragment(self):
        lfc = rdMolStandardize.LargestFragmentChooser()

        self.df["largest_fragment_products"] = [
            Chem.MolToSmiles(lfc.choose(Chem.MolFromSmiles(smi.replace("~", "."))))
            for smi in tqdm(self.df.products)
        ]

        print("Size od dataset: {}".format(len(self.df)))
        self.df["tokenized_largest_fragment_products"] = [
            tokenize_smiles(x) for x in self.df.largest_fragment_products
        ]
        print("Kept largest fragment in product")
        return

    def remove_purification(self):
        self.df["is_in_precursors"] = [
            is_in_precursors(row) for idx, row in tqdm(self.df.iterrows())
        ]
        print("Size od dataset: {}".format(len(self.df)))
        self.df = self.df.loc[self.df["is_in_precursors"] == False]
        self.df = self.df.drop(columns=["is_in_precursors"])
        print("Filtered purification rxns")
        print("Size od dataset: {}".format(len(self.df)))
        return

    def do_set_on_precursors_and_products(self):
        # the column 'precursors_multi' is the one that keeps multiple appearances of a molecule
        print("Performing set on precursors and products ...")
        self.df["precursors"] = [
            ".".join(sorted(list(set(x.split(".")))))
            for x in self.df["precursors_multi"].values
        ]
        self.df["products"] = [
            ".".join(sorted(list(set(x.split(".")))))
            for x in self.df["products_multi"].values
        ]
        print("Generating rxns with precursors and products after 'set'...")
        self.df["reactions_mixed"] = self.df.precursors + ">>" + self.df.products
        return

    def remove_residual_reagents(self):
        print("Removing residual reagents from products ...")
        self.df["products_no_residual"] = [
            remove_res_reagents(x)
            for x in tqdm(self.df[["precursors", "products"]].values)
        ]
        self.df["products"] = self.df["products_no_residual"]
        self.df = self.df.drop(columns=["products_no_residual"])
        self.df = self.df.loc[self.df["products"] != ""]
        self.df["reactions_mixed"] = self.df.precursors + ">>" + self.df.products
        print("Size od dataset: {}".format(len(self.df)))

    def remove_single_precursor_rxns(self):
        print("Removing sigle precursor rxns (except cycloaddition) ...")
        saveset = ["cycloaddition"]

        self.df["filtering"] = [len(x.split(".")) > 1 for x in self.df.precursors]

        self.df["splitted_classes_names"] = [
            x.split(" ")[-1] for x in self.df.classes_names
        ]
        print(self.df[["splitted_classes_names", "classes_names"]].loc[0:10])
        self.df.loc[
            (self.df["splitted_classes_names"] == saveset[0])
            & (self.df["filtering"] == False),
            "filtering",
        ] = True
        self.df.loc[self.df["filtering"] == False].to_csv(
            self.destpath + f"/dump_step_4_one_precursor_removed.csv"
        )
        self.df.loc[
            (self.df["splitted_classes_names"] == saveset[0])
            & (self.df["filtering"] == False)
        ].to_csv(self.destpath + f"/dump_step_4_one_precursor_{saveset[0]}_kept.csv")

        self.df = self.df.loc[self.df["filtering"] == True]
        self.df = self.df.drop(columns=["filtering", "splitted_classes_names"])
        print("Size of dataset: {}".format(len(self.df)))

    def split_train_test_by_id(self, split_ratio, id_column):
        self.df["crc_values"] = self.df[id_column].apply(lambda x: hash_key(x))

        in_test_set = self.df["crc_values"].apply(
            lambda id_: id_ < split_ratio * 2 ** 64
        )
        in_valid_set = self.df["crc_values"].apply(
            lambda id_: (split_ratio * 2 ** 64 <= id_)
            and (id_ < split_ratio * 2 * 2 ** 64)
        )
        in_train_set = self.df["crc_values"].apply(
            lambda id_: id_ >= split_ratio * 2 * 2 ** 64
        )

        return (
            self.df.loc[in_train_set],
            self.df.loc[in_valid_set],
            self.df.loc[in_test_set],
        )

    def generate_training_files_stable(self, split_ratio):
        train, valid, test = self.split_train_test_by_id(split_ratio, "products")
        return train, valid, test

    def generate_training_files_random(self, split_ratio):
        random.seed(42)
        np.random.seed(42)
        number_of_unique_products = len(self.df["tokenized_products"].unique())
        sample = np.split(
            pd.Series(self.df["tokenized_products"].unique()).sample(
                frac=1, random_state=42
            ),
            [
                int((1.0 - split_ratio * 2) * number_of_unique_products),
                int((1.0 - split_ratio) * number_of_unique_products),
            ],
        )

        train = self.df.loc[self.df["tokenized_products"].isin(sample[0])]
        valid = self.df.loc[self.df["tokenized_products"].isin(sample[1])]
        test = self.df.loc[self.df["tokenized_products"].isin(sample[2])]

        return train, valid, test

    def generate_training_files(self, tipo: str = "stable", split_ratio=0.05):

        if tipo == "stable":
            train, valid, test = self.generate_training_files_stable(split_ratio)
        elif tipo == "random":
            train, valid, test = self.generate_training_files_random(split_ratio)
        else:
            raise ValueError(f"{type} not implemented! ")

        train = train.sample(frac=1, random_state=42)
        valid = valid.sample(frac=1, random_state=42)
        test = test.sample(frac=1, random_state=42)

        print(len(train), len(valid), len(test))
        print("First lines of the test dataset:\n")
        print(test[:5])

        print("Generating train, test, validation files ...")
        for name, data in tqdm([("train", train), ("valid", valid), ("test", test)]):
            with open(self.destpath + "precursors-{}.txt".format(name), "w") as f:
                f.write("\n".join(data.tokenized_precursors.values))
            with open(self.destpath + "product-{}.txt".format(name), "w") as f:
                f.write("\n".join(data.tokenized_products.values))

        # Generate the needed files for the classes
        print("Generating classes files ...")
        with open(self.destpath + "class-single-train.txt", "w") as f:
            f.write("\n".join(train.classes_ids.values))
        with open(self.destpath + "class-multi-train.txt", "w") as f:
            classes_splitted = [
                (x.split(".")[0], x.rsplit(".", 1)[0], x)
                for x in train.classes_ids.values
            ]
            incremental_classes = [" ".join(x) for x in classes_splitted]
            f.write("\n".join(incremental_classes))

        with open(self.destpath + "class-single-test.txt", "w") as f:
            f.write("\n".join(test.classes_ids.values))
        with open(self.destpath + "class-multi-test.txt", "w") as f:
            classes_splitted = [
                (x.split(".")[0], x.rsplit(".", 1)[0], x)
                for x in test.classes_ids.values
            ]
            incremental_classes = [" ".join(x) for x in classes_splitted]
            f.write("\n".join(incremental_classes))

        with open(self.destpath + "class-single-valid.txt", "w") as f:
            f.write("\n".join(valid.classes_ids.values))
        with open(self.destpath + "class-multi-valid.txt", "w") as f:
            classes_splitted = [
                (x.split(".")[0], x.rsplit(".", 1)[0], x)
                for x in valid.classes_ids.values
            ]
            incremental_classes = [" ".join(x) for x in classes_splitted]
            f.write("\n".join(incremental_classes))

        print("Generating rxns files ...")
        # Save full rxns
        with open(self.destpath + "rxns-train.txt", "w") as f:
            f.write("\n".join(train.tokenized_reactions_mixed.values))
        with open(self.destpath + "rxns-test.txt", "w") as f:
            f.write("\n".join(test.tokenized_reactions_mixed.values))
        with open(self.destpath + "rxns-valid.txt", "w") as f:
            f.write("\n".join(valid.tokenized_reactions_mixed.values))

        print("Generating patent files ...")
        # Save patents
        with open(self.destpath + "rxn-patents-train.txt", "w") as f:
            f.write("\n".join(train.US_classes.values))
        with open(self.destpath + "rxn-patents-test.txt", "w") as f:
            f.write("\n".join(test.US_classes.values))
        with open(self.destpath + "rxn-patents-valid.txt", "w") as f:
            f.write("\n".join(valid.US_classes.values))

        print("Generating rxn name files ...")

        # Save class name
        with open(self.destpath + "rxn-name-train.txt", "w") as f:
            f.write("\n".join(train.classes_names.values))
        with open(self.destpath + "rxn-name-test.txt", "w") as f:
            f.write("\n".join(test.classes_names.values))
        with open(self.destpath + "rxn-name-valid.txt", "w") as f:
            f.write("\n".join(valid.classes_names.values))

        print("Generating atom-mapping files ...")

        # Save atom mapping
        if self.atom_mapping:
            train.atom_mapping.to_csv(self.destpath + "rxn-atom-mapping-train.csv")
            test.atom_mapping.to_csv(self.destpath + "rxn-atom-mapping-test.csv")
            test.atom_mapping.to_csv(self.destpath + "rxn-atom-mapping-valid.csv")

        """
        with open(self.destpath + 'rxn-atom-mapping-train.txt','w') as f:
            f.write('\n'.join(train.atom_mapping.values))
        with open(self.destpath + 'rxn-atom-mapping-test.txt','w') as f:
            f.write('\n'.join(test.atom_mapping.values))
        with open(self.destpath + 'rxn-atom-mapping-valid.txt','w') as f:
            f.write('\n'.join(valid.atom_mapping.values))
        """


if __name__ == "__main__":
    print("No errors for now!")
