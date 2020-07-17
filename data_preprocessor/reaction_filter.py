import re
from rdkit import Chem
from action_sequences.chemistry.utils import tokenize_smiles


class MixedReactionFilter():
    def __init__(self, rxn,
                        maximum_number_of_precursors=10,
                        maximum_number_of_products=1,
                        maximum_number_of_precursor_tokens=300,
                        maximum_number_of_product_tokens=200,
                        max_absolute_formal_charge=2,
                        filter_single_atom_products=True,
                        filter_purifications=True,
                        filter_alchemy=True,
                        fragment_bond='.'
                        ): # atoms not in precursors
        super().__init__()
        self.precursors_smi, self.product_smi = rxn.split('>>')
        self.maximum_number_of_products = maximum_number_of_products
        self.maximum_number_of_precursors = maximum_number_of_precursors
        self.maximum_number_of_product_tokens = maximum_number_of_product_tokens
        self.maximum_number_of_precursor_tokens = maximum_number_of_precursor_tokens
        self.max_absolute_formal_charge = max_absolute_formal_charge

        self.precursor_mol = None
        self.product_mol = None
        self.fragment_bond=fragment_bond

    def apply_max_number_of_products(self):
        if not self.less_or_equal_N_molecules(self.product_smi, self.maximum_number_of_products):
            return False
        return True

    def apply_max_number_of_precursors(self):
        if not self.less_or_equal_N_molecules(self.precursors_smi, self.maximum_number_of_precursors):
            return False
        return True

    def apply_all_products_in_precursors(self):
        if self.all_products_in_precursors():
            return False
        return True

    def apply_product_is_single_atom(self):
        if self.product_is_single_atom():
            return False
        return True

    def apply_max_number_of_tokens(self):
        if len(tokenize_smiles(self.product_smi)) > self.maximum_number_of_product_tokens:
            return False
        elif len(tokenize_smiles(self.precursors_smi)) > self.maximum_number_of_precursor_tokens:
            return False
        return True

    def apply_formal_charge_precursors_and_products(self):

        self.precursor_mol = Chem.MolFromSmiles(self.precursors_smi.replace(self.fragment_bond, '.'))
        precursor_charge = self.get_charge(self.precursor_mol)
        if (precursor_charge is None) or (abs(precursor_charge) > self.max_absolute_formal_charge):
            return False

        self.product_mol = Chem.MolFromSmiles(self.product_smi.replace(self.fragment_bond, '.'))
        product_charge = self.get_charge(self.product_mol)
        if (product_charge is None) or (abs(product_charge) > self.max_absolute_formal_charge):
            return False
        return True

    def apply_same_atom_types(self):
        if not self.products_and_reactants_same_atom_types():
            return False
        return True

    def apply_all_filters(self):
        """Return True if all filters are passed
        """

        if not self.less_or_equal_N_molecules(self.product_smi, self.maximum_number_of_products):
            return False
        elif not self.less_or_equal_N_molecules(self.precursors_smi, self.maximum_number_of_precursors):
            return False
        elif self.product_in_precursors():
            return False
        elif self.product_is_single_atom():
            return False
        elif len(tokenize_smiles(self.product_smi)) > self.maximum_number_of_product_tokens:
            return False
        elif len(tokenize_smiles(self.precursors_smi)) > self.maximum_number_of_precursor_tokens:
            return False

        self.precursor_mol = Chem.MolFromSmiles(self.precursors_smi.replace(self.fragment_bond, '.'))

        precursor_charge = self.get_charge(self.precursor_mol)
        if (precursor_charge is None) or (abs(precursor_charge) > self.max_absolute_formal_charge):
            return False

        self.product_mol = Chem.MolFromSmiles(self.product_smi.replace(self.fragment_bond, '.'))
        product_charge = self.get_charge(self.product_mol)

        if (product_charge is None) or (abs(product_charge) > self.max_absolute_formal_charge):
            return False

        if not self.products_and_reactants_same_atom_types():
            return False
        return True


    def less_or_equal_N_molecules(self, smi, N):
        return len(smi.split('.')) <= N

    def all_products_in_precursors(self):
        for elem in self.product_smi.split('.'):
            if elem not in self.precursors_smi.split('.'):
                return False
        return True

    def product_in_precursors(self):
        return self.product_smi in self.precursors_smi.split('.')

    def product_is_single_atom(self):
        single_atom_regex = re.compile(r'^\[\w{1,2}[+-]\d?\]$|^\w{1,2}$')
        match = re.match(single_atom_regex, self.product_smi)
        if match is not None:
            return True
        else:
            return False

    def get_charge(self, mol):
        if mol is None:
            return None
        try:
            return Chem.GetFormalCharge(mol)
        except:
            return None

    # Check if all product atom types in reactants (+reagents)
    def products_and_reactants_same_atom_types(self):
        if self.product_mol is None:
            self.product_mol = Chem.MolFromSmiles(self.product_smi.replace(self.fragment_bond, '.'))
        if self.precursor_mol is None:
            self.precursor_mol = Chem.MolFromSmiles(self.precursors_smi.replace(self.fragment_bond, '.'))

        product_atoms = set([_.GetSymbol() for _ in self.product_mol.GetAtoms()])
        precursor_atoms = set([_.GetSymbol() for _ in self.precursor_mol.GetAtoms()])

        product_atoms -= precursor_atoms

        return len(product_atoms) == 0


if __name__ == '__main__':
    rxn = 'O=[N+]([O-])c1cc(-c2nc3ccccc3o2)ccc1F.C~C>>Nc1cc(-c2nc3ccccc3o2)ccc1NCC(=O)N1CCOCC1'

    rxn_filter = MixedReactionFilter(rxn, fragment_bond='~')

    passes_filter = rxn_filter.apply_all_filters()

    print(rxn, passes_filter)