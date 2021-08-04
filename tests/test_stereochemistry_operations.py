# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED

from rxn_reaction_preprocessing.stereochemistry_operations import remove_chiral_centers


def test_remove_chiral_centers():
    input_expected_dict = {
        'O[C@](Br)(C)N': 'OC(Br)(C)N',
        'C[C@@](Br)(O)N': 'CC(Br)(O)N',
        'N[C@H](O)C': 'NC(O)C',
        'N1[C@H](Cl)[C@@H](Cl)C(Cl)CC1': 'N1C(Cl)C(Cl)C(Cl)CC1',
        'F/C=C/F': 'F/C=C/F'
    }

    assert all(
        [remove_chiral_centers(smi) == expected for smi, expected in input_expected_dict.items()]
    )
