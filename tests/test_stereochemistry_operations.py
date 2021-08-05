# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED

from rxn_reaction_preprocessing.stereochemistry_operations import remove_chiral_centers


def test_remove_chiral_centers():
    input_expected_dict = {
        'O[C@](Br)(C)N': 'O[C](Br)(C)N',
        'C[C@@](Br)(O)N': 'C[C](Br)(O)N',
        'N[C@H](O)C': 'N[CH](O)C',
        'N1[C@H](Cl)[C@@H](Cl)C(Cl)CC1': 'N1[CH](Cl)[CH](Cl)C(Cl)CC1',
        'F/C=C/F': 'F/C=C/F',
        '[N@+]': '[N+]',
        '[N@@+]': '[N+]',
        '[Si@@]': '[Si]',
        '[Si@H]': '[SiH]'
    }

    assert all(
        [remove_chiral_centers(smi) == expected for smi, expected in input_expected_dict.items()]
    )
