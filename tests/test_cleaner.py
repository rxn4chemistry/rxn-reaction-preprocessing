from rxn_reaction_preprocessing.cleaner import remove_isotope_information


def test_remove_isotope_information():
    assert remove_isotope_information('[13C].CCC.[Na+]~[Fe4+]') == '[C].CCC.[Na+]~[Fe4+]'
    assert remove_isotope_information('[C].CCC.[Na+]~[Fe4+]') == '[C].CCC.[Na+]~[Fe4+]'
