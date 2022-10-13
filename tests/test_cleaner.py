from rxn.reaction_preprocessing.cleaner import remove_isotope_information


def test_remove_isotope_information() -> None:
    assert remove_isotope_information("[13CH3][13CH3]") == "[CH3][CH3]"
    assert (
        remove_isotope_information("[13C].CCC.[Na+]~[Fe4+]") == "[C].CCC.[Na+]~[Fe4+]"
    )
    assert remove_isotope_information("[C].CCC.[Na+]~[Fe4+]") == "[C].CCC.[Na+]~[Fe4+]"
