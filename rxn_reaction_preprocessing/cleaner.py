import re

_ISOTOPE_REMOVAL_REGEX = re.compile(r"(?<=\[)([0-9]+)(?=[A-Za-z])")


def remove_isotope_information(rxn: str) -> str:
    """
    Function that removes the isotope information from a reaction SMILES.

    For example [13CH3][13CH3] ---> [CH3][CH3].
    """
    return _ISOTOPE_REMOVAL_REGEX.sub("", rxn.strip())
