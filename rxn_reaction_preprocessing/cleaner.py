import re


def remove_isotope_information(rxn: str) -> str:
    """
    Function that removes the isotope information from a reaction SMILES. For example [13C] ---> [C].
    """
    return re.sub(
        r'(?<=\[)([0-9]+)(?=[A-Za-z])',  # Remove isotopes
        '',
        rxn.strip()
    )
