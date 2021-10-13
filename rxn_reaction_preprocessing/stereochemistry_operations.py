# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
""" Contains functions to handle stereochemistry in molecules. """

import re

CHIRAL_CENTER_PATTERN = re.compile(
    r'\[([^],@]+)@+([^]]*)]'
)  # Matches stereo centres, and groups what comes before @


def remove_chiral_centers(smiles: str) -> str:
    """
    Return smiles where all chiral centres are removed.
    Args:
        smiles: non-atom-mapped smiles string.
    """
    return re.sub(CHIRAL_CENTER_PATTERN, r'[\g<1>\g<2>]', smiles)
