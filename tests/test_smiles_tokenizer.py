# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2022
# ALL RIGHTS RESERVED
import pytest

from rxn.reaction_preprocessing import SmilesTokenizer


@pytest.fixture
def tokenizer() -> SmilesTokenizer:
    return SmilesTokenizer()


def test_tokenize(tokenizer: SmilesTokenizer) -> None:
    # Should be a more complete / complex test reaction
    assert (
        tokenizer.tokenize("[14C]Cl.O[Na]>O>[Na]Cl.[14C]O")
        == "[14C] Cl . O [Na] > O > [Na] Cl . [14C] O"
    )
