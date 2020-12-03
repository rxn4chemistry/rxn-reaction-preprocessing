# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
import pytest

from data_preprocessor import SmilesTokenizer


@pytest.fixture
def tokenizer():
    return SmilesTokenizer()


def test_tokenize(tokenizer):
    # Should be a more complete / complex test reaction
    assert (
        tokenizer.tokenize('[14C]Cl.O[Na]>O>[Na]Cl.[14C]O') ==
        '[14C] Cl . O [Na] > O > [Na] Cl . [14C] O'
    )
