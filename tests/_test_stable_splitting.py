from pathlib import Path
import pytest
import pandas as pd
import numpy as np

PATH = str(Path(__file__).parent.parent)


def test_overlap_of_examples():
    with open(PATH + "/data/example1/output/rxns-train.txt") as f:
        rxn_train1 = [line.strip() for line in f.readlines()]
    with open(PATH + "/data/example2/output/rxns-train.txt") as f:
        rxn_train2 = [line.strip() for line in f.readlines()]
    df = pd.DataFrame({"rxn_train1": rxn_train1})

    assert np.all(df.isin(np.array(rxn_train2)).values == True)

    with open(PATH + "/data/example1/output/rxns-test.txt") as f:
        rxn_test1 = [line.strip() for line in f.readlines()]
    with open(PATH + "/data/example2/output/rxns-test.txt") as f:
        rxn_test2 = [line.strip() for line in f.readlines()]
    df = pd.DataFrame({"rxn_test1": rxn_test1})

    assert np.all(df.isin(np.array(rxn_test2)).values == True)

    with open(PATH + "/data/example1/output/rxns-valid.txt") as f:
        rxn_valid1 = [line.strip() for line in f.readlines()]
    with open(PATH + "/data/example2/output/rxns-valid.txt") as f:
        rxn_valid2 = [line.strip() for line in f.readlines()]
    df = pd.DataFrame({"rxn_valid1": rxn_valid1})

    assert np.all(df.isin(np.array(rxn_valid2)).values == True)