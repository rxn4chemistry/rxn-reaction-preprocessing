# rxn-reaction-preprocessing

Repository to apply many preprocessing transformations including stable train/test/validation splits.

## Install / Build
### Create Environment and Install
```bash
conda create -c rdkit -n rxn-reaction-preprocessing rdkit python=3.6
pip install -e .
```

### Development
```bash
pip install -e .[dev]
```
Escape brackets in zsh
```bash
pip install -e .\[dev\]
```

### Running tests
```bash
python -m pytest
```

## Before you commit
Please apply pre-commit before commiting. This makes life easier for everyone (including you, as chances for CI failing will be smaller afterwards).
```bash
pre-commit install
```

## Usage
The following command line scripts are installed with the package

### rxn-standardize
Usage:
```
rxn-standardize [OPTIONS] INPUT OUTPUT

  The entry point for this cli script.

  Args:     input (str):  The input file path (at least one column with one reaction SMILES per
            line). The column, if not the only one, should be named `rxn`.
            output (str): The output file path.

Options:
  --fragment_bond STRING [default='.']
  --help                          Show this message and exit.
```
### rxn-preprocess
Usage:
```
rxn-preprocess [OPTIONS] INPUT OUTPUT

  The entry point for this cli script.

  Args:     input (str):  The input file path (at least one column with one reaction SMILES per
            line). The column, if not the only one, should be named `rxn`.
            output (str): The output file path.

Options:
  --fragment_bond STRING [default='.']
  --max-reactants INTEGER
  --max-agents INTEGER
  --max-products INTEGER
  --min-reactants INTEGER
  --min-agents INTEGER
  --min-products INTEGER
  --max-reactants_tokens INTEGER
  --max-agents_tokens INTEGER
  --max-products-tokens INTEGER
  --max-absolute-formal-charge INTEGER
  --help                          Show this message and exit.
```
### rxn-split
Usage:
```
rxn-split [OPTIONS] INPUT OUTPUT

  The entry point for this cli script.

  Args:     input (str):  The input file path.     output (str): The output
  directory (where train, validation, and test files will be written to).

Options:
  --split-ratio FLOAT
  --index-column STRING
```
### rxn-tokenize
Usage:
```
rxn-tokenize [OPTIONS] INPUT_OUTPUT_PAIRS...

  The entry point for this cli script.

  Args:     input_output_pairs (Tuple[str]):  Paths to the input and output
  files in the form <input_a> <output_a> <input_b> <output_b> [...].

Options:
  --help  Show this message and exit.
```

## Running in Docker
The package is readily dockerized using the supplied ```Dockerfile```. 
It is necessary to provide a GitHub token as the `GHE_TOKEN_ARG` argument.
```bash
docker build --build-arg GHE_TOKEN_ARG=${GHE_TOKEN} -t rxn_reaction_preprocessing .
```

`demo.sh` contains an example of running the preprocessing pipeline using the resulting docker image
```bash
#!/bin/sh

docker run \
--mount type=bind,source=/mnt/c/Users/DanielProbst/Code/rxn_reaction_preprocessing/data,target=/data \
uk.icr.io/rxn-test/rxn-reaction-preprocessing \
rxn-preprocess \
/data/example/raw/rxns-small.txt \
/data/example/processed/rxns-small.processed.csv

docker run \
--mount type=bind,source=/mnt/c/Users/DanielProbst/Code/rxn_reaction_preprocessing/data,target=/data \
uk.icr.io/rxn-test/rxn-reaction-preprocessing \
rxn-split \
/data/example/processed/rxns-small.processed.csv \
/data/example/processed

docker run \
--mount type=bind,source=/mnt/c/Users/DanielProbst/Code/rxn_reaction_preprocessing/data,target=/data \
uk.icr.io/rxn-test/rxn-reaction-preprocessing \
rxn-tokenize \
/data/example/processed/rxns-small.processed.train.csv \
/data/example/processed/rxns-small.processed.train \
/data/example/processed/rxns-small.processed.validation.csv \
/data/example/processed/rxns-small.processed.validation \
/data/example/processed/rxns-small.processed.test.csv \
/data/example/processed/rxns-small.processed.test
```

## Note about reading CSV files
Pandas appears not to always be able to write a CSV and re-read it if it contains Windows carriage returns.
In order for the scripts to work despite this, all the `pd.read_csv` function calls should include the argument `lineterminator='\n'`.
