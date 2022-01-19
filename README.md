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

## Usage
The following command line scripts are installed with the package.

### rxn-data-pipeline
Wrapper for all other scripts. Allows constructing flexible data pipelines. Entrypoint for Hydra structured configuration.

For an overview of all available configuration parameters and default values, run: `rxn-data-pipeline --cfg job`.

Configuration using YAML (see the file `config.py` for more options and their meaning):
```yaml
defaults:
  - base_config

data:
  path: /tmp/inference/input.csv
  proc_dir: /tmp/rxn-preproc/exp
common:
  sequence:
    # Define which steps and in which order to execute:
    - IMPORT
    - STANDARDIZE
    - PREPROCESS
    - SPLIT
    - TOKENIZE
  fragment_bond: TILDE
preprocess:
  min_products: 0
split:
  split_ratio: 0.05
  max_in_valid: 10000
tokenize:
  input_output_pairs:
    - inp: ${data.proc_dir}/${data.name}.processed.train.csv
      out: ${data.proc_dir}/${data.name}.processed.train
    - inp: ${data.proc_dir}/${data.name}.processed.validation.csv
      out: ${data.proc_dir}/${data.name}.processed.validation
    - inp: ${data.proc_dir}/${data.name}.processed.test.csv
      out: ${data.proc_dir}/${data.name}.processed.test
```
```bash
rxn-data-pipeline --config-dir . --config-name example_config
```

Configuration using command line arguments (example):
```bash
rxn-data-pipeline \
  data.path=/path/to/data/rxns-small.csv \
  data.proc_dir=/path/to/proc/dir \
  common.fragment_bond=TILDE \
  rxn_import.data_format=TXT \
  tokenize.input_output_pairs.0.out=train.txt \
  tokenize.input_output_pairs.1.out=validation.txt \
  tokenize.input_output_pairs.2.out=test.txt
```

## Note about reading CSV files
Pandas appears not to always be able to write a CSV and re-read it if it contains Windows carriage returns.
In order for the scripts to work despite this, all the `pd.read_csv` function calls should include the argument `lineterminator='\n'`.
