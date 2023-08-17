# RXN reaction preprocessing

[![Actions tests](https://github.com/rxn4chemistry/rxn-reaction-preprocessing/actions/workflows/tests.yaml/badge.svg)](https://github.com/rxn4chemistry/rxn-reaction-preprocessing/actions)

This repository is devoted to preprocessing chemical reactions: standardization, filtering, etc. 
It also includes code for stable train/test/validation splits and data augmentation.

Links:
* [GitHub repository](https://github.com/rxn4chemistry/rxn-reaction-preprocessing)
* [Documentation](https://rxn4chemistry.github.io/rxn-reaction-preprocessing/)
* [PyPI package](https://pypi.org/project/rxn-reaction-preprocessing/)

## System Requirements

This package is supported on all operating systems.
It has been tested on the following systems:
* macOS: Big Sur (11.1)
* Linux: Ubuntu 18.04.4

A Python version of 3.7 or greater is recommended.

## Installation guide

The package can be installed from Pypi:
```bash
pip install rxn-reaction-preprocessing[rdkit]
```
You can leave out `[rdkit]` if you prefer to install `rdkit` manually (via Conda or Pypi).

For local development, the package can be installed with:
```bash
pip install -e ".[dev]"
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

## Examples

### A pipeline supporting augmentation

A config supporting augmentation of the training split called `train-augmentation-config.yaml`:
```yaml
defaults:
  - base_config

data:
  name: pipeline-with-augmentation
  path: /tmp/file-with-reactions.txt
  proc_dir: /tmp/rxn-preprocessing/experiment
common:
  sequence:
    # Define which steps and in which order to execute:
    - IMPORT
    - STANDARDIZE
    - PREPROCESS
    - SPLIT
    - AUGMENT
    - TOKENIZE
  fragment_bond: TILDE
rxn_import:
  data_format: TXT
preprocess:
  min_products: 1
split:
  input_file_path: ${preprocess.output_file_path}
  split_ratio: 0.05
augment:
  input_file_path: ${data.proc_dir}/${data.name}.processed.train.csv
  output_file_path: ${data.proc_dir}/${data.name}.augmented.train.csv
  permutations: 10
  tokenize: false
  random_type: rotated
tokenize:
  input_output_pairs:
    - inp: ${data.proc_dir}/${data.name}.augmented.train.csv
      out: ${data.proc_dir}/${data.name}.augmented.train
      reaction_column_name: rxn_rotated
    - inp: ${data.proc_dir}/${data.name}.processed.validation.csv
      out: ${data.proc_dir}/${data.name}.processed.validation
    - inp: ${data.proc_dir}/${data.name}.processed.test.csv
      out: ${data.proc_dir}/${data.name}.processed.test
```
```bash
rxn-data-pipeline --config-dir . --config-name train-augmentation-config
```
