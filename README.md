# data_preprocessor

Repository to apply many preprocessing transformations including stable train/test/validation splits.

## Installation


```python
conda create -n data_preprocessor python=3.6
conda activate data_preprocessor
conda install rdkit -c rdkit
conda install sentencepiece==0.1.83
conda install -c powerai pytorch=1.3.1
pip install -e .

