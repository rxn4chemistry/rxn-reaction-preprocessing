#!/usr/bin/env python
# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED
from dataclasses import dataclass
from dataclasses import field
from enum import auto
from enum import Enum
from pathlib import Path
from typing import Any
from typing import List

from hydra.core.config_store import ConfigStore
from omegaconf import MISSING
from omegaconf import OmegaConf
from omegaconf import SI

from rxn_reaction_preprocessing.utils import RandomType

OmegaConf.register_new_resolver('stem', lambda p: Path(p).stem)


@dataclass
class DataConfig:
    """Configuration of data sources and intermediate storage.

    Fields:
        path: Absolute path to input data file.
        name: Name of the input data file (without extension).
        proc_dir: Directory for storing intermediate and final output files.
    """
    path: str = MISSING
    name: str = SI('${stem:${data.path}}')
    proc_dir: str = MISSING


class FragmentBond(Enum):
    DOT = '.'
    TILDE = '~'


class Step(Enum):
    STANDARDIZE = auto()
    PREPROCESS = auto()
    AUGMENT = auto()
    SPLIT = auto()
    TOKENIZE = auto()


@dataclass
class CommonConfig:
    """Configuration used by multiple steps.

    Fields:
        sequence: Ordered sequence of data transformation steps to perform.
        fragment_bond: Token used to denote a fragment bond in the SMILES of the reactions to process.
    """
    sequence: List[Step] = field(
        default_factory=lambda: [Step.STANDARDIZE, Step.PREPROCESS, Step.SPLIT, Step.TOKENIZE]
    )
    fragment_bond: FragmentBond = FragmentBond.DOT


@dataclass
class StandardizeConfig:
    """Configuration for the standardization transformation step.

    Fields:
        input_file_path: The input file path (one SMILES per line).
        output_file_path: The output file path containing the result after standardization.
    """
    input_file_path: str = SI('${data.path}')
    output_file_path: str = SI('${data.proc_dir}/${data.name}.standardized.csv')
    fragment_bond: FragmentBond = SI('${common.fragment_bond}')


@dataclass
class PreprocessConfig:
    """Configuration for the preprocess transformation step.

    Fields:
        input_file_path: The input file path (one reaction SMARTS per line).
        output_file_path: The output file path containing the result after preprocessing.
        min_reactants: The minimum number of reactants.
        max_reactants: The maximum number of reactants.
        max_reactants_tokens: The maximum number of reactants tokens.
        min_agents: The minimum number of agents.
        max_agents: The maximum number of agents.
        max_agents_tokens: The maximum number of agents tokens.
        min_products: The minimum number of products.
        max_products: The maximum number of products.
        max_products_tokens: The maximum number of products tokens.
        max_absolute_formal_charge: The maximum absolute formal charge.
    """
    input_file_path: str = SI('${standardize.output_file_path}')
    output_file_path: str = SI('${data.proc_dir}/${data.name}.processed.csv')
    min_reactants: int = 2
    max_reactants: int = 10
    max_reactants_tokens: int = 300
    min_agents: int = 0
    max_agents: int = 0
    max_agents_tokens: int = 0
    min_products: int = 1
    max_products: int = 1
    max_products_tokens: int = 200
    max_absolute_formal_charge: int = 2
    fragment_bond: FragmentBond = SI('${common.fragment_bond}')


@dataclass
class AugmentConfig:
    """Configuration for the augmentation transformation step.

    Fields:
        input_file_path:  The input file path (one SMILES per line).
        output_file_path: The output file path.
        tokenize: if tokenization is to be performed
        random_type: The randomization type to be applied
        permutations: number of randomic permutations for input SMILES
    """
    input_file_path: str = SI('${preprocess.output_file_path}')
    output_file_path: str = SI('${data.proc_dir}/${data.name}.augmented.csv')
    tokenize: bool = True
    random_type: RandomType = RandomType.unrestricted
    permutations: int = 1
    fragment_bond: FragmentBond = SI('${common.fragment_bond}')


@dataclass
class SplitConfig:
    """Configuration for the split transformation step.

    Fields:
        input_file_path: The input file path.
        output_directory: The directory containing the files after splitting.
        split_ratio: The split ratio between training, and test and validation sets.
        seed: Seed for the hashing function used for splitting.
    """
    input_file_path: str = SI('${preprocess.output_file_path}')
    output_directory: str = SI('${data.proc_dir}')
    split_ratio: float = 0.05
    index_column: str = 'rxn'
    seed: int = 42


@dataclass
class InputOutputTuple:
    inp: str = MISSING
    out: str = MISSING


@dataclass
class TokenizeConfig:
    """Configuration for the tokenization tranformation step.

    Fields:
        input_output_pairs: Paths to the input and output files.
    """
    input_output_pairs: List[InputOutputTuple] = field(
        default_factory=lambda: [
            InputOutputTuple(
                SI('${data.proc_dir}/${data.name}.processed.train.csv'),
                SI('${data.proc_dir}/${data.name}.processed.train')
            ),
            InputOutputTuple(
                SI('${data.proc_dir}/${data.name}.processed.validation.csv'),
                SI('${data.proc_dir}/${data.name}.processed.validation')
            ),
            InputOutputTuple(
                SI('${data.proc_dir}/${data.name}.processed.test.csv'),
                SI('${data.proc_dir}/${data.name}.processed.test')
            ),
        ]
    )


@dataclass
class Config:
    defaults: List[Any] = field(
        default_factory=lambda: [
            {
                'data': 'base_data'
            }, {
                'common': 'base_common'
            }, {
                'standardize': 'base_standardize'
            }, {
                'preprocess': 'base_preprocess'
            }, {
                'augment': 'base_augment'
            }, {
                'split': 'base_split'
            }, {
                'tokenize': 'base_tokenize'
            }
        ]
    )

    data: DataConfig = MISSING
    common: CommonConfig = MISSING
    standardize: StandardizeConfig = MISSING
    preprocess: PreprocessConfig = MISSING
    augment: AugmentConfig = MISSING
    split: SplitConfig = MISSING
    tokenize: TokenizeConfig = MISSING


cs = ConfigStore.instance()
cs.store(group='data', name='base_data', node=DataConfig)
cs.store(group='common', name='base_common', node=CommonConfig)
cs.store(group='standardize', name='base_standardize', node=StandardizeConfig)
cs.store(group='preprocess', name='base_preprocess', node=PreprocessConfig)
cs.store(group='augment', name='base_augment', node=AugmentConfig)
cs.store(group='tokenize', name='base_tokenize', node=TokenizeConfig)
cs.store(group='split', name='base_split', node=SplitConfig)
cs.store(name='base_config', node=Config)
