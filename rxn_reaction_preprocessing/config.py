#!/usr/bin/env python
# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED
from enum import auto
from enum import Enum
from pathlib import Path
from typing import Any, Optional
from typing import List

from dataclasses import dataclass
from dataclasses import field
from hydra.core.config_store import ConfigStore
from omegaconf import MISSING
from omegaconf import OmegaConf
from omegaconf import SI

from rxn_reaction_preprocessing.utils import RandomType
from rxn_reaction_preprocessing.utils import ReactionSection
from rxn_reaction_preprocessing.utils import standardization_files_directory

OmegaConf.register_new_resolver('stem', lambda p: Path(p).stem)

DEFAULT_ANNOTATION_FILES = [
    str(standardization_files_directory() / 'pistachio-210428.json'),
    str(standardization_files_directory() / 'catalyst-annotation-210428.json'),
    str(standardization_files_directory() / 'catalyst-annotation-210826.json')
]


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


class InitialDataFormat(Enum):
    TXT = auto()
    CSV = auto()


class Step(Enum):
    IMPORT = auto()
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
        reaction_column_name: Name of the reaction column for the data file.
    """
    sequence: List[Step] = field(
        default_factory=lambda: [
            Step.IMPORT,
            Step.STANDARDIZE,
            Step.PREPROCESS,
            Step.SPLIT,
            Step.TOKENIZE,
        ]
    )
    fragment_bond: FragmentBond = FragmentBond.DOT
    reaction_column_name: str = 'rxn'


@dataclass
class RxnImportConfig:
    """Configuration for the initial import of the reaction data.

    Fields:
        input_file: the input file path (.txt, .csv).
        output_csv: the ouptut file path.
        initial_data_format: whether the input file is in TXT or CSV format.
        reaction_column_name: name the column containing the reactions if the input
            is in CSV format. The value is ignored if the input is not in CSV format.
        reaction_column_name: how to name the column with the reaction SMILES in
            the output CSV.
        fragment_bond: token denoting a fragment bond in the reaction SMILES (after import).
        remove_atom_mapping: whether to remove the atom mapping from the input reaction SMILES.
        column_for_light: name of the column containing boolean values that indicate
            whether the reaction happens under light. If specified, the light token will
            be added to the precursors of the corresponding reactions.
        column_for_heat: name of the column containing boolean values that indicate
            whether the reaction happens under heat. If specified, the heat token will
            be added to the precursors of the corresponding reactions.
    """
    input_file: str = SI('${data.path}')
    output_csv: str = SI('${data.proc_dir}/${data.name}.imported.csv')
    data_format: InitialDataFormat = InitialDataFormat.CSV
    input_csv_column_name: str = SI('${common.reaction_column_name}')
    reaction_column_name: str = SI('${common.reaction_column_name}')
    fragment_bond: FragmentBond = SI('${common.fragment_bond}')
    remove_atom_mapping: bool = True
    column_for_light: Optional[str] = None
    column_for_heat: Optional[str] = None


@dataclass
class StandardizeConfig:
    """Configuration for the standardization transformation step.

    Fields:
        input_file_path: The input CSV (one SMILES per line).
        output_file_path: The output file path containing the result after standardization.
        annotation_file_paths: The files to load the annotated molecules from.
        discard_unannotated_metals: whether reactions containing unannotated
            molecules with transition metals must be rejected.
        fragment_bond: Token used to denote a fragment bond in the reaction SMILES.
        reaction_column_name: Name of the reaction column for the data file.
        remove_stereo_if_not_defined_in_precursors: Remove chiral centers from product.
    """
    input_file_path: str = SI('${rxn_import.output_csv}')
    annotation_file_paths: List[str] = field(default_factory=lambda: DEFAULT_ANNOTATION_FILES)
    discard_unannotated_metals: bool = True
    output_file_path: str = SI('${data.proc_dir}/${data.name}.standardized.csv')
    fragment_bond: FragmentBond = SI('${common.fragment_bond}')
    reaction_column_name: str = SI('${common.reaction_column_name}')
    remove_stereo_if_not_defined_in_precursors: bool = False


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
        fragment_bond: Token used to denote a fragment bond in the reaction SMILES.
        reaction_column_name: Name of the reaction column for the data file.
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
    reaction_column_name: str = SI('${common.reaction_column_name}')


@dataclass
class AugmentConfig:
    """Configuration for the augmentation transformation step.

    Fields:
        input_file_path:  The input file path (one SMILES per line).
        output_file_path: The output file path.
        tokenize: if tokenization is to be performed
        random_type: The randomization type to be applied
        permutations: number of randomic permutations for input SMILES
        reaction_column_name: Name of the reaction column for the data file.
        rxn_section_to_augment: The section of the rxn SMILES to augment.
            "precursors" for augmenting only the precursors
            "products" for augmenting only the products
        fragment_bond: Token used to denote a fragment bond in the reaction SMILES.
    """
    input_file_path: str = SI('${preprocess.output_file_path}')
    output_file_path: str = SI('${data.proc_dir}/${data.name}.augmented.csv')
    tokenize: bool = True
    random_type: RandomType = RandomType.unrestricted
    permutations: int = 1
    reaction_column_name: str = SI('${common.reaction_column_name}')
    rxn_section_to_augment: ReactionSection = ReactionSection.precursors
    fragment_bond: FragmentBond = SI('${common.fragment_bond}')


@dataclass
class SplitConfig:
    """Configuration for the split transformation step.

    Fields:
        input_file_path: The input file path.
        output_directory: The directory containing the files after splitting.
        split_ratio: The split ratio between training, and test and validation sets.
        reaction_column_name: Name of the reaction column for the data file.
        index_column: Name of the column used to generate the hash ensuring
            stable splitting.
        max_in_valid: Maximal number of reactions to keep in the validation set.
            This can be useful to avoid unnecessarily prolonging training. It
            is considered as an approximate limit (due to randomness in hashing).
            Defaults to no restriction.
        seed: Seed for the hashing function used for splitting.
    """
    input_file_path: str = SI('${preprocess.output_file_path}')
    output_directory: str = SI('${data.proc_dir}')
    split_ratio: float = 0.05
    reaction_column_name: str = SI('${common.reaction_column_name}')
    index_column: str = SI('${split.reaction_column_name}')
    max_in_valid: Optional[int] = None
    seed: int = 42


@dataclass
class InputOutputTuple:
    inp: str = MISSING
    out: str = MISSING


@dataclass
class TokenizeConfig:
    """Configuration for the tokenization transformation step.

    Fields:
        input_output_pairs: Paths to the input and output files.
        reaction_column_name: Name of the reaction column for the data file.
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
    reaction_column_name: str = SI('${common.reaction_column_name}')


@dataclass
class Config:
    data: DataConfig = field(default_factory=DataConfig)
    common: CommonConfig = field(default_factory=CommonConfig)
    rxn_import: RxnImportConfig = field(default_factory=RxnImportConfig)
    standardize: StandardizeConfig = field(default_factory=StandardizeConfig)
    preprocess: PreprocessConfig = field(default_factory=PreprocessConfig)
    augment: AugmentConfig = field(default_factory=AugmentConfig)
    split: SplitConfig = field(default_factory=SplitConfig)
    tokenize: TokenizeConfig = field(default_factory=TokenizeConfig)

    @classmethod
    def from_generic_config(cls, config: Any) -> 'Config':
        cfg_dict = OmegaConf.merge(OmegaConf.structured(Config), config)
        cfg = OmegaConf.to_object(cfg_dict)
        return cfg  # type: ignore


cs = ConfigStore.instance()
cs.store(group='data', name='base_data', node=DataConfig)
cs.store(group='common', name='base_common', node=CommonConfig)
cs.store(group='rxn_import', name='base_rxn_import', node=RxnImportConfig)
cs.store(group='standardize', name='base_standardize', node=StandardizeConfig)
cs.store(group='preprocess', name='base_preprocess', node=PreprocessConfig)
cs.store(group='augment', name='base_augment', node=AugmentConfig)
cs.store(group='tokenize', name='base_tokenize', node=TokenizeConfig)
cs.store(group='split', name='base_split', node=SplitConfig)
cs.store(name='base_config', node=Config)
