#!/usr/bin/env python
# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED
from pathlib import Path

import hydra
from omegaconf import OmegaConf

from rxn_reaction_preprocessing.augmenter import augment
from rxn_reaction_preprocessing.config import Config
from rxn_reaction_preprocessing.config import Step
from rxn_reaction_preprocessing.preprocessor import preprocess
from rxn_reaction_preprocessing.smiles_tokenizer import tokenize
from rxn_reaction_preprocessing.stable_data_splitter import split
from rxn_reaction_preprocessing.standardizer import standardize


@hydra.main(config_name='base_config', config_path=None)
def data_pipeline(cfg: Config) -> None:
    """Preprocess data to generate a dataset for training transformer models."""
    print(f'Running with the following configuration:\n\n{OmegaConf.to_yaml(cfg, resolve=True)}\n')

    # Enforce config schema. Will also convert strings to Enums when necessary.
    cfg = Config.from_generic_config(cfg)

    # make sure that the required output directories exist
    Path(cfg.data.proc_dir).mkdir(parents=True, exist_ok=True)

    for step in cfg.common.sequence:
        print(f'Running step: {step.name}')

        if step is Step.STANDARDIZE:
            standardize(cfg.standardize)
        elif step is Step.PREPROCESS:
            preprocess(cfg.preprocess)
        elif step is Step.AUGMENT:
            augment(cfg.augment)
        elif step is Step.SPLIT:
            split(cfg.split)
        elif step is Step.TOKENIZE:
            tokenize(cfg.tokenize)


if __name__ == '__main__':
    data_pipeline()
