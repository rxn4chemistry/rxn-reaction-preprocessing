#!/usr/bin/env python
# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED
import logging
from pathlib import Path

import hydra
from omegaconf import OmegaConf

from rxn.reaction_preprocessing import __version__
from rxn.reaction_preprocessing.augmenter import augment
from rxn.reaction_preprocessing.config import Config, Step
from rxn.reaction_preprocessing.importer import rxn_import
from rxn.reaction_preprocessing.preprocessor import preprocess
from rxn.reaction_preprocessing.smiles_tokenizer import tokenize
from rxn.reaction_preprocessing.stable_data_splitter import split
from rxn.reaction_preprocessing.standardizer import standardize
from rxn.reaction_preprocessing.utils import (
    add_custom_logger_to_file,
    overwrite_logging_format,
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def preprocess_data(cfg: Config) -> None:
    """Preprocess data to generate a dataset for training transformer models."""

    # Enforce config schema. Will also convert strings to Enums when necessary.
    cfg = Config.from_generic_config(cfg)

    logger.info(
        f"Preprocessing reaction data with rxn-reaction-preprocessing, "
        f"version {__version__}."
    )
    logger.info(
        "Running with the following configuration:\n"
        f"{OmegaConf.to_yaml(cfg, resolve=True)}"
    )

    # Create the output directory (may exist already if this function was
    # called from the main script).
    processing_dir = Path(cfg.data.proc_dir)
    processing_dir.mkdir(parents=True, exist_ok=True)

    # Save the config
    with open(processing_dir / "preprocessing_config.yml", "wt") as f:
        f.write(OmegaConf.to_yaml(cfg, resolve=True))

    for step in cfg.common.sequence:
        logger.info(f"Running step: {step.name}")

        if step is Step.IMPORT:
            rxn_import(cfg.rxn_import)
        elif step is Step.STANDARDIZE:
            standardize(cfg.standardize)
        elif step is Step.PREPROCESS:
            preprocess(cfg.preprocess)
        elif step is Step.AUGMENT:
            augment(cfg.augment)
        elif step is Step.SPLIT:
            split(cfg.split)
        elif step is Step.TOKENIZE:
            tokenize(cfg.tokenize)


@hydra.main(config_name="base_config", config_path=None)
def data_pipeline(cfg: Config) -> None:
    """Preprocess data to generate a dataset for training transformer models."""

    # Enforce config schema. Will also convert strings to Enums when necessary.
    cfg = Config.from_generic_config(cfg)

    # Setup logging to file, and overwrite the log format
    processing_dir = Path(cfg.data.proc_dir)
    processing_dir.mkdir(parents=True, exist_ok=True)
    add_custom_logger_to_file(processing_dir / "log.txt")
    overwrite_logging_format()

    preprocess_data(cfg)


if __name__ == "__main__":
    data_pipeline()
