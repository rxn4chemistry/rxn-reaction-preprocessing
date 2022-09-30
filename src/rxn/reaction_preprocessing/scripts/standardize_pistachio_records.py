#!/usr/bin/env python
# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED
import json
import logging
from pathlib import Path

import hydra
from rxn.utilities.files import iterate_lines_from_file

from rxn.reaction_preprocessing.config import Config, Step
from rxn.reaction_preprocessing.pistachio_record_standardizer import (
    PistachioRecordStandardizer,
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@hydra.main(config_name="base_config", config_path=None)
def main(cfg: Config) -> None:
    """Standardize reactions in a Pistachio JSON.

    This will do steps similar to the STANDARDIZE and PREPROCESS steps from the
    normal (CSV-based) reaction preprocessing.

    Relies on the same config format as the main pipeline, although not all the
    values from there are needed.
    """

    # Enforce config schema. Will also convert strings to Enums when necessary.
    cfg = Config.from_generic_config(cfg)

    # make sure that the required output directories exist
    Path(cfg.data.proc_dir).mkdir(parents=True, exist_ok=True)

    if cfg.common.sequence != [Step.STANDARDIZE, Step.PREPROCESS]:
        raise SystemExit(
            "Standardization of Pistachio: steps must be STANDARDIZE and PREPROCESS exactly"
        )

    pistachio_standardizer = PistachioRecordStandardizer(
        cfg_standardize=cfg.standardize, cfg_preprocess=cfg.preprocess
    )

    jsonl_file = cfg.data.path
    output_jsonl = Path(cfg.data.proc_dir) / "processed.jsonl"

    with open(output_jsonl, "wt") as f:
        for json_line in iterate_lines_from_file(jsonl_file):
            reaction_record = json.loads(json_line)

            try:
                updated_record = pistachio_standardizer.standardize(reaction_record)
            except Exception as e:
                logger.info(f"Ignoring record: {e}")
                continue

            f.write(json.dumps(updated_record))
            f.write("\n")


if __name__ == "__main__":
    main()
