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
