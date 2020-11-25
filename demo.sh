#!/bin/sh

docker run -v /mnt/c/Users/DanielProbst/Code/data_preprocessor/data/example/raw/rxns-small.txt:/data/input.txt -v /mnt/c/Users/DanielProbst/Code/data_preprocessor/data/example/processed/:/data/output demo/rxn-reaction-preprocessor
docker run -v /mnt/c/Users/DanielProbst/Code/data_preprocessor/data/example/processed/processed.csv:/data/input.txt -v /mnt/c/Users/DanielProbst/Code/data_preprocessor/data/example/processed/:/data/output demo/rxn-reaction-splitter
docker run -v /mnt/c/Users/DanielProbst/Code/data_preprocessor/data/example/processed:/data/input demo/rxn-reaction-tokenizer