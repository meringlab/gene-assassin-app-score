#!/bin/bash

source activate ga-scores
export PYTHONPATH=.:src

for genome in `ls conf`
do

    # Step 1. download
    ## Download the raw files from ENSEMBL ....
    python ./src/download/main.py conf/$genome


    # Step 2. create guides info
    python ./src/guides_info/main.py conf/$genome


    # Step 3. Calculate scores for the guides
    ## Uses guides info from the Step 2
    python ./src/scores/main.py conf/$genome

done