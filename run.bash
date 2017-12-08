#!/bin/bash

# Step 1. download
## Download the raw files from ENSEMBL ....
python ./src/download/main.py conf/drerio.json


# Step 2. create guides info
python ./src/guides_info/main.py conf/drerio.json


# Step 3. Calculate scores for the guides
## Uses guides info from the Step 2
python ./src/scores/main.py conf/drerio.json

