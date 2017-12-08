#!/bin/bash

# download package
## Step 1. Download the raw files from ENSEMBL ....
python ./src/download.py conf/drerio.json

## Step 2. Parse to create transcript cds info file
python ./src/download/parse.py conf/drerio.json

# guide info package

## Step 3 create guides info

python src/info/making_guide_file_per_gene_with_info.py conf/drerio.json


################ Calculating scores for the guides
### ##### Step 7The first argument should be path for the guide_with all info (step 6)
#### Give the base path to new folder where danio_rereio with right version;this will store guide with scores (output from this step)
#### Change the path whereber neccesary

python calculate_scores_guides.py conf/drerio.json

