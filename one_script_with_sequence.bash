#!/bin/bash

################ Go to Script/Downloading and creating Data files
##### Step 1. Execute the script to download the raw files from ENSEMBL .... You need to change temfile path to whicever path you wouldlike
python ./Download_create_data_files/download_raw_datafiles.py conf/drerio.json

#### Step 2. Make parsed GTF file  ...... Give the path for gtf file
python ./Download_create_data_files/parsing_gtf_files.py conf/drerio.json

##### Step 3 Using the parsed file to create transcript cds info file  ..... give the path for parsed gtf file form step 2
python ./Download_create_data_files/making_transcript_cds_info_file.py conf/drerio.json


################ Go to Script/Making Guide file with all info

########### Make all the necessary dict
#### # Step 4 Make a chromosome dict  ####### ##### Give the path for Danio_rerio.GRCz10.dna.toplevel.fa.gz
#python chromosome_sequence_dict.py
#
##### Step 5 Make exon dict
#python exon_dict.py  ######### ################ Give the path for parsed_Danio_rerio.GRCz10.85.gtf.txt

#### Step 6 Run the script for all the info for guides
##### The first argument should be path for the guide 
#### Give the base path to new folder where danio_rereio with right version;this will store guide with info (output from this step)

python making_guide_file_per_gene_with_info.py conf/drerio.json


################ Calculating scores for the guides
### ##### Step 7The first argument should be path for the guide_with all info (step 6)
#### Give the base path to new folder where danio_rereio with right version;this will store guide with scores (output from this step)
#### Change the path whereber neccesary

python calculate_scores_guides.py conf/drerio.json

