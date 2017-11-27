#!/bin/bash

################ Go to Script/Downloading and creating Data files
##### Step 1. Execute the script to download the raw files from ENSEMBL .... You need to change temfile path to whicever path you wouldlike
python ./Download_create_data_files/download_raw_datafiles.py
#### Step 2. Make parsed GTF file  ...... Give the path for gtf file
python ./Download_create_data_files/parsing_gtf_files.py  
##### Step 3 Using the parsed file to create transcript cds info file  ..... give the path for parsed gtf file form step 2
python ./Download_create_data_files/making_transcript_cds_info_file.py


################ Go to Script/Making Guide file with all info

########### Make all the necessary dict
### # Step 4 Maake a chromosome dict  ####### ##### Give the path for Danio_rerio.GRCz10.dna.toplevel.fa.gz
python chromosome_sequence_dict.py 

#### Step 5 Make exon dict
python exon_dict.py  ######### ################ Give the path for parsed_Danio_rerio.GRCz10.85.gtf.txt

#### Step 6 Run the script for all the info for guides
##### The first argument should be path for the guide 
#### Give the base path to new folder where danio_rereio with right version;this will store guide with info (output from this step)
###### also give the path of where the script is stored sys.path.append

python making_guide_file_per_gene_with_info_d_rereio.py /Users/neha/Desktop/Crispr_dummy/Guide_files/ENSDARG00000079029.guides.txt


################ Calculating scores for the guides
### ##### Step 7The first argument should be path for the guide_with all info (step 6)
#### Give the base path to new folder where danio_rereio with right version;this will store guide with scores (output from this step)
#### Change the path whereber neccesary

python calculate_scores_guides.py /home/neha/Projects/Christian_Mosimann/Crispr_Project/v85_2017/danio_rerio/v_85/Guide_files_with_information/v0/ENSDARG00000079029_guides_info.txt

