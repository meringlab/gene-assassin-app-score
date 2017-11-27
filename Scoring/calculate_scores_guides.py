import os, sys
import ast
import math

sys.path.append("/home/neha/Projects/Christian_Mosimann/Crispr_Project/v85_2017/Crispr_dummy_1/")
from Script.Universal import universal_function as fn_universal

import all_scoring_function as fn_scoring
import scoring_function_main as fn_scoring_main


################## Function calling

base_path = "/home/neha/Projects/Christian_Mosimann/Crispr_Project/v85_2017/danio_rerio/"
dir_to_be_created = "Guide_files_with_scores"
species = "danio_rerio"
ensembl_relase = "v_85"
exec_version = "v0"

#####Step 1
guide_file_score_directory = fn_universal.make_dir_structure(base_path, ensembl_relase,dir_to_be_created,exec_version)

input_file_path = sys.argv[1]
output_file_path =  guide_file_score_directory
output_file_descript = "_scores.txt"

### Step 2
#### Transcript cds info dict

transcript_cds_file_name = "transcript_cds_info_parsed_Danio_rerio.GRCz10.85.gtf.txt"
transcript_cds_dir_path = "/home/neha/Projects/Christian_Mosimann/Crispr_Project/v85_2017/danio_rerio/v_85"
transcript_cds_info_dict = fn_scoring.make_transcript_cds_info_dict (transcript_cds_file_name, transcript_cds_dir_path)


###### Protein dict
name_protein_dir = "proteins"
protein_dir_path = "/home/neha/Projects/Christian_Mosimann/Crispr_Project/v85_2017/zebrafish_v85_new/Data_files/v_85/"
protein_domain_info_dict =  fn_scoring.make_protein_dict (name_protein_dir, protein_dir_path)


###### Variation dict
variation_file_name = "Danio_rerio.gvf.gz"
var_dir_path = "/home/neha/Projects/Christian_Mosimann/Crispr_Project/v85_2017/danio_rerio/v_85"
snv_dict = fn_scoring.make_var_dict(variation_file_name,var_dir_path)

################ callling the main function

fn_scoring_main.calculate_score_guide_main (input_file_path, output_file_path, output_file_descript,transcript_cds_info_dict, protein_domain_info_dict,snv_dict)




                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
            
            
            
            
            

