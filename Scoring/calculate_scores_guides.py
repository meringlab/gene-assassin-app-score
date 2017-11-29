import os

import timeit

start = timeit.default_timer()
from Universal import universal_function as fn_universal

import all_scoring_function as fn_scoring
import scoring_function_main as fn_scoring_main

################## Function calling

base_path = "output/v85_2017/danio_rerio/"
dir_to_be_created = "Guide_files_with_scores"
species = "danio_rerio"
ensembl_relase = "v_85"
exec_version = "v0"

#####Step 1
guide_file_score_directory = fn_universal.make_dir_structure(base_path, ensembl_relase,dir_to_be_created,exec_version)

output_file_path =  guide_file_score_directory
output_file_descript = "_scores.txt"

### Step 2
#### Transcript cds info dict

transcript_cds_file_name = "transcript_cds_info_parsed_Danio_rerio.GRCz10.85.gtf.txt"
transcript_cds_dir_path = "output/v85_2017/danio_rerio/v_85"
transcript_cds_info_dict = fn_scoring.make_transcript_cds_info_dict (transcript_cds_file_name, transcript_cds_dir_path)


###### Protein dict
name_protein_dir = "proteins"
protein_dir_path = "input/"
protein_domain_info_dict =  fn_scoring.make_protein_dict (name_protein_dir, protein_dir_path)


###### Variation dict
variation_file_name = "Danio_rerio.gvf.gz"
var_dir_path = "output/v85_2017/danio_rerio/v_85"
snv_dict = fn_scoring.make_var_dict(variation_file_name,var_dir_path)
stop = timeit.default_timer()
print('time to prepare for computation %dsec' % (stop - start))

start = timeit.default_timer()


################ callling the main function
# input_file_path = sys.argv[1]
num_processed = 0

guides_info_dir = 'output/v85_2017/danio_rerio/v_85/Guide_files_with_information/v0/'
for input_file in os.listdir(guides_info_dir):
    input_file_path = os.path.join(guides_info_dir, input_file)
    # input_file_path = os.path.join(guides_info_dir, 'ENSDARG00000079029_guides_info.txt')
    fn_scoring_main.calculate_score_guide_main (input_file_path, output_file_path, output_file_descript, transcript_cds_info_dict, protein_domain_info_dict,snv_dict)
    num_processed += 1
    if num_processed % 10 == 0:
        stop = timeit.default_timer()
        print('%d processed in %dsec' % (num_processed, stop - start))

stop = timeit.default_timer()
print('time to compute scores %dsec' % (stop - start))

                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
            
            
            
            
            

