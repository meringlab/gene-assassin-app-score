#!/usr/bin/python
import os
import sys
import ast



sys.path.append("/home/neha/Projects/Christian_Mosimann/Crispr_Project/v85_2017/Crispr_dummy_1/")
from Script.Universal import universal_function as fn_universal

import making_guide_file_info_function_main as fn_guideinfo_main

######## making a Folder for storing the guides (Output_file)


base_path = "/home/neha/Projects/Christian_Mosimann/Crispr_Project/v85_2017/danio_rerio/"
dir_to_be_created = "Guide_files_with_information"
species = "danio_rerio"
ensembl_relase = "v_85"
exec_version = "v0"

#####Step 1
guide_file_info_directory = fn_universal.make_dir_structure(base_path, ensembl_relase, dir_to_be_created,exec_version)

### Step 2

guide_file_path = sys.argv[1]
fn_guideinfo_main.making_guide_file_with_info(guide_file_path,guide_file_info_directory)
















