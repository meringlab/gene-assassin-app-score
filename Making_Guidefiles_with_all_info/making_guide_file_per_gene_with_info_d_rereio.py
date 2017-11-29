#!/usr/bin/python

import os
import timeit

start = timeit.default_timer()
from Universal import universal_function as fn_universal
import making_guide_file_info_function_main as fn_guideinfo_main
stop = timeit.default_timer()
print('time to load modules %dsec' % (stop - start))

######## making a Folder for storing the guides (Output_file)


base_path = "output/v85_2017/danio_rerio/"
dir_to_be_created = "Guide_files_with_information"
species = "danio_rerio"
ensembl_relase = "v_85"
exec_version = "v0"

#####Step 1
guide_file_info_directory = fn_universal.make_dir_structure(base_path, ensembl_relase, dir_to_be_created,exec_version)

### Step 2

# guide_file_path = sys.argv[1]
# for guide_file_path in ['input/ENSDARG00000032959.guides.txt', 'input/ENSDARG00000070986.guides.txt', 'input/ENSDARG00000079029.guides.txt']:
# for guide_file_path in ['input/ENSDARG00000079756.guides.txt', 'input/ENSDARG00000089217.guides.txt']:
# for guide_file_path in ['input/ENSDARG00000089217.guides.txt']:
# for guide_file_path in ['input/ENSDARG00000097979.guides.txt']:

start = timeit.default_timer()
num_processed = 0
for guide_file_path in os.listdir('input'):
    fn_guideinfo_main.making_guide_file_with_info(os.path.join('input', guide_file_path), guide_file_info_directory)
    num_processed += 1
    if num_processed % 100 == 0:
        stop = timeit.default_timer()
        print('%d processed in %dsec' % (num_processed, stop - start))

stop = timeit.default_timer()
print('time to generate guides %dsec' % (stop - start))















