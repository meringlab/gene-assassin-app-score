#!/usr/bin/python
import os
import sys
from Queue import *
from threading import *



def make_dir_structure (base_path_with_species, ensembl_relase, dir_to_be_created, execution_version):
    
    
    
    if os.path.exists(base_path_with_species):
        
        new_dir_path_version = os.path.join(base_path_with_species, str(ensembl_relase))  
        new_dir_path = os.path.join(new_dir_path_version, dir_to_be_created)
        
        #new_dir_path_species = os.path.join(new_dir_path_version, species)
        new_dir_path_version_exe_version= os.path.join(new_dir_path, execution_version)
        try:
            if not os.path.exists(new_dir_path_version):
                os.makedirs(new_dir_path_version)
            
            if not os.path.exists(new_dir_path):
                os.makedirs(new_dir_path)
            
            #if not os.path.exists(new_dir_path_species):
                #os.makedirs(new_dir_path_species)
            
            if not os.path.exists(new_dir_path_version_exe_version):
                os.makedirs(new_dir_path_version_exe_version)
        
        except Exception as e:
            exit("\t ... Cannot create intermediary directories %s" % new_dir_path_version_exe_version)
    
    else:
        exit("\t .... Input Base path does not exits %s" % base_path_with_species)
        
        
    return (new_dir_path_version_exe_version)
    
    

def making_output_log_file_names (input_file_path,output_file_path, output_file_descript):
    
    ######### Input_file
    guide_file = input_file_path
    
    
    #########  making of  Output_file and log directory 
    
    gene_name = os.path.basename(guide_file).split("_")[0]
    gene_output_file_name = gene_name + output_file_descript   ####### like "_guide_info.txt"
    
    if os.path.exists(output_file_path):
        
        gene_output_file_name_path = os.path.join(output_file_path, gene_output_file_name)
        #gene_output_file_handle = open (gene_output_file_name_path, "w")        
        
        
        ####### Eception file
        
        log_file_name = os.path.join(output_file_path,"log_file_exceptions%s " %output_file_descript)
        #log_file_handle = open(log_file_name, "a")
        
    else:
        
        exit("Outpur file path does not exist %s" %output_file_path)
        
        
    return (gene_output_file_name_path,log_file_name)
    
    
def search_file_path (name_of_file, search_path_dir):
    file_abs_path = ""
    
    for r,d,f in os.walk(search_path_dir):
        for files in f:
            if files == name_of_file :
                file_abs_path = os.path.join(r,files)
    
    return (file_abs_path)
    
    
def search_dir_path (name_of_dir, dir_path):
    
    folder_abs_path = ""
    
    for r,d,f in os.walk(dir_path):
        for folders in d:
            if folders == name_of_dir :
                folder_abs_path = os.path.join(r,folders)
               
    
    return (folder_abs_path)
  
    
    

    
    
    
    
    
    
    
    
