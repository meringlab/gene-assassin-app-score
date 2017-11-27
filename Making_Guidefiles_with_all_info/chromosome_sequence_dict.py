#!/usr/bin/env python
import os
import sys
import gzip


##### Give the path for Danio_rerio.GRCz10.dna.toplevel.fa.gz

chromosome_data_path = "/home/neha/Projects/Christian_Mosimann/Crispr_Project/v85_2017/danio_rerio/v_85/Raw_data_files/Danio_rerio.GRCz10.dna.toplevel.fa.gz"



################# Processing

try :
    os.path.exists(chromosome_data_path)
except Exception as e:
    print e
    exit("\t ... chromosome_data_path does not exist %s."% chromosome_data_path)
    




input_file_handle = gzip.open(chromosome_data_path,"rb")

seq= ""
chromosome_sequence_dict = {}

for line in input_file_handle:
    
    if line.startswith(">"):
        
        if seq!= "":
            
            #print l ,chromosome_name, len(seq)  ####### Good check to see if everything is rightly stored
            
            if chromosome_name not in chromosome_sequence_dict:
                chromosome_sequence_dict[chromosome_name] = seq
                
            else:
                print "found"
                
        seq = ""
        
        l = line.strip("\n").strip(">").split(" ")
        chromosome_name = l[0]
            
    else:
        seq = seq + line.strip("\n")




chromosome_sequence_dict[chromosome_name] = seq

