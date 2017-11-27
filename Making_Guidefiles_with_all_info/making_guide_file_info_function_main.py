#!/usr/bin/python
import os
import sys
import ast


import making_guide_files_per_gene_functions as fn_guideinfo


########## Processing the file funcion

def making_guide_file_with_info (guide_file_path,guide_file_info_directory):
    
    ######### Input_file
    guide_file = guide_file_path
    
    
    #########  making of  Output_file and log directory 
    
    gene_name = os.path.basename(guide_file).split(".")[0]
    gene_output_file_name = gene_name + "_guides_info.txt"
    
    if os.path.exists(guide_file_info_directory):
        
        output_file_path = guide_file_info_directory
        gene_output_file_name_path = os.path.join(output_file_path, gene_output_file_name)
        gene_output_file_handle = open (gene_output_file_name_path, "w")        
        
        
        ####### Eception file
        
        log_file_name = os.path.join(guide_file_info_directory,"log_file_exceptions.txt")
        log_file_handle = open(log_file_name, "a")
        
    else:
        
        exit("Guide file info dir doesnot exist %s" %guide_file_info_directory)
        
        
        
    ######### Header
    
    header = ["Gene_id", "Guide_seq_with_ngg", "Guide_seq", "Guide_chr_no", "Guide_start","Guide_stop", "Guide_strand", "Guide_uniqness", \
              "Guide_exon_list", "Exon_biotype_list", \
            "cutsite18", "dist_cds_start_cutsite", "dist_cds_stop_cutsite",\
              "Total_transcript_count", "Transcript_list", "Total_exons_in_transcript_list", "Exon_rank_in_transcript_list",\
               "microhomology_seq"]
    
    header_output = "\t".join(header) + "\n"

    ###### Writing header
    gene_output_file_handle.write(header_output)
    
    
    ################## Processing the file

    guide_input_file = open(guide_file)
    
    for line in guide_input_file:
        
        l = line.strip("\n").split("\t")
        guide_chr = l [0]; guide_start = l[1]; guide_stop = l[2]; guide_strand = l[3]; guide_seq = l[4]; guide_offtarget_profile = l[10]
        
        
    
    ###### Calculate cutsite18
    
        try:
            guide_cutsite18 = fn_guideinfo.calculate_cutsite18_guide(guide_strand,guide_start,guide_stop)  
        except Exception as e:
            log_output = gene_name+ "\t" + str(e) + "\t" + "guide_cutsite18" + "\n"
            log_file_handle.write(log_output)
            
    
    
    ##### Searching for exon_id
        
        guide_exon_list = []
        for exon_location in fn_guideinfo.exon_dict.location_exon_id:
            exon_chr = exon_location.split(":")[0]
            exon_start_stop = exon_location.split(":")[1]
            
            if guide_chr == exon_chr: ##### 1st step
                exon_start = float(exon_start_stop.split("-")[0])
                exon_stop = float(exon_start_stop.split("-")[1])
                ## Guide should be between exon start/stop
                if exon_start <= exon_stop:
                    if float(guide_cutsite18) >= exon_start and float(guide_cutsite18) <= exon_stop:
                        guide_exon = fn_guideinfo.exon_dict.location_exon_id[exon_location]
                        guide_exon_list = guide_exon_list + guide_exon
                    else: ###### 
                        pass        
                else:
                    e = "Exon_start > Exon_stop"
                    log_output = gene_name + "\t" + str(e) + "\t" + "Exon_id_search_problem" + "\n"
                    log_file_handle.write(log_output)
            else:
                pass
    
        ####### If exon list is empty because there is no protein-coding transcript then the gene_id_list would also be empty
        ######## there are two cases here: one when exon list is empty becaus ethere is no protein coding transcript and also when exon list is not empty there ,
        #### after overlapping removal it still would contain gene from input file
        
        ##### there are three cases: 1. Non-protein coding transcript (guide exon list  = 0 )2. Overlapping genes 3. Both  (guide_exon_list_updated = 0)
        
        
        if len(guide_exon_list) !=0 :                         #### It could belong to case 2 and case 3
            gene_id_list = [fn_guideinfo.exon_dict.exon_gene_info[x]["gene_id"] for x in guide_exon_list]  ###### To detect problem in gene_id
            ######### To remove overlapping genes
            guide_exon_list_remove = []
            
            for exon in guide_exon_list:
                gene_id = fn_guideinfo.exon_dict.exon_gene_info[exon]["gene_id"]
                if gene_id != gene_name:
                    guide_exon_list_remove.append(exon)
                    
            guide_exon_list_updated = [x for x in guide_exon_list if x not in guide_exon_list_remove]
            gene_id_list_updated = [fn_guideinfo.exon_dict.exon_gene_info[x]["gene_id"] for x in guide_exon_list_updated]
            unique_gene_id = list(set(gene_id_list_updated))
            
            if len(gene_id_list_updated) != 0 :
                if len(unique_gene_id) !=1 and unique_gene_id[0] != gene_name:
                    e = "Gene_id not match input file after overlapping"
                    log_output = gene_name+ "\t" + str(e) + "\t" + "gene_id_match_after_overlap" 
                    log_file_handle.write(log_output)
        else:
            guide_exon_list_updated = guide_exon_list
            
            
    ######### Getting guide info
    
        try:
        ####### Guide_uniqueness
                    
            offtarget_profile_list = guide_offtarget_profile.split(",")
            guide_uniq = offtarget_profile_list[1]
                    
            ###### seq_with_ngg
            guide_seq_with_ngg = fn_guideinfo.get_guide_seq_with_ngg(guide_start, guide_stop, guide_strand, guide_chr)
            
             ########## Extracting feature
            cds_start_stop_cut_site = fn_guideinfo.calculate_cutsite18_dist_fom_exon_cds_start_stop_for_exon_list_modified(guide_exon_list_updated,guide_cutsite18)
            exon_genomic_features = fn_guideinfo.extract_exon_features_from_gtf_for_exonlist_modified(guide_exon_list_updated)
            
            if len(guide_exon_list_updated) != 0:
                microhomology_sequence = fn_guideinfo.generate_seq_for_microhomology_scoring(guide_start, guide_stop, guide_strand, guide_chr)
                exon_biotype_list = [fn_guideinfo.exon_dict.exon_id_biotype[x]["biotype"] for x in guide_exon_list_updated]
            else:
                microhomology_sequence = "nan"
                exon_biotype_list  = "nan"
                
        ############ printing
        
            output =  unique_gene_id[0] + "\t" + guide_seq_with_ngg  + "\t" + guide_seq  + "\t" + guide_chr  + "\t" + guide_start  + "\t" + guide_stop  + "\t" + guide_strand  + "\t" + guide_uniq  \
                      + "\t" + str(guide_exon_list_updated) + "\t" + str(exon_biotype_list) + "\t" + str(guide_cutsite18) + "\t" + \
            "\t".join(cds_start_stop_cut_site)  + "\t" + "\t".join(exon_genomic_features)  + "\t" + microhomology_sequence  + "\n"
            
            gene_output_file_handle.write(output)
        
            #print output
        
        except Exception as e:
            
            log_output = gene_name + "\t" + str(e) + "\t" + "Guide_info_problem" + "\n"
            log_file_handle.write(log_output)
            
            
    log_file_handle.close()   
    gene_output_file_handle.close()


                
            
                
         
        
    
    
    
        
