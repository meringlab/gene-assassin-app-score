import os, sys
import ast
import math
import gzip
from math import exp
from re import findall

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


def make_var_dict (variation_file_name,dir_path):
    
    variation_file_path = search_file_path(variation_file_name,dir_path)
    snv_dict = {}
    
    if variation_file_path!= "":

        variation_input_file_handle = gzip.open(variation_file_path,"rb")
        
        for line in variation_input_file_handle:
            if line[0]!= "#":
                l = line.strip("\n").split("\t")
                
                chr_no = l[0]; vartiation_type = l[2]; vartiation_start = l[3]; vartiation_stop = l[4]; descript = l[-1]
                
                if vartiation_type == "SNV":
                    if chr_no not in snv_dict:
                        snv_dict[chr_no] = set()
                    #although it's a FU SNP, it might have different start/stop positions ?!?
                    snv_dict[chr_no].add(vartiation_start)
                    snv_dict[chr_no].add(vartiation_stop)
    else:
    
        print "Error : Variation file path not found"
    
    
    return (snv_dict)


def make_transcript_cds_info_dict(transcript_cds_file_name, dir_path):
    transcript_cds_file_path = search_file_path(transcript_cds_file_name, dir_path)
    if not transcript_cds_file_path:
        raise Exception("Error: Transcript CDS file path not found %s" % transcript_cds_file_name)

    transcript_cds_info = {}

    with open(transcript_cds_file_path) as transcript_cds_input_file_handle:
        transcript_cds_input_file_handle.readline()  # throw away the header

        for line in transcript_cds_input_file_handle:
            l = line.strip().split('\t')
            transcript_id = l[0]
            first_exon_cds = l[1]
            last_exon_cds = l[2]
            cds_start = l[4]
            cds_stop = l[5]
            cds_length = l[-1]
            # transcript_strand = l[-2]

            if transcript_id not in transcript_cds_info:
                transcript_cds_info[transcript_id] = {}

            transcript_cds_info[transcript_id]["cds_start_stop_len"] = [cds_start, cds_stop, cds_length]
            transcript_cds_info[transcript_id]["cds_first_last_exon_rank"] = [first_exon_cds, last_exon_cds]

    return transcript_cds_info


def make_protein_dict (name_protein_dir, dir_path):
    
    protein_path = search_dir_path (name_protein_dir, dir_path)
    protein_domain_info_dict = {}
    
    if protein_path != "":
        
        all_gene_list_domains = os.listdir(protein_path)
        for input_file in all_gene_list_domains:
    
            input_file_path = os.path.join(protein_path, input_file)
            input_file_handle = open(input_file_path)
            
            for line in input_file_handle:
                
                l = line.strip("\n").split("\t")
                if l[0] == '#>P':
                    gene_id = l[1]
                ################## To test the presence of domains, if it is not found in dict that means no domains   
                if not l[0].startswith("#"):
                    chr_name = l[0]; domain_start = l[1]; domain_stop = l[2]; domain_name = l[4]; domain_id = l[-2]
                    #print chr_name, domain_start, domain_stop, domain_name, domain_id
                    if int(domain_start) <= int(domain_stop): 
                        pass
                    else:
                        print "domain_start > domain_stop "
                    domain_location = chr_name + ":" + domain_start + "-" + domain_stop
                    domain_id_desc = domain_name + ":" + domain_id
                    if gene_id not in protein_domain_info_dict:
                        protein_domain_info_dict[gene_id] = {}
                        
                    if domain_location not in protein_domain_info_dict[gene_id]:
                        protein_domain_info_dict[gene_id][domain_location] = domain_id_desc
                        
                        
    else:
        print "Error: Protein dir path not found"
        
    return (protein_domain_info_dict)
    
###################### Genome Context_scores

#######1

def calculate_proximity_to_CDS (transcript_id, cutsite18,transcript_cds_info_dict): ################# if within 10%CDS_start/stop then penalised
    
    transcript_cds_data = transcript_cds_info_dict[transcript_id]["cds_start_stop_len"]
    
    cds_start = float(transcript_cds_data[0])
    cds_stop = float(transcript_cds_data[1])
    cds_length = float(transcript_cds_data[2])
    
    cutsite18_float = float(cutsite18)
    
    cds_len_10_percent = round(0.1*cds_length)
    
    if cds_start < cds_stop:     #################### forward strand
        
        lower_limit = cds_start + cds_len_10_percent
        upper_limit = cds_stop - cds_len_10_percent
        
    else:  ############ reverse strand ()
        
        lower_limit = cds_stop + cds_len_10_percent
        upper_limit = cds_start - cds_len_10_percent
        
        
        
    if  cutsite18_float > lower_limit and cutsite18_float < upper_limit: ############ safe region to be not too near cds start or stop
        
        penalty_score = 0
        
    else:
        
        penalty_score = -1
        
    
    return (penalty_score)

def calculate_proximity_to_CDS_for_transcript_list_modified (transcript_list_input, cutsite18,transcript_cds_info_dict): ################# if within 10%CDS_start/stop then penalised
    

    
    transcript_list = ast.literal_eval(transcript_list_input)
    
    if len(transcript_list) == 0:
        avg_penalty_score_prox_CDS = -1
    
    else:
        
        penalty_score_prox_CDS = []
        
        
        for i in range(0,len(transcript_list)):
                
            transcript = transcript_list[i]
            
            guide_cutsite_proximity_CDS = calculate_proximity_to_CDS(transcript_id = transcript, cutsite18 = cutsite18, transcript_cds_info_dict = transcript_cds_info_dict)
            
            penalty_score_prox_CDS.append(guide_cutsite_proximity_CDS)
            
        avg_penalty_score_prox_CDS = round( float(sum(penalty_score_prox_CDS))/float(len(penalty_score_prox_CDS)), 2)
        
    
    return (avg_penalty_score_prox_CDS)




######### 2
def calculate_proximity_splice_site (dist_cutsite_exon_cds_start, dist_cutsite_exon_cds_stop, exon_rank_in_transcript, transcript_id, transcript_cds_info_dict):
    
        
    transcript_cds_exon_ranks = transcript_cds_info_dict[transcript_id]["cds_first_last_exon_rank"]
    
    input_exon_rank = str(exon_rank_in_transcript)
    
    if input_exon_rank in transcript_cds_exon_ranks:       ########### first and last exon
        
        if float(dist_cutsite_exon_cds_start) <= 20 or float(dist_cutsite_exon_cds_stop) <= 20:
            
            proximity_penalty = -1.5
        
        else:
            
            proximity_penalty = 0
            
    else:
        
        if float(dist_cutsite_exon_cds_start) <= 15 or float(dist_cutsite_exon_cds_stop) <= 15:
            
            proximity_penalty = -1
        
        else:
            
            proximity_penalty = 0
            
            
            
    return (proximity_penalty)

def calculate_proximity_splice_site_for_exon_rank_list_modified (exon_rank_list_input, dist_cutsite_exon_cds_start_list_input, dist_cutsite_exon_cds_stop_list_input, transcript_list_input, transcript_cds_info_dict):
    
    
    exon_list = ast.literal_eval(exon_rank_list_input)
    transcript_list = ast.literal_eval(transcript_list_input)
    dist_cutsite_exon_cds_start_list  = ast.literal_eval(dist_cutsite_exon_cds_start_list_input)
    dist_cutsite_exon_cds_stop_list = ast.literal_eval(dist_cutsite_exon_cds_stop_list_input)
    
    
    if len(exon_list) == len(transcript_list): ###### no, non-coding exon
        
        penalty_score_prox_splicesite = []
    
        for i in range(0,len(exon_list)):
                
            exon_rank = exon_list[i]
            transcript = transcript_list[i]
            dist_cutsite_exon_cds_start = dist_cutsite_exon_cds_start_list[i]
            dist_cutsite_exon_cds_stop = dist_cutsite_exon_cds_stop_list [i]
                
            guide_cutsite_proximity_splicesite = calculate_proximity_splice_site(dist_cutsite_exon_cds_start = dist_cutsite_exon_cds_start, dist_cutsite_exon_cds_stop = dist_cutsite_exon_cds_stop,\
                                                  exon_rank_in_transcript = exon_rank, transcript_id = transcript, transcript_cds_info_dict = transcript_cds_info_dict)
            
            
            penalty_score_prox_splicesite.append(guide_cutsite_proximity_splicesite)
            
        avg_penalty_score_prox_splicesite =  round(float(sum(penalty_score_prox_splicesite))/float(len(penalty_score_prox_splicesite)) ,2)
            
        
        
    else:
        
        non_noncoding_exons = len(exon_list) - len(transcript_list)
        
        exon_list_modified = [x for x in exon_list if x!= 0]
        dist_cutsite_exon_cds_start_list_modified  = [x for x in dist_cutsite_exon_cds_start_list if x!= "nan"]
        dist_cutsite_exon_cds_stop_list_modified  = [x for x in dist_cutsite_exon_cds_stop_list if x!= "nan"]
    
        
        penalty_score_prox_splicesite = []
        
        penalty_score_prox_splicesite_nc_exons = [-1.5 for x in range(non_noncoding_exons)] ###### penalising for no. of non-coding exons
        penalty_score_prox_splicesite = penalty_score_prox_splicesite + penalty_score_prox_splicesite_nc_exons
        
        
        if len(transcript_list) != 0:   ######## Case where only one exon is hit by guide and that is non-coding
            
            for i in range(0,len(exon_list_modified)):
                    
                exon_rank = exon_list_modified[i]
                transcript = transcript_list[i]
                dist_cutsite_exon_cds_start = dist_cutsite_exon_cds_start_list_modified[i]
                dist_cutsite_exon_cds_stop = dist_cutsite_exon_cds_stop_list_modified[i]
                    
                guide_cutsite_proximity_splicesite = calculate_proximity_splice_site(dist_cutsite_exon_cds_start = dist_cutsite_exon_cds_start, dist_cutsite_exon_cds_stop = dist_cutsite_exon_cds_stop,\
                                                      exon_rank_in_transcript = exon_rank, transcript_id = transcript)
                
                
                penalty_score_prox_splicesite.append(guide_cutsite_proximity_splicesite)
            
        avg_penalty_score_prox_splicesite =  round(float(sum(penalty_score_prox_splicesite))/float(len(penalty_score_prox_splicesite)) ,2)
            
    
    return (avg_penalty_score_prox_splicesite)
    
    
##### 3

def calculate_exon_ranking_score_modified (exon_ranks) :
    
    
    ##### there can be no cases of nan
    
    exon_ranks_unstring = ast.literal_eval(exon_ranks)
    exon_rank_list = [float(x) for x in exon_ranks_unstring]
    
    
    for x in exon_rank_list:
        if x < 0 or math.isnan(x):
            raise Exception 
    
    
    
    combined_exon_rank = round(sum(exon_rank_list)/len(exon_rank_list),2)
    
    if combined_exon_rank == 1:
        exon_rank_score = 1
    
    if combined_exon_rank >1 and combined_exon_rank <= 3:
        exon_rank_score = 2
    
    if combined_exon_rank > 3 and combined_exon_rank <= 5:
        exon_rank_score = 1
        
    if combined_exon_rank >5 and combined_exon_rank <= 7:
        exon_rank_score = 0.5
        
    if combined_exon_rank >7 :
        exon_rank_score = 0.25
    
    if combined_exon_rank < 1:
        exon_rank_score = 0
    
    
    return (exon_rank_score)
    
    
#### 4

def calculate_transcript_coverage_score_modified (transcript_ids,transcript_count):
    
    
    transcript_ids_unstring = ast.literal_eval(transcript_ids)
    
    transcript_covered = len(transcript_ids_unstring)
    transcript_count = float(transcript_count)
    
    
    assert(transcript_covered > 0) , "Transcript_list_empty"
    assert(transcript_covered <= transcript_count) , "Transcript_count < transcipt_id_list"
    
    if "nan" in transcript_ids_unstring:
        raise Exception
        
    
    transcript_covered_ratio = round(transcript_covered/transcript_count,2)
    
    if transcript_covered_ratio == 0:   ###### cases where transcript list is empty
        transcript_covered_score = 0
    
    if transcript_covered_ratio > 0 and transcript_covered_ratio <= 0.25:
        
        transcript_covered_score = 0.5
        
    if transcript_covered_ratio > 0.25 and transcript_covered_ratio <= 0.5:
        
        transcript_covered_score = 1
        
    if transcript_covered_ratio > 0.5 and transcript_covered_ratio <= 0.75:
        
        transcript_covered_score = 1.5
        
    if transcript_covered_ratio > 0.75 and transcript_covered_ratio <= 1:
        
        transcript_covered_score = 2.5
     
    
    return (transcript_covered_score)

########################################### Protein Domain

def calculate_score_protein_domains (cutsite_18, gene_id,protein_domain_info_dict):
    
    cutsite_18_int = float(cutsite_18)
    
    if gene_id not in protein_domain_info_dict:
        
        domain_score = 0
    
    
    #### to find if the guide covers the domain   
    else:
        
        domain_location = protein_domain_info_dict[gene_id].keys()
        
        domain_start_stop_list = [x.split(":")[1] for x in domain_location]
        
        
        for location in domain_start_stop_list:
            
            domain_start = float(location.split("-")[0])
            domain_stop = float(location.split("-")[1])
            
            
            if domain_start > domain_stop:
                
                print "Error, start > stop"
                
                domain_score = "Error"
                
            else:
                
                
                if cutsite_18_int >= domain_start and cutsite_18_int <= domain_stop:
                    
                    domain_score = 2.5
                     
                    break   ############## to stop searching once it finds it otherwise score becomes 0 again
                
                else:
                    domain_score = 0
                    
                         
    
    return (domain_score)
    
    
    
###### Microhomology Scores


def  get_microhomology_out_of_frame_score (seq):

    length_weight=20.0
    left=30 # Insert the position expected to be broken.
    right=len(seq)-int(left)
    
    
    output_list_before_duplication = []
    
    for k in range(2,left)[::-1]:
        for j in range(left,left+right-k+1): 
            for i in range(0,left-k+1):
                if seq[i:i+k]==seq[j:j+k]:
                    length=j-i
                    
                    sequence_list_before_duplication = [seq[i:i+k],str(i),str(i+k),str(j),str(j+k),str(length)]
                    
                    output_list_before_duplication.append(sequence_list_before_duplication)
                    
    
    
    ## After searching out all microhomology patterns, duplication should be removed!!
    if output_list_before_duplication !="":
        
        sum_score_3=0
        sum_score_not_3=0
        
        for i in range(len(output_list_before_duplication)):
            
            n=0
            score_3=0
            score_not_3=0
            
            line= output_list_before_duplication[i]   ########### first line
            
            scrap=line[0]
            left_start=int(line[1])
            left_end=int(line[2])
            right_start=int(line[3])
            right_end=int(line[4])
            length=int(line[5])
    
        
            for j in range(i):
            
                line_ref= output_list_before_duplication[j]   ######### line to be comapred
                
                left_start_ref=int(line_ref[1])
                left_end_ref=int(line_ref[2])
                right_start_ref=int(line_ref[3])
                right_end_ref=int(line_ref[4])
                
                
                if (left_start >= left_start_ref) and (left_end <= left_end_ref) and (right_start >= right_start_ref) and (right_end <= right_end_ref):
                    
                    if (left_start - left_start_ref)==(right_start - right_start_ref) and (left_end - left_end_ref)==(right_end -right_end_ref):
                        
                        n+=1
                    
                else: pass
            
            
            if n == 0:
                
                length_factor = round(1/exp((length)/(length_weight)),3)
                num_GC=len(findall('G',scrap))+len(findall('C',scrap))
                common_score = 100*length_factor*((len(scrap)-num_GC)+(num_GC*2))
                
                if (length % 3)==0:
            
                    score_3 = common_score
                    
                elif (length % 3)!=0:
                
                    
                    score_not_3 = common_score
                    
                
                
            
            sum_score_3+=score_3
            sum_score_not_3+=score_not_3
            
        
    
        microhomology_score = sum_score_3+sum_score_not_3
        
        
        if microhomology_score!= 0:
            out_of_frame_score = round((sum_score_not_3)*100/(microhomology_score),2)
            
        else:
            out_of_frame_score = "nan"
        
        
        
    return (microhomology_score,out_of_frame_score)



def calculate_microhomology_score (seq):
    
    microhomology_out_frame_score =  get_microhomology_out_of_frame_score(seq = seq)
    microhomology_score = microhomology_out_frame_score[0]
    out_of_frame_score = microhomology_out_frame_score[1]
    
    if out_of_frame_score == "nan":
        calculated_micrhomology_score = "nan"
    else:
        if  out_of_frame_score > 66.0 :
            calculated_micrhomology_score = 0.5
        else:
            calculated_micrhomology_score = 0
            
    return (calculated_micrhomology_score)





######### SNP Score

def calculate_snp_score(guide_chr_input, guide_cutsite_18_input,snv_dict):
    
    guide_chr = str(guide_chr_input)
    guide_cutsite_18 = str(guide_cutsite_18_input)

    if guide_cutsite_18 in snv_dict[guide_chr]:
        s = 0   ##### "snp found" ##### penalised
    else:
        s = 2.5  ##### no_snp_so_good
         
    return s

    
    

    
