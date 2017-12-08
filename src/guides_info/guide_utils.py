import re
import ast

def ReverseComplement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])

def generate_seq_for_microhomology_scoring (guide_start, guide_stop, guide_strand, guide_chr, chromosome_sequence_dict):
    chr_no = guide_chr
    sequence_start = int(guide_start) - 13; sequence_end = int(guide_stop) + 27; 
    
    chr_seq = chromosome_sequence_dict[chr_no]

    if guide_strand == "+" or guide_strand == "1":
        sequence_for_microhomology = chr_seq[sequence_start-1:sequence_end]
        
    if guide_strand == "-" or guide_strand == "-1":
        sequence_for_microhomology = ReverseComplement(chr_seq[sequence_start-1:sequence_end])

    return (sequence_for_microhomology)


def get_guide_seq_with_ngg (guide_start, guide_stop, guide_strand, guide_chr, chromosome_sequence_dict):
    chr_no = guide_chr
    sequence_start = int(guide_start)
    sequence_end = int(guide_stop)
    chr_seq = chromosome_sequence_dict[chr_no]
    
    if guide_strand == "+" or guide_strand == "1":
        new_seq_end = sequence_end + 3
        seq = chr_seq[sequence_start -1: new_seq_end]
        
    if guide_strand == "-" or guide_strand == "-1":
        new_seq_start = sequence_start -3
        seq = ReverseComplement(chr_seq[new_seq_start-1 :sequence_end])
        
    ##### test
    seq_ngg_regex = "[ATGC]{20}[ATGC]{1}GG"
    search_object_ngg  = re.search(seq_ngg_regex, seq)
    if search_object_ngg:
        ngg_outcome = seq
    else:
        ngg_outcome = "nan"

    return (ngg_outcome)


def get_guide_seq_with_ngg_test (guide_start, guide_stop, guide_strand, guide_chr, chromosome_sequence_dict):
    chr_no = guide_chr
    sequence_start = int(guide_start); sequence_end = int(guide_stop) ; 
    chr_seq = chromosome_sequence_dict[chr_no]
    
    if guide_strand == "+" or guide_strand == "1":
        new_seq_end = sequence_end + 3
        seq = chr_seq[sequence_start -1: new_seq_end]

    if guide_strand == "-" or guide_strand == "-1":
        new_seq_start = sequence_start -3
        seq = ReverseComplement(chr_seq[new_seq_start-1 :sequence_end])

    ##### test
    seq_ngg_regex = "[ATGC]{20}[ATGC]{1}GG"
    search_object_ngg  = re.search(seq_ngg_regex, seq)
    if search_object_ngg:
        ngg_outcome = seq
    else:
        ngg_outcome = "nan"

    return (seq, ngg_outcome)


def extract_exon_features_from_gtf (exon_id, exon_dict):
    if exon_id in exon_dict.exon_transcript_info:
        total_transcript_count =  exon_dict.exon_gene_info[exon_id]["total_transcript_per_gene"]   ###### gene_feature

        ###### exon_transcript info
        transcript_list = exon_dict.exon_transcript_info[exon_id]["transcript_list"]
        exon_rank_list = exon_dict.exon_transcript_info[exon_id]["exon_rank_list"]
        total_no_exon_transcript_list =  exon_dict.exon_transcript_info[exon_id]["total_exons_in_transcript"]
        
        return (total_transcript_count, str(transcript_list) , str(total_no_exon_transcript_list) ,str(exon_rank_list))
    else:
        return ("nan", "nan", "nan", "nan")
        
def extract_exon_features_from_gtf_for_exonlist (exon_list, exon_dict):
    transcript_list_all = []
    exon_rank_list_all = []
    total_no_exon_transcript_list_all = []

    for exon in exon_list:
        
        exon_biotype = exon_dict.exon_id_biotype[exon]["biotype"]
        exon_feautre = extract_exon_features_from_gtf(exon, exon_dict)
        
        
        if exon_biotype == "coding":
        
            
            transcript_list_all = transcript_list_all + ast.literal_eval(exon_feautre[1])
            total_no_exon_transcript_list_all = total_no_exon_transcript_list_all + ast.literal_eval(exon_feautre[2])
            exon_rank_list_all = exon_rank_list_all + ast.literal_eval(exon_feautre[3])
        
        else:
            no_of_transcript_covered = len(ast.literal_eval(exon_feautre[1]))
            list_exon_rank_insert  = [0 for x in range(no_of_transcript_covered)]
            exon_rank_list_all = exon_rank_list_all + list_exon_rank_insert
            
    
    total_transcipt_count = exon_feautre[0]  ####### since it is gene feature remains same for all the exons
    
    
    return (total_transcipt_count,str(transcript_list_all), str(total_no_exon_transcript_list_all), str(exon_rank_list_all)) 
        

def extract_exon_features_from_gtf_for_exonlist_modified (exon_list, exon_dict):
    
    
    ######## In case of no protein coding trsnacript
    if len(exon_list) == 0:
        
        total_transcipt_count = "nan"
        transcript_list_all = []
        exon_rank_list_all = []
        total_no_exon_transcript_list_all = []
        
    
    else:   
    
        transcript_list_all = []
        exon_rank_list_all = []
        total_no_exon_transcript_list_all = []
        
        
        
        
        for exon in exon_list:
            
            exon_biotype = exon_dict.exon_id_biotype[exon]["biotype"]
            exon_feautre = extract_exon_features_from_gtf(exon, exon_dict)
            
            
            if exon_biotype == "coding":
            
                
                transcript_list_all = transcript_list_all + ast.literal_eval(exon_feautre[1])
                total_no_exon_transcript_list_all = total_no_exon_transcript_list_all + ast.literal_eval(exon_feautre[2])
                exon_rank_list_all = exon_rank_list_all + ast.literal_eval(exon_feautre[3])
            
            else:
                no_of_transcript_covered = len(ast.literal_eval(exon_feautre[1]))
                list_exon_rank_insert  = [0 for x in range(no_of_transcript_covered)]
                exon_rank_list_all = exon_rank_list_all + list_exon_rank_insert
                
        
        total_transcipt_count = exon_feautre[0]  ####### since it is gene feature remains same for all the exons
    
    
    return (total_transcipt_count,str(transcript_list_all), str(total_no_exon_transcript_list_all), str(exon_rank_list_all)) 
        




def test_for_overlapping_exon_guide (exon_list, exon_dict):
    
    
    transcript_list_all = []
    
    for exon in exon_list:
        exon_feautre = extract_exon_features_from_gtf(exon, exon_dict)
        transcript_list_all = transcript_list_all + ast.literal_eval(exon_feautre[1])
        
        
    ######## To find common transcripts in transcript_list_all
    
    len_transcript_list_all = len(transcript_list_all)
    len_unique_transcript_list_all = len(list(set(transcript_list_all)))
    
    if len_unique_transcript_list_all == len_unique_transcript_list_all:
        
        outcome = "no_exon_overlap"
        
    else:
        
        outcome = "exon_overlap"
        
    
    return (outcome)
    
    
        
        
#print extract_exon_features_from_zv10_gtf85 ("ENSDARE00001181680")


def calculate_cutsite18_guide(guide_strand, guide_start, guide_stop):
    if guide_strand == "+" or guide_strand == "1":
        return (float(guide_stop) - 2)
    return (float(guide_start) + 2)


def calculate_cutsite18_dist_fom_exon_cds_start_stop_list (guide_exon,cutsite18, exon_dict) :

    if guide_exon in exon_dict.exon_transcript_info:
    
    
        exon_cds_start_stop_list = exon_dict.exon_transcript_info[guide_exon]["exon_cds_list"]
        exon_strand = exon_dict.exon_gene_info[guide_exon]["gene_strand"]
           
        
        
        ###### To remove nan
        
        exon_cds_start_stop_list_coding  = [x for x in exon_cds_start_stop_list if x!= "nan"]
        
        ####### To split cds for coding exons
        
        
        exon_cds_start_list = [float(x.split(":")[1].split("-")[0]) for x in exon_cds_start_stop_list_coding]
        exon_cds_stop_list = [float(x.split(":")[1].split("-")[-1]) for x in exon_cds_start_stop_list_coding]
        
        
        ###### Cutsite18
        cut_site_18 = float(cutsite18)
        
        
        ##### To add nan to the distance
        
        if "nan" in exon_cds_start_stop_list:
            
            no_transcript_covered = len(exon_cds_start_stop_list)
            
            dist_cds_start_cutsite_list_noncoding = ["nan" for x in range(no_transcript_covered)]
            dist_cds_stop_cutsite_list_noncoding = ["nan" for x in range(no_transcript_covered)]
            
            dist_cds_start_cutsite_list_coding = []
            dist_cds_stop_cutsite_list_coding = []
            
        
        else:
            
            dist_cds_start_cutsite_list_noncoding = []
            dist_cds_stop_cutsite_list_noncoding  = []
           
        
            ####### Calculating distances
                
            if exon_strand == "-" :  ####### real cds start is stop according to annotation 
                    
                dist_cds_start_cutsite_list_coding =  [abs(cut_site_18 - x) for x in exon_cds_stop_list]
                
                dist_cds_stop_cutsite_list_coding  =  [abs(cut_site_18 - x) for x in exon_cds_start_list]
                
            if exon_strand == "+" :
                    
                dist_cds_start_cutsite_list_coding  =  [abs(cut_site_18 - x) for x in exon_cds_start_list]
                dist_cds_stop_cutsite_list_coding  =  [abs(cut_site_18 - x) for x in exon_cds_stop_list]
                
        
        dist_cds_start_cutsite_list = dist_cds_start_cutsite_list_coding + dist_cds_start_cutsite_list_noncoding
        dist_cds_stop_cutsite_list = dist_cds_stop_cutsite_list_coding + dist_cds_stop_cutsite_list_noncoding
             
    
    return (dist_cds_start_cutsite_list, dist_cds_stop_cutsite_list)
    
        
def calculate_cutsite18_dist_fom_exon_cds_start_stop_for_exon_list (exon_list, cutsite18, exon_dict):
    
    dist_cds_start_cutsite_list_all = []
    dist_cds_stop_cutsite_list_all = []
    
    for exon in exon_list:
        
        distance_feature = calculate_cutsite18_dist_fom_exon_cds_start_stop_list(exon,cutsite18, exon_dict)
        dist_cds_start_cutsite_list_all = dist_cds_start_cutsite_list_all + distance_feature[0]
        dist_cds_stop_cutsite_list_all = dist_cds_stop_cutsite_list_all + distance_feature[1]
        
    
    return (str(dist_cds_start_cutsite_list_all), str(dist_cds_stop_cutsite_list_all))



def calculate_cutsite18_dist_fom_exon_cds_start_stop_for_exon_list_modified (exon_list, cutsite18, exon_dict):
    
    if len(exon_list) == 0:
        
        dist_cds_start_cutsite_list_all = []
        dist_cds_stop_cutsite_list_all = []
        
    else:
    
    
        dist_cds_start_cutsite_list_all = []
        dist_cds_stop_cutsite_list_all = []
        
        for exon in exon_list:
            
            distance_feature = calculate_cutsite18_dist_fom_exon_cds_start_stop_list(exon,cutsite18, exon_dict)
            dist_cds_start_cutsite_list_all = dist_cds_start_cutsite_list_all + distance_feature[0]
            dist_cds_stop_cutsite_list_all = dist_cds_stop_cutsite_list_all + distance_feature[1]
            
    
    return (str(dist_cds_start_cutsite_list_all), str(dist_cds_stop_cutsite_list_all))


##################### Notes

###### Exon list empty means no protein-coding transcript
####### Transcript list empty means two things : No protein-coding transcript and also that there is no coding exons so it is empty
###### Dist cds start-stop empty: no protein-coding transcript; Dist cds start-stop "nan" : no coding exons

####### Total score becomes "-20": when Exon list empty and transcript list is empty so in cases where there no coding transcript and also there is only a single non-coding exon for a transcript.
     
#def making_guide_file_with_info (guide_file_path,guide_file_info_directory):
#    
#    ######### Input_file
#    guide_file = guide_file_path
#    
#    
#    #########  making of  Output_file and log directory 
#    
#    gene_name = os.path.basename(guide_file).split(".")[0]
#    gene_output_file_name = gene_name + "_guides_info.txt"
#    
#    if os.path.exists(guide_file_info_directory):
#        
#        output_file_path = guide_file_info_directory
#        gene_output_file_name_path = os.path.join(output_file_path, gene_output_file_name)
#        gene_output_file_handle = open (gene_output_file_name_path, "w")        
#        
#        
#        ####### Eception file
#        
#        log_file_name = os.path.join(guide_file_info_directory,"log_file_exceptions.txt")
#        log_file_handle = open(log_file_name, "a")
#        
#    else:
#        
#        exit("Guide file info dir doesnot exist %s" %guide_file_info_directory)
#        
#        
#        
#    ######### Header
#    
#    header = ["Gene_id", "Guide_seq_with_ngg", "Guide_seq", "Guide_chr_no", "Guide_start","Guide_stop", "Guide_strand", "Guide_uniqness", \
#              "Guide_exon_list", "Exon_biotype_list", \
#            "cutsite18", "dist_cds_start_cutsite", "dist_cds_stop_cutsite",\
#              "Total_transcript_count", "Transcript_list", "Total_exons_in_transcript_list", "Exon_rank_in_transcript_list",\
#               "microhomology_seq"]
#    
#    header_output = "\t".join(header) + "\n"
#
#    ###### Writing header
#    gene_output_file_handle.write(header_output)
#    
#    
#    ################## Processing the file
#
#    guide_input_file = open(guide_file)
#    
#    for line in guide_input_file:
#        
#        l = line.strip("\n").split("\t")
#        guide_chr = l [0]; guide_start = l[1]; guide_stop = l[2]; guide_strand = l[3]; guide_seq = l[4]; guide_offtarget_profile = l[10]
#        
#        
#    
#    ###### Calculate cutsite18
#    
#        try:
#            guide_cutsite18 =  calculate_cutsite18_guide(guide_strand,guide_start,guide_stop)  
#        except Exception as e:
#            log_output = gene_name+ "\t" + str(e) + "\t" + "guide_cutsite18" + "\n"
#            log_file_handle.write(log_output)
#            
#    
#    
#    ##### Searching for exon_id
#        
#        guide_exon_list = []
#        for exon_location in  exon_dict.location_exon_id:
#            exon_chr = exon_location.split(":")[0]
#            exon_start_stop = exon_location.split(":")[1]
#            
#            if guide_chr == exon_chr: ##### 1st step
#                exon_start = float(exon_start_stop.split("-")[0])
#                exon_stop = float(exon_start_stop.split("-")[1])
#                ## Guide should be between exon start/stop
#                if exon_start <= exon_stop:
#                    if float(guide_cutsite18) >= exon_start and float(guide_cutsite18) <= exon_stop:
#                        guide_exon =  exon_dict.location_exon_id[exon_location]
#                        guide_exon_list = guide_exon_list + guide_exon
#                    else: ###### 
#                        pass        
#                else:
#                    e = "Exon_start > Exon_stop"
#                    log_output = gene_name + "\t" + str(e) + "\t" + "Exon_id_search_problem" + "\n"
#                    log_file_handle.write(log_output)
#            else:
#                pass
#    
#        ####### If exon list is empty because there is no protein-coding transcript then the gene_id_list would also be empty
#        ######## there are two cases here: one when exon list is empty becaus ethere is no protein coding transcript and also when exon list is not empty there ,
#        #### after overlapping removal it still would contain gene from input file
#        
#        ##### there are three cases: 1. Non-protein coding transcript (guide exon list  = 0 )2. Overlapping genes 3. Both  (guide_exon_list_updated = 0)
#        
#        
#        if len(guide_exon_list) !=0 :                         #### It could belong to case 2 and case 3
#            gene_id_list = [ exon_dict.exon_gene_info[x]["gene_id"] for x in guide_exon_list]  ###### To detect problem in gene_id
#            ######### To remove overlapping genes
#            guide_exon_list_remove = []
#            
#            for exon in guide_exon_list:
#                gene_id =  exon_dict.exon_gene_info[exon]["gene_id"]
#                if gene_id != gene_name:
#                    guide_exon_list_remove.append(exon)
#                    
#            guide_exon_list_updated = [x for x in guide_exon_list if x not in guide_exon_list_remove]
#            gene_id_list_updated = [ exon_dict.exon_gene_info[x]["gene_id"] for x in guide_exon_list_updated]
#            unique_gene_id = list(set(gene_id_list_updated))
#            
#            if len(gene_id_list_updated) != 0 :
#                if len(unique_gene_id) !=1 and unique_gene_id[0] != gene_name:
#                    e = "Gene_id not match input file after overlapping"
#                    log_output = gene_name+ "\t" + str(e) + "\t" + "gene_id_match_after_overlap" 
#                    log_file_handle.write(log_output)
#        else:
#            guide_exon_list_updated = guide_exon_list
#            
#            
#    ######### Getting guide info
#    
#        try:
#        ####### Guide_uniqueness
#                    
#            offtarget_profile_list = guide_offtarget_profile.split(",")
#            guide_uniq = offtarget_profile_list[1]
#                    
#            ###### seq_with_ngg
#            guide_seq_with_ngg =  get_guide_seq_with_ngg(guide_start, guide_stop, guide_strand, guide_chr)
#            
#             ########## Extracting feature
#            cds_start_stop_cut_site =  calculate_cutsite18_dist_fom_exon_cds_start_stop_for_exon_list_modified(guide_exon_list_updated,guide_cutsite18)
#            exon_genomic_features =  extract_exon_features_from_gtf_for_exonlist_modified(guide_exon_list_updated)
#            
#            if len(guide_exon_list_updated) != 0:
#                microhomology_sequence =  generate_seq_for_microhomology_scoring(guide_start, guide_stop, guide_strand, guide_chr)
#                exon_biotype_list = [ exon_dict.exon_id_biotype[x]["biotype"] for x in guide_exon_list_updated]
#            else:
#                microhomology_sequence = "nan"
#                exon_biotype_list  = "nan"
#                
#        ############ printing
#        
#            output =  unique_gene_id[0] + "\t" + guide_seq_with_ngg  + "\t" + guide_seq  + "\t" + guide_chr  + "\t" + guide_start  + "\t" + guide_stop  + "\t" + guide_strand  + "\t" + guide_uniq  \
#                      + "\t" + str(guide_exon_list_updated) + "\t" + str(exon_biotype_list) + "\t" + str(guide_cutsite18) + "\t" + \
#            "\t".join(cds_start_stop_cut_site)  + "\t" + "\t".join(exon_genomic_features)  + "\t" + microhomology_sequence  + "\n"
#            
#            gene_output_file_handle.write(output)
#        
#            #print output
#        
#        except Exception as e:
#            
#            log_output = gene_name + "\t" + str(e) + "\t" + "Guide_info_problem" + "\n"
#            log_file_handle.write(log_output)
#            
#            
#    log_file_handle.close()   
#    gene_output_file_handle.close()


                
            
                
