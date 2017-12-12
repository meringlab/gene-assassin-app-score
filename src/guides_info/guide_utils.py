import re
import ast


def extract_exon_features_from_gtf (exon_id, exon_dict):
    if exon_id in exon_dict.exon_transcript_info:
        total_transcript_count =  exon_dict.exon_gene_info[exon_id]["total_transcript_per_gene"]   ###### gene_feature

        transcript_list = exon_dict.exon_transcript_info[exon_id]["transcript_list"]
        exon_rank_list = exon_dict.exon_transcript_info[exon_id]["exon_rank_list"]
        total_no_exon_transcript_list =  exon_dict.exon_transcript_info[exon_id]["total_exons_in_transcript"]
        
        return (total_transcript_count, str(transcript_list) , str(total_no_exon_transcript_list) ,str(exon_rank_list))
    else:
        return ("nan", "nan", "nan", "nan")


def extract_exon_features_from_gtf_for_exonlist_modified (exon_list, exon_dict):
    if len(exon_list) == 0:  # no protein coding trsnacript
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

    # find common transcripts in transcript_list_all
    len_transcript_list_all = len(transcript_list_all)
    len_unique_transcript_list_all = len(list(set(transcript_list_all)))
    
    if len_unique_transcript_list_all == len_unique_transcript_list_all:
        outcome = "no_exon_overlap"
    else:
        outcome = "exon_overlap"

    return (outcome)


#print extract_exon_features_from_zv10_gtf85 ("ENSDARE00001181680")

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