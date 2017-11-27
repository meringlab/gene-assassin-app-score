import os
import sys
import gzip
import re


gtf_file_path = "/home/neha/Projects/Christian_Mosimann/Crispr_Project/v85_2017/danio_rerio/v_85/Raw_data_files/Danio_rerio.GRCz10.85.gtf.gz"

################################## Making new output file name and directory 
try :
    os.path.exists(gtf_file_path)
except Exception as e:
    print e
    exit("\t ... gtf file path does not exist %s."% gtf_file_path)
    

        
gtf_file_path_modified = gtf_file_path.replace("Raw_data_files", "Proesssed_data_files")
new_basefile_path =  "/".join(gtf_file_path_modified.split("/")[0:-1])

try:
    if  not os.path.exists(new_basefile_path):
        os.makedirs(new_basefile_path)

except Exception as e:
    print e
    exit("\t ... Cannot create a directory %s."% new_basefile_path)

        
gtf_file_name = os.path.basename(gtf_file_path)
new_file_name = "parsed" + "_" +gtf_file_name.strip(".gz") + ".txt"
new_file_path = os.path.join(new_basefile_path, new_file_name)


############################################## Processing gtf file

###### Input gtf file

input_file_handle = gzip.open(gtf_file_path, 'rb') 



########### Regular expression to extract gene features

reg_ex_string_gene = ".*?gene_id\s*\"(.*?)\".*"
reg_ex_string_transcript = ".*?gene_id\s*\"(.*?)\".*?transcript_id\s*\"(.*?)\".*?transcript_biotype\s*\"(.*?)\".*"

reg_ex_string_exon = ".*?transcript_id\s*\"(.+?)\".*?exon_number\D*(\d*).*?transcript_biotype\s*\"(.*?)\".*?exon_id\s*\"(.+?)\".*"
reg_ex_string_cds = ".*?transcript_id\s*\"(.+?)\".*?exon_number\D*(\d*).*"

######################## Reading the gtf file 
   
transcript_exon_no_exon_cds_dict = {}
gene_id_transcript_dict = {}
transcript_id_gene_dict = {}
gene_location_dict = {}


for line in input_file_handle:
    
    if line[0]!= "#":
        
        
        l = line.strip("\n").split("\t")
      
            
        gene_feature = l[2]
        chr_no = l[0]; start_position = l[3] ; stop_position = l[4]
        location = chr_no + ":" + start_position + "-" +stop_position
        strand = l[6]
        
        key_value_list = l[8].split(";")
        key_value_list_cat = "\t".join(key_value_list)
        
        
        
    ############## Making gene location_dict
        
        if gene_feature == "gene":
            
            search_obj_gene = re.search(reg_ex_string_gene,key_value_list_cat)
            gene_name = search_obj_gene.group(1)
            
            #print gene_name
            
            if gene_name not in gene_location_dict:
                 gene_location_dict[gene_name] = {}
                
            gene_location_dict[gene_name]["gene_location"] = location
            gene_location_dict[gene_name]["strand"] = strand
        
        
        
        
            
        
        ######## Searching for gene_features
        
        
        if gene_feature == "transcript":
            
            search_obj_transcript = re.search(reg_ex_string_transcript,key_value_list_cat)
                
            #print key_value_list_cat
            
            if search_obj_transcript:
               
                gene_id = search_obj_transcript.group(1)
                transcript_id = search_obj_transcript.group(2)
                transcript_type = search_obj_transcript.group(3)
                 
                #print gene_id, transcript_id, transcript_type
                
            
                if transcript_type == "protein_coding":
                    
                    
                    ############# Making transcript_gene_dict
                   
                    if transcript_id not in transcript_id_gene_dict:
                        transcript_id_gene_dict[transcript_id] = gene_id
                    
                    ############ Making a gene transcriptdict
                    
                    if gene_id not in gene_id_transcript_dict:
                        
                        gene_id_transcript_dict[gene_id] = []
                        
                    if transcript_id not in gene_id_transcript_dict[gene_id]:
                        gene_id_transcript_dict[gene_id].append(transcript_id)
                        
                
                #################### Doing exon search
                
                
        if gene_feature == "exon":
            
            #print key_value_list_cat
            
            search_obj_exon = re.search(reg_ex_string_exon,key_value_list_cat)
            transcript_id_exon = search_obj_exon.group(1)
            exon_number_exon = search_obj_exon.group(2)
            transcript_biotype_exon = search_obj_exon.group(3)
            exon_id = search_obj_exon.group(4)
            
            #print transcript_id_exon, exon_number_exon, exon_id, transcript_biotype_exon#
            if transcript_biotype_exon == "protein_coding":
            
                if transcript_id_exon not in transcript_exon_no_exon_cds_dict:
                    transcript_exon_no_exon_cds_dict[transcript_id_exon] = {}
                    
                if exon_number_exon not in transcript_exon_no_exon_cds_dict[transcript_id_exon]:
                    
                    transcript_exon_no_exon_cds_dict[transcript_id_exon][exon_number_exon] = {}
                    transcript_exon_no_exon_cds_dict[transcript_id_exon][exon_number_exon]["exon_id"] = {}
                    transcript_exon_no_exon_cds_dict[transcript_id_exon][exon_number_exon]["cds_location"] = []
                
                transcript_exon_no_exon_cds_dict[transcript_id_exon][exon_number_exon]["exon_id"][exon_id] = location
                
                    
                ################# Doing CDS search
                
                    
        if gene_feature == "CDS":
            
            #print key_value_list_cat
            
            search_obj_cds = re.search(reg_ex_string_cds,key_value_list_cat)
              
            transcript_id_cds = search_obj_cds.group(1)
            exon_number_cds = search_obj_cds.group(2)
            
            #print transcript_id_cds, exon_number_cds
            
            
            
            if transcript_id_cds not in transcript_exon_no_exon_cds_dict:
                transcript_exon_no_exon_cds_dict[transcript_id_cds] = {}
            
            if exon_number_cds not in transcript_exon_no_exon_cds_dict[transcript_id_cds]:
                
                transcript_exon_no_exon_cds_dict[transcript_id_cds][exon_number_cds] = {}
                transcript_exon_no_exon_cds_dict[transcript_id_cds][exon_number_cds]["exon_id"] = {}
                transcript_exon_no_exon_cds_dict[transcript_id_cds][exon_number_cds]["cds_location"] = []
            
            transcript_exon_no_exon_cds_dict[transcript_id_cds][exon_number_cds]["cds_location"].append(location)
            
            
##################### Printing the output 


output_file_handle = open(new_file_path, "w")                   
                    
                    
header = ["Gene_id", "Gene_location", "Gene_strand","Total_transcript", "Transcript_id", "Total_Exons","Exon_id", "Exon_location", "Exon_rank","CDS_location"]
output = "\t".join(header) + "\n"

output_file_handle.write(output)

for transcript in transcript_id_gene_dict:
    
    gene = transcript_id_gene_dict[transcript]
    gene_location = gene_location_dict[gene]["gene_location"]
    gene_strand = gene_location_dict[gene]["strand"]
    
    total_transcript_count = len(gene_id_transcript_dict[gene])
    total_exons_in_transcript = len(transcript_exon_no_exon_cds_dict[transcript])
    
    for exon_rank in transcript_exon_no_exon_cds_dict[transcript]:
        
        
        exon = transcript_exon_no_exon_cds_dict[transcript][exon_rank]["exon_id"].keys()[0]
        exon_location = transcript_exon_no_exon_cds_dict[transcript][exon_rank]["exon_id"].values()[0]
        
        #print gene, transcript, exon, total_transcript_count, exon_rank, total_exons_in_transcript
    
    
        
        cds_location_list = transcript_exon_no_exon_cds_dict[transcript][exon_rank]["cds_location"]
        
        if len(cds_location_list) > 0:
            
            cds_location = cds_location_list[0]
        
        else:
            cds_location = "nan"
            
            
        output = gene + "\t" + gene_location + "\t" + gene_strand + "\t" + str(total_transcript_count) + "\t" +transcript + "\t" + str(total_exons_in_transcript)  + "\t" +  exon + "\t" +exon_location + "\t" + str(exon_rank) + "\t" + cds_location + "\n"
        
        
        output_file_handle.write(output)
    
output_file_handle.close()

sys.stdout.write('\t... saved at %s \n' % new_file_path)
        
        
        
       
        

