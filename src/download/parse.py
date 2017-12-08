import os
import sys
import json
import gzip
import re
import logging


def parse_gtf(params):
    gtf_file_name = os.path.basename(params['base_url_gtf'])
    gtf_file_path = os.path.join('output', params['ensembl_release'], params['species_name'], 'Raw_data_files',
                                 gtf_file_name)
    logging.info('parsing %s', gtf_file_path)

    if not os.path.exists(gtf_file_path):
        exit("gtf file does not exist %s." % gtf_file_path)

    gtf_file_path_modified = gtf_file_path.replace("Raw_data_files", "Processed_data_files")
    new_basefile_path = os.path.dirname(gtf_file_path_modified)
    new_file_name = 'parsed_' + gtf_file_name.replace('.gz', '.txt')
    new_file_path = os.path.join(new_basefile_path, new_file_name)

    try:
        if not os.path.exists(new_basefile_path):
            os.makedirs(new_basefile_path)
    except Exception as e:
        logging.exception("Cannot create a directory %s. %s" % new_basefile_path, e)
        exit()

    # Processing gtf file

    ## Input gtf file

    input_file_handle = gzip.open(gtf_file_path, 'rb')

    ### Regular expression to extract gene features

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

        if line[0] != "#":

            l = line.strip("\n").split("\t")
            # example: '4	ensembl_havana	gene	6733	52120	.	-	.	gene_id "ENSDARG00000104632"; gene_version "1"; gene_name "si:ch73-252i11.3"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; havana_gene "OTTDARG00000037780"; havana_gene_version "1";'

            gene_feature = l[2]
            chr_no = l[0];
            start_position = l[3];
            stop_position = l[4]
            location = chr_no + ":" + start_position + "-" + stop_position
            strand = l[6]

            key_value_list = l[8].split(";")
            key_value_list_cat = "\t".join(key_value_list)

            ############## Making gene location_dict

            if gene_feature == "gene":

                search_obj_gene = re.search(reg_ex_string_gene, key_value_list_cat)
                gene_name = search_obj_gene.group(1)

                # print gene_name

                if gene_name not in gene_location_dict:
                    gene_location_dict[gene_name] = {}

                gene_location_dict[gene_name]["gene_location"] = location
                gene_location_dict[gene_name]["strand"] = strand

            ######## Searching for gene_features


            if gene_feature == "transcript":

                search_obj_transcript = re.search(reg_ex_string_transcript, key_value_list_cat)

                # print key_value_list_cat

                if search_obj_transcript:

                    gene_id = search_obj_transcript.group(1)
                    transcript_id = search_obj_transcript.group(2)
                    transcript_type = search_obj_transcript.group(3)

                    # print gene_id, transcript_id, transcript_type


                    if transcript_type == "protein_coding":

                        ############# Making transcript_gene_dict

                        if transcript_id not in transcript_id_gene_dict:
                            transcript_id_gene_dict[transcript_id] = gene_id

                        ############ Making a gene transcriptdict

                        if gene_id not in gene_id_transcript_dict:
                            gene_id_transcript_dict[gene_id] = []

                        if transcript_id not in gene_id_transcript_dict[gene_id]:
                            gene_id_transcript_dict[gene_id].append(transcript_id)

                else:
                    logging.warning('found no transcript in %s, record: %s', l[8], line)
                    #################### Doing exon search

            if gene_feature == "exon":

                # print key_value_list_cat

                search_obj_exon = re.search(reg_ex_string_exon, key_value_list_cat)
                transcript_id_exon = search_obj_exon.group(1)
                exon_number_exon = search_obj_exon.group(2)
                transcript_biotype_exon = search_obj_exon.group(3)
                exon_id = search_obj_exon.group(4)

                # print transcript_id_exon, exon_number_exon, exon_id, transcript_biotype_exon#
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
                search_obj_cds = re.search(reg_ex_string_cds, key_value_list_cat)

                transcript_id_cds = search_obj_cds.group(1)
                exon_number_cds = search_obj_cds.group(2)

                if transcript_id_cds not in transcript_exon_no_exon_cds_dict:
                    transcript_exon_no_exon_cds_dict[transcript_id_cds] = {}

                if exon_number_cds not in transcript_exon_no_exon_cds_dict[transcript_id_cds]:
                    transcript_exon_no_exon_cds_dict[transcript_id_cds][exon_number_cds] = {}
                    transcript_exon_no_exon_cds_dict[transcript_id_cds][exon_number_cds]["exon_id"] = {}
                    transcript_exon_no_exon_cds_dict[transcript_id_cds][exon_number_cds]["cds_location"] = []

                transcript_exon_no_exon_cds_dict[transcript_id_cds][exon_number_cds]["cds_location"].append(location)

    # Printing the output
    output_file_handle = open(new_file_path, "w")

    header = ["Gene_id", "Gene_location", "Gene_strand", "Total_transcript", "Transcript_id", "Total_Exons", "Exon_id",
              "Exon_location", "Exon_rank", "CDS_location"]
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
            # print gene, transcript, exon, total_transcript_count, exon_rank, total_exons_in_transcript
            cds_location_list = transcript_exon_no_exon_cds_dict[transcript][exon_rank]["cds_location"]

            if len(cds_location_list) > 0:
                cds_location = cds_location_list[0]
            else:
                cds_location = "nan"

            output = gene + "\t" + gene_location + "\t" + gene_strand + "\t" + str(
                total_transcript_count) + "\t" + transcript + "\t" + str(
                total_exons_in_transcript) + "\t" + exon + "\t" + exon_location + "\t" + str(
                exon_rank) + "\t" + cds_location + "\n"

            output_file_handle.write(output)

    output_file_handle.close()

    logging.info('parsed gtf saved at %s', new_file_path)

if __name__ == '__main__':
    if len(sys.argv) < 2 or not os.path.exists(sys.argv[1]):
        exit('missing config file!')
    params = json.load(open(sys.argv[1]))
    parse_gtf(params)
