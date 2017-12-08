import os
import sys
import json


def load_transcript_dicts(parsed_gtf_file_path):
    transcript_exon_info = {}
    transcript_exon_cds_length = {}
    input_file_handle = open(parsed_gtf_file_path)
    for line in input_file_handle:
        l = line.strip().split("\t")
        if l[0] != "Gene_id":
            gene_id = l[0];
            gene_strand = l[2];
            transcript_id = l[4];
            total_exons = l[5]
            exon_id = l[6];
            exon_cds = l[-1];
            exon_rank = int(l[-2])

            # print gene_id, gene_strand, transcript_id, total_exons, exon_id, exon_cds, exon_rank

            if exon_cds != "nan":

                exon_cds_start = int(exon_cds.split(":")[1].split("-")[0])
                exon_cds_stop = int(exon_cds.split(":")[1].split("-")[1])
                exon_chr = exon_cds.split(":")[0]

                if gene_strand == "+":
                    exon_cds_start_1 = exon_cds_start
                    exon_cds_stop_1 = exon_cds_stop

                if gene_strand == "-":
                    exon_cds_start_1 = exon_cds_stop
                    exon_cds_stop_1 = exon_cds_start

                exon_cds_length = abs(exon_cds_start - exon_cds_stop) + 1

                if transcript_id not in transcript_exon_info:
                    transcript_exon_info[transcript_id] = {}

                transcript_exon_info[transcript_id][exon_rank] = [exon_cds_start_1, exon_cds_stop_1, exon_chr,
                                                                  gene_strand]

                if transcript_id not in transcript_exon_cds_length:
                    transcript_exon_cds_length[transcript_id] = []

                transcript_exon_cds_length[transcript_id].append(exon_cds_length)
    return (transcript_exon_info,transcript_exon_cds_length)



if __name__ == "__main__":

    if len(sys.argv) < 2 or not os.path.exists(sys.argv[1]):
        exit('missing config file!')
    params = json.load(open(sys.argv[1]))

    gtf_file_name = os.path.basename(params['base_url_gtf'])
    parset_gtf_file_name = 'parsed_' + gtf_file_name.replace('.gz', '.txt')
    # "output/v85_2017/danio_rerio/v_85/Proesssed_data_files/parsed_Danio_rerio.GRCz10.85.gtf.txt"
    parsed_gtf_file_path = os.path.join('output', params['ensembl_release'], params['species_name'], 'Processed_data_files',
                                        parset_gtf_file_name)
    if not os.path.exists(parsed_gtf_file_path):
        exit("\t ... parsed gtf file path does not exist %s."% parsed_gtf_file_path)


    transcript_exon_info,transcript_exon_cds_length = load_transcript_dicts(parsed_gtf_file_path)
    #for transcript in transcript_exon_cds_length:
    #    print transcript, transcript_exon_cds_length[transcript], sum(transcript_exon_cds_length[transcript])

    base_output_file_path = os.path.dirname(parsed_gtf_file_path)
    output_file_name = "transcript_cds_info_" + os.path.basename(parsed_gtf_file_path)
    output_file_path = os.path.join(base_output_file_path, output_file_name)

    output_file_handle = open(output_file_path, "w")


    header = ["Transcript_id", "First_cds_exon_rank_in_transcript", "Last_cds_exon_rank_in_transcript", "Transcript_chr", "Transcript_cds_start", "Transcript_cds_stop","Transcript_strand", "Transcipt_cds_length"]
    output = "\t".join(header) + "\n"
    output_file_handle.write(output)

    for transcript in transcript_exon_info:

        #if transcript == 'ENSDART00000003544':

        sorted_dict = sorted(transcript_exon_info[transcript].items())

        transcript_first_cds_exon_rank_info = sorted_dict[0]
        transcript_last_cds_exon_rank_info = sorted_dict[-1]

        transcript_first_cds_exon_rank = transcript_first_cds_exon_rank_info[0]
        transcript_last_cds_exon_rank = transcript_last_cds_exon_rank_info[0]

        transcript_cds_start = transcript_first_cds_exon_rank_info[1][0]
        transcript_cds_stop = transcript_last_cds_exon_rank_info[1][1]

        transcript_strand = transcript_first_cds_exon_rank_info[1][-1]
        transcript_chr = transcript_first_cds_exon_rank_info[1][-2]

        #print transcript, transcript_chr, transcript_strand, transcript_first_cds_exon_rank, transcript_last_cds_exon_rank, transcript_cds_start, transcript_cds_stop,

        output= transcript +"\t" + str(transcript_first_cds_exon_rank) +"\t" + str(transcript_last_cds_exon_rank)  +"\t" + \
        str(transcript_chr)+"\t" + str(transcript_cds_start)  +"\t" + str(transcript_cds_stop)  +"\t" + transcript_strand  +"\t" +\
        str(sum(transcript_exon_cds_length[transcript])) + "\n"

        output_file_handle.write(output)


    output_file_handle.close()
