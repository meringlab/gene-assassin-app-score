import os
import sys
import json
import timeit
import ast

import all_scoring_function as fn_scoring
from Universal import universal_function as fn_universal


######  Main function to calculate the scores
def calculate_score_guide_main(input_file_path, output_file_path, output_file_descript, transcript_cds_info_dict,
                               protein_domain_info_dict, snv_dict):
    ###### Handling files
    output_file_name = \
    fn_universal.making_output_log_file_names(input_file_path, output_file_path, output_file_descript)[0]
    log_file_name = fn_universal.making_output_log_file_names(input_file_path, output_file_path, output_file_descript)[
        1]

    input_file_handle = open(input_file_path)
    output_file_handle = open(output_file_name, "w")
    output_log_file_handle = open(log_file_name, "a")

    ####### Gene name

    gene_name = os.path.basename(input_file_path).split("_")[0]

    ###### Printing the header
    header = ["Gene_id", "Exon_list", "Guide_with_ngg", "Guide_chr", "Guide_start", "Guide_stop", "Guide_strand", \
              "SNP_score", "Domain_score", "Microhomology_score", "CDS_penalty_score", "Splicesite_penalty_score",
              "Transcript_coverage_score", "Exon_ranking_score", "Total_score"]

    output_header = "\t".join(header) + "\n"
    output_file_handle.write(output_header)

    ######## Processing the file

    for line in input_file_handle:
        l = line.strip("\n").split("\t")
        if l[0] != "Gene_id":

            gene_id = l[0];
            guide_with_ngg = l[1];
            guide_seq = l[2];
            guide_chr = l[3];
            guide_start = l[4];
            guide_stop = l[5];
            guide_strand = l[6];
            guide_uniqueness = l[7];
            exon_list = l[8];
            exon_biotype = l[9];
            cutsite18 = l[10];
            dist_cutsite_cds_start = l[11];
            dist_cutsite_cds_stop = l[12]
            transcript_count = l[13];
            transcript_id_list = l[14];
            exon_rank_list = l[-2]
            microhomogy_guide = l[-1]

            ######### Making the string as lisT
            exon_list_unstring = ast.literal_eval(exon_list)
            transcript_list_unstring = ast.literal_eval(transcript_id_list)

            ###### in case of there is no proten-coding transcript or single non-coding exons
            if len(exon_list_unstring) == 0 or len(transcript_list_unstring) == 0:
                total_score = -20
                snp_score = "nan";
                domain_score = "nan";
                calculated_micrhomology_score = "nan";
                guide_transcripts_prox_CDS_penalty = "nan"
                guide_exons_prox_splicesite_penalty = "nan";
                guide_transcript_covered_score = "nan";
                guide_exon_ranking_score = "nan"
            else:

                try:
                    #########  Genomic Context Score

                    guide_transcripts_prox_CDS_penalty = fn_scoring.calculate_proximity_to_CDS_for_transcript_list_modified(
                        transcript_id_list, cutsite18, transcript_cds_info_dict)

                    guide_exons_prox_splicesite_penalty = fn_scoring.calculate_proximity_splice_site_for_exon_rank_list_modified(
                        exon_rank_list, dist_cutsite_exon_cds_start_list_input=dist_cutsite_cds_start, \
                        dist_cutsite_exon_cds_stop_list_input=dist_cutsite_cds_stop,
                        transcript_list_input=transcript_id_list, transcript_cds_info_dict=transcript_cds_info_dict)

                    guide_exon_ranking_score = fn_scoring.calculate_exon_ranking_score_modified(
                        exon_ranks=exon_rank_list)

                    guide_transcript_covered_score = fn_scoring.calculate_transcript_coverage_score_modified(
                        transcript_id_list, transcript_count)

                    ################################### Protein Domain

                    domain_score = fn_scoring.calculate_score_protein_domains(cutsite18, gene_id,
                                                                              protein_domain_info_dict)

                    ############### Microhomology
                    calculated_micrhomology_score = fn_scoring.calculate_microhomology_score(seq=microhomogy_guide)

                    ######## Snp scoring

                    snp_score = fn_scoring.calculate_snp_score(guide_chr, cutsite18, snv_dict)

                    ######### Total scores

                    scoring_list = list((guide_transcripts_prox_CDS_penalty, guide_exons_prox_splicesite_penalty,
                                         guide_exon_ranking_score, guide_transcript_covered_score, domain_score,
                                         calculated_micrhomology_score, snp_score))
                    scoring_list_no_nan = [x for x in scoring_list if x != "nan"]
                    scoring_list_no_nan_float = [float(x) for x in scoring_list_no_nan]
                    total_score = sum(scoring_list_no_nan_float) / 10

                except Exception as e:
                    output_exception = gene_name + "\t" + str(e) + "\t" + guide_seq + "\n"
                    output_log_file_handle.write(output_exception)
                    total_score = "nan"

                    ############# Printing

            output = gene_id + "\t" + exon_list + "\t" + guide_with_ngg + "\t" + guide_chr + "\t" + guide_start + "\t" + guide_stop + "\t" + guide_strand + "\t" + \
                     str(snp_score) + "\t" + str(domain_score) + "\t" + str(calculated_micrhomology_score) + \
                     "\t" + str(guide_transcripts_prox_CDS_penalty) + "\t" + str(
                guide_exons_prox_splicesite_penalty) + "\t" + str(guide_transcript_covered_score) + "\t" + str(
                guide_exon_ranking_score) + \
                     "\t" + str(total_score) + "\n"

            output_file_handle.write(output)

    output_file_handle.close()
    output_log_file_handle.close()


if __name__ == "__main__":
    if len(sys.argv) < 2 or not os.path.exists(sys.argv[1]):
        exit('missing config file!')
    params = json.load(open(sys.argv[1]))

    ensembl_relase = params['ensembl_release']
    species = params['species_name']

    base_path = os.path.join('output', ensembl_relase, species)
    dir_to_be_created = "Guide_files_with_scores"

    guide_file_score_directory = os.path.join(base_path, dir_to_be_created)
    try:
        if not os.path.exists(guide_file_score_directory):
            os.makedirs(guide_file_score_directory)
    except Exception as e:
        exit("\t ... Cannot create output directories %s" % guide_file_score_directory)


    output_file_path =  guide_file_score_directory
    output_file_descript = "_scores.txt"

    #####Step 1
    ### Step 2
    #### Transcript cds info dict

    start = timeit.default_timer()

    transcript_cds_file_name = params['transcript_cds_file_name']
    transcript_cds_info_dict = fn_scoring.make_transcript_cds_info_dict (transcript_cds_file_name, base_path)


    ###### Protein dict
    name_protein_dir = "proteins"
    protein_dir_path = os.path.join('input', ensembl_relase, species)
    protein_domain_info_dict =  fn_scoring.make_protein_dict (name_protein_dir, protein_dir_path)


    ###### Variation dict
    variation_file_name = os.path.basename(params['GVF_file'])
    var_dir_path = os.path.join('output', ensembl_relase, species, 'Raw_data_files')
    snv_dict = fn_scoring.make_var_dict(variation_file_name,var_dir_path)
    stop = timeit.default_timer()
    print('time to prepare for computation %dsec' % (stop - start))

    start = timeit.default_timer()


    ################ callling the main function
    # input_file_path = sys.argv[1]
    num_processed = 0

    guides_info_dir = 'output/v85_2017/danio_rerio/v_85/Guide_files_with_information/v0/'
    for input_file in os.listdir(guides_info_dir):
        input_file_path = os.path.join(guides_info_dir, input_file)
        # input_file_path = os.path.join(guides_info_dir, 'ENSDARG00000079029_guides_info.txt')
        calculate_score_guide_main (input_file_path, output_file_path, output_file_descript, transcript_cds_info_dict, protein_domain_info_dict,snv_dict)
        num_processed += 1
        if num_processed % 10 == 0:
            stop = timeit.default_timer()
            print('%d processed in %dsec' % (num_processed, stop - start))

    stop = timeit.default_timer()
    print('time to compute scores %dsec' % (stop - start))

