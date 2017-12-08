import os
import sys
import json
import timeit
import logging
from scores import score_utils
import download.main as downloads


def get_output_filepath(guide_file, output_file_path, output_file_descript):
    gene_name = os.path.basename(guide_file).split("_")[0]
    gene_output_file_name = gene_name + output_file_descript  # like "_guide_info.txt"

    if os.path.exists(output_file_path):
        gene_output_file_name_path = os.path.join(output_file_path, gene_output_file_name)
    else:
        exit("Output file path does not exist %s" % output_file_path)

    return gene_output_file_name_path


def calculate_score_guide_main(input_file_path, output_file_path, output_file_descript, transcript_cds_info_dict,
                               protein_domain_info_dict, snv_dict):
    output_file_name = get_output_filepath(input_file_path, output_file_path, output_file_descript)

    gene_name = os.path.basename(input_file_path).split("_")[0]
    logging.debug('gene %s', gene_name)
    input_file_handle = open(input_file_path)

    line = input_file_handle.readline() # throw away the header
    if not line.startswith('Gene_id'):
        logging.warning('seems like header is missing %s', line)

    output_buffer = []
    for line in input_file_handle:
        logging.debug('line: %s', line)
        l = line.strip().split('\t')
        if len(l) < 2:
            continue
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
        exon_list_unstring = score_utils.parse_as_list(exon_list) # ast.literal_eval(exon_list)
        transcript_list_unstring = score_utils.parse_as_list(transcript_id_list) # ast.literal_eval(transcript_id_list)

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

                guide_transcripts_prox_CDS_penalty = score_utils.calculate_proximity_to_CDS_for_transcript_list_modified(transcript_id_list, cutsite18, transcript_cds_info_dict)

                guide_exons_prox_splicesite_penalty = score_utils.calculate_proximity_splice_site_for_exon_rank_list_modified(exon_rank_list, dist_cutsite_exon_cds_start_list_input=dist_cutsite_cds_start,                    dist_cutsite_exon_cds_stop_list_input=dist_cutsite_cds_stop,                    transcript_list_input=transcript_id_list, transcript_cds_info_dict=transcript_cds_info_dict)

                guide_exon_ranking_score = score_utils.calculate_exon_ranking_score_modified(exon_ranks=exon_rank_list)

                guide_transcript_covered_score = score_utils.calculate_transcript_coverage_score_modified(transcript_id_list, transcript_count)

                domain_score = score_utils.calculate_score_protein_domains(cutsite18, gene_id, protein_domain_info_dict)

                calculated_micrhomology_score = score_utils.calculate_microhomology_score(seq=microhomogy_guide)

                snp_score = score_utils.calculate_snp_score(guide_chr, cutsite18, snv_dict)

                # Total scores
                scoring_list = list((guide_transcripts_prox_CDS_penalty, guide_exons_prox_splicesite_penalty,
                                     guide_exon_ranking_score, guide_transcript_covered_score, domain_score,
                                     calculated_micrhomology_score, snp_score))
                scoring_list_no_nan = [x for x in scoring_list if x != "nan"]
                scoring_list_no_nan_float = [float(x) for x in scoring_list_no_nan]
                total_score = sum(scoring_list_no_nan_float) / 10

                output = gene_id + "\t" + exon_list + "\t" + guide_with_ngg + "\t" + guide_chr + "\t" + guide_start + "\t" + guide_stop + "\t" + guide_strand + "\t" + \
                         str(snp_score) + "\t" + str(domain_score) + "\t" + str(calculated_micrhomology_score) + \
                         "\t" + str(guide_transcripts_prox_CDS_penalty) + "\t" + str(
                    guide_exons_prox_splicesite_penalty) + "\t" + str(guide_transcript_covered_score) + "\t" + str(
                    guide_exon_ranking_score) + \
                         "\t" + str(total_score) + "\n"

                output_buffer.append(output)

            except Exception as e:
                logging.exception('failed to score %s %s, %s', gene_name, guide_seq, e)


    if output_buffer:
        header = ["Gene_id", "Exon_list", "Guide_with_ngg", "Guide_chr", "Guide_start", "Guide_stop", "Guide_strand", \
                  "SNP_score", "Domain_score", "Microhomology_score", "CDS_penalty_score", "Splicesite_penalty_score",
                  "Transcript_coverage_score", "Exon_ranking_score", "Total_score"]

        output_header = "\t".join(header) + "\n"
        with open(output_file_name, "w") as output_file_handle:
            output_file_handle.write(output_header)
            output_file_handle.writelines(output_buffer)
    else:
        logging.warning('no scores for %s', gene_name)

    input_file_handle.close()

if __name__ == "__main__":
    if len(sys.argv) < 2 or not os.path.exists(sys.argv[1]):
        exit('missing config file!')
    logging.basicConfig(filename=None, level=getattr(logging, 'INFO', None),
                        format='%(asctime)s %(levelname)s %(funcName)s %(message)s')

    params = json.load(open(sys.argv[1]))
    logging.info("computing scores, parameters: %s", params)

    ensembl_relase = params['ensembl_release']
    species = params['species_name']

    base_path = os.path.join('output', ensembl_relase, species)
    dir_to_be_created = "Guide_files_with_scores"

    output_file_path = os.path.join(base_path, dir_to_be_created)
    try:
        if not os.path.exists(output_file_path):
            os.makedirs(output_file_path)
    except Exception as e:
        exit("\t ... Cannot create output directories %s" % output_file_path)

    output_file_descript = "_scores.txt"

    start = timeit.default_timer()

    transcript_cds_filepath = downloads.get_transcript_cds_filepath(params)
    logging.info('loading transcripts from %s', transcript_cds_filepath)
    transcript_cds_info_dict = score_utils.make_transcript_cds_info_dict(transcript_cds_filepath)

    protein_dir_path = os.path.join('input', ensembl_relase, species, "proteins")
    logging.info('loading proteins from %s', protein_dir_path)
    protein_domain_info_dict = score_utils.make_protein_dict(protein_dir_path)

    snv_dict = {}
    if 'GVF_file' in params:
        compressed_variation_filepath = os.path.join('output', ensembl_relase, species, 'Raw_data_files',
                                                     os.path.basename(params['GVF_file']))
        variation_filepath = os.path.splitext(compressed_variation_filepath)[0]
        logging.info('loading variation from %s', compressed_variation_filepath)
        snv_dict = score_utils.make_var_dict(variation_filepath)

    stop = timeit.default_timer()
    logging.info('time to prepare for computation %dsec' % (stop - start))

    start = timeit.default_timer()

    num_processed = 0
    guides_info_dir = os.path.join(base_path, 'Guide_files_with_information/')
    logging.info('guide info directory: %s', guides_info_dir)
    # for input_file in ['ENSG00000005001_guides_info.txt']:
    # for input_file in ['ENSDARG00000000966_guides_info.txt']:
    # for input_file in ['ENSDARG00000000189_guides_info.txt']:
    for input_file in sorted(os.listdir(guides_info_dir)):
        input_file_path = os.path.join(guides_info_dir, input_file)
        logging.info('scoring %s', input_file_path)

        calculate_score_guide_main(input_file_path, output_file_path, output_file_descript, transcript_cds_info_dict,
                                   protein_domain_info_dict, snv_dict)
        num_processed += 1
        if num_processed % 10 == 0:
            stop = timeit.default_timer()
            logging.info('%d processed in %dsec', num_processed, stop - start)

    stop = timeit.default_timer()
    logging.info('time to compute scores %dsec', stop - start)
