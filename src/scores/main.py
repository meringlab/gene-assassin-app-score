import os
import sys
import json
import timeit
import logging
from scores import score_utils
import download.main as downloads

output_file_descript = "_scores.txt"


def _extract_gene_name(guide_filepath):
    return os.path.basename(guide_filepath).split("_")[0]


def get_output_filepath(guide_file, output_file_path):
    gene_name = _extract_gene_name(guide_file)
    gene_output_file_name = gene_name + output_file_descript

    if os.path.exists(output_file_path):
        gene_output_file_name_path = os.path.join(output_file_path, gene_output_file_name)
    else:
        exit("Output file path does not exist %s" % output_file_path)

    return gene_output_file_name_path


def calculate_score_guide_main(input_file_path, output_file_path, transcript_cds_info_dict,
                               protein_domain_info_dict, snv_dict):
    logging.debug('scoring %s', input_file_path)

    output_file_name = get_output_filepath(input_file_path, output_file_path)

    gene_name = _extract_gene_name(input_file_path)
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
        if exon_list_unstring and transcript_list_unstring:
            continue

        try:
            prox_CDS_penalty = score_utils.calculate_proximity_to_CDS_for_transcript_list_modified(transcript_id_list,
                                                                                                   cutsite18,
                                                                                                   transcript_cds_info_dict)
            prox_splicesite_penalty = score_utils.calculate_proximity_splice_site_for_exon_rank_list_modified(
                exon_rank_list, dist_cutsite_cds_start, dist_cutsite_cds_stop, transcript_id_list,
                transcript_cds_info_dict)
            exon_ranking_score = score_utils.calculate_exon_ranking_score_modified(exon_rank_list)
            transcript_covered_score = score_utils.calculate_transcript_coverage_score_modified(transcript_id_list,
                                                                                                transcript_count)
            domain_score = score_utils.calculate_score_protein_domains(cutsite18, gene_id, protein_domain_info_dict)
            micrhomology_score = score_utils.calculate_microhomology_score(seq=microhomogy_guide)
            snp_score = score_utils.calculate_snp_score(guide_chr, cutsite18, snv_dict)

            scoring_list = list((prox_CDS_penalty, prox_splicesite_penalty,
                                 exon_ranking_score, transcript_covered_score, domain_score,
                                 micrhomology_score, snp_score))
            scoring_list_no_nan = [x for x in scoring_list if x != "nan"]
            scoring_list_no_nan_float = [float(x) for x in scoring_list_no_nan]
            total_score = sum(scoring_list_no_nan_float) / 10

            output = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\n".format(
                gene_id, exon_list, guide_with_ngg, guide_chr, guide_start, guide_stop, guide_strand,
                str(snp_score), str(domain_score), str(micrhomology_score), str(prox_CDS_penalty),
                str(prox_splicesite_penalty), str(transcript_covered_score), str(exon_ranking_score),
                str(total_score))

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


def prepare_output_directory(params):
    base_path = os.path.join('output', params['ensembl_release'], params['species_name'])
    dir_to_be_created = "Guide_files_with_scores"
    output_file_path = os.path.join(base_path, dir_to_be_created)
    try:
        if not os.path.exists(output_file_path):
            os.makedirs(output_file_path)
        return output_file_path
    except Exception as e:
        exit("\t ... Cannot create output directories %s" % output_file_path)


class Runner(object):
    def __init__(self, params):
        self.output_file_path = prepare_output_directory(params)

        start = timeit.default_timer()
        transcript_cds_filepath = downloads.get_transcript_cds_filepath(params)
        logging.info('loading transcripts from %s', transcript_cds_filepath)
        self.transcript_cds_info_dict = score_utils.make_transcript_cds_info_dict(transcript_cds_filepath)

        protein_dir_path = os.path.join('input', params['ensembl_release'], params['species_name'], "proteins")
        logging.info('loading proteins from %s', protein_dir_path)
        self.protein_domain_info_dict = score_utils.make_protein_dict(protein_dir_path)

        self.snv_dict = {}
        if 'GVF_file' in params:
            compressed_variation_filepath = os.path.join('output', params['ensembl_release'], params['species_name'],
                                                         'Raw_data_files',
                                                         os.path.basename(params['GVF_file']))
            variation_filepath = os.path.splitext(compressed_variation_filepath)[0]
            logging.info('loading variation from %s', variation_filepath)
            self.snv_dict = score_utils.make_var_dict(variation_filepath)
        logging.info('time to load data %dsec' % (stop - start))
        self.start = timeit.default_timer()
        self.num_processed = 0

    def calculate_scores(self, input_file):
        calculate_score_guide_main(input_file_path, self.output_file_path, self.transcript_cds_info_dict,
                                   self.protein_domain_info_dict, self.snv_dict)
        self.num_processed += 1
        if self.num_processed % 100 == 0:
            now = timeit.default_timer()
            logging.info('%d processed in %dsec', self.num_processed, now - self.start)


def find_guide_files(params):
    base_path = os.path.join('output', params['ensembl_release'], params['species_name'])
    guides_info_dir = os.path.join(base_path, 'Guide_files_with_information/')
    logging.info('guide info directory: %s', guides_info_dir)
    guides = [os.path.join(guides_info_dir, f) for f in sorted(os.listdir(guides_info_dir))]
    return guides


if __name__ == "__main__":
    if len(sys.argv) < 2 or not os.path.exists(sys.argv[1]):
        exit('missing config file!')
    logging.basicConfig(filename=None, level=getattr(logging, 'INFO', None),
                        format='%(asctime)s %(levelname)s %(funcName)s %(message)s')

    params = json.load(open(sys.argv[1]))
    logging.info("computing scores, parameters: %s", params)

    guide_files = find_guide_files()

    if not guide_files:
        logging.warning('no input files')
        exit(0)

    runner = Runner(params)
    start = timeit.default_timer()
    for input_file_path in guide_files:
        runner.calculate_scores(input_file_path)

    stop = timeit.default_timer()
    logging.info('time to compute scores %dsec', stop - start)
