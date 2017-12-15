import os
import sys
import json
import timeit
import logging

from helper.guide import Guide
from scores import score_utils
import download.main as downloads
import guides_info.main as guides_info
from helper.progress import ProgressLogger

OUTPUT_FOLDER = "scores"

def prepare_output_directory(params):
    base_path = os.path.join('output', params['ensembl_release'], params['species_name'])
    output_file_path = os.path.join(base_path, OUTPUT_FOLDER)
    try:
        if not os.path.exists(output_file_path):
            os.makedirs(output_file_path)
        return output_file_path
    except Exception as e:
        exit("\t ... Cannot create output directories %s" % output_file_path)


class Runner(object):
    def __init__(self, params, output_file_descript="_scores.txt"):
        self.output_file_path = prepare_output_directory(params)
        self.output_file_descript = output_file_descript

        begin = timeit.default_timer()
        transcript_cds_filepath = downloads.get_transcript_cds_filepath(params)
        logging.info('loading transcripts from %s', transcript_cds_filepath)
        self.transcript_cds_info_dict = score_utils.make_transcript_cds_info_dict(transcript_cds_filepath)

        protein_dir_path = os.path.join('input', params['ensembl_release'], params['species_name'], "proteins")
        logging.info('loading proteins from %s', protein_dir_path)
        self.protein_domain_info_dict = score_utils.make_protein_dict(protein_dir_path)

        self.snv_dict = {}
        if 'GVF_file' in params:
            compressed_variation_filepath = os.path.join('output', params['ensembl_release'], params['species_name'],
                                                         'Raw_data_files', os.path.basename(params['GVF_file']))
            variation_filepath = os.path.splitext(compressed_variation_filepath)[0]
            logging.info('loading variation from %s', variation_filepath)
            self.snv_dict = score_utils.make_var_dict(variation_filepath)

        now = timeit.default_timer()
        logging.info('time to load data %dsec' % (now - begin))
        self.start = timeit.default_timer()
        self.num_processed = 0
        self.progress = ProgressLogger(100)

    def _extract_gene_name(self, guide_filepath):
        return os.path.basename(guide_filepath).split("_")[0]

    def _get_output_filepath(self, guide_file, output_file_path):
        gene_name = self._extract_gene_name(guide_file)
        gene_output_file_name = gene_name + self.output_file_descript
        gene_output_file_name_path = os.path.join(output_file_path, gene_output_file_name)
        return gene_output_file_name_path

    def calculate_scores(self, input_file):
        guides = self._load_guides(input_file)
        if not guides:
            logging.warning('no guides for %s', input_file)

        scores = self._score_guides(guides)

        if scores:
            self._write(input_file_path, scores)

        self.progress.log()

    def _write(self, input_file_path, scores):
        output_file_name = self._get_output_filepath(input_file_path, self.output_file_path)
        header = ["Gene_id", "Exon_list", "Guide_with_ngg", "Guide_chr", "Guide_start", "Guide_stop", "Guide_strand", \
                  "SNP_score", "Domain_score", "Microhomology_score", "CDS_penalty_score", "Splicesite_penalty_score",
                  "Transcript_coverage_score", "Exon_ranking_score", "Total_score"]

        output_header = "\t".join(header) + "\n"
        with open(output_file_name, "w") as output_file_handle:
            output_file_handle.write(output_header)
            output_file_handle.writelines(scores)

    def _score_guides(self, guides):
        output_buffer = []
        for guide in guides:
            scores = self.score_guide(guide)
            if scores:
                output_buffer.append(scores)
        return output_buffer

    def _load_guides(self, guide_file):
        gene_name = self._extract_gene_name(input_file_path)
        logging.debug('gene %s', gene_name)

        with open(guide_file) as f:
            line = f.readline()  # throw away the header
            if not line.startswith('Gene_id'):
                logging.warning('seems like header is missing %s, %s', line, input_file_path)
            records = f.readlines()

        guides = [Guide.from_full_tsv(record) for record in records]
        guides = list(filter(None, guides))

        if len(records) != len(guides):
            logging.warning('not all guides loaded from %s', guide_file)


        return guides

    def score_guide(self, guide):
        exon_list_unstring = score_utils.parse_as_list(guide.guide_exons)  # ast.literal_eval(exon_list)
        transcript_list_unstring = score_utils.parse_as_list(guide.transcripts)  # ast.literal_eval(transcript_id_list)

        # in case of there is no proten-coding transcript or single non-coding exons
        if not exon_list_unstring or not transcript_list_unstring:
            return None

        try:
            prox_CDS_penalty = score_utils.calculate_proximity_to_CDS_for_transcript_list_modified(guide.transcripts,
                                                                                                   guide._compute_cutsite(),
                                                                                                   self.transcript_cds_info_dict)
            prox_splicesite_penalty = score_utils.calculate_proximity_splice_site_for_exon_rank_list_modified(
                guide.exon_rank_list, guide.cds_start_cutsite, guide.cds_stop_cutsite, guide.transcripts,
                self.transcript_cds_info_dict)
            exon_ranking_score = score_utils.calculate_exon_ranking_score_modified(guide.exon_rank_list)
            transcript_covered_score = score_utils.calculate_transcript_coverage_score_modified(guide.transcripts,
                                                                                                guide.transcript_count)
            domain_score = score_utils.calculate_score_protein_domains(guide._compute_cutsite(), guide.gene_id,
                                                                       self.protein_domain_info_dict)
            micrhomology_score = score_utils.calculate_microhomology_score(guide.microhomology_sequence)
            snp_score = score_utils.calculate_snp_score(guide.chromosome, guide._compute_cutsite(), self.snv_dict)

            scoring_list = list((prox_CDS_penalty, prox_splicesite_penalty,
                                 exon_ranking_score, transcript_covered_score, domain_score,
                                 micrhomology_score, snp_score))
            scoring_list_no_nan = [x for x in scoring_list if x != "nan"]
            scoring_list_no_nan_float = [float(x) for x in scoring_list_no_nan]
            total_score = sum(scoring_list_no_nan_float) / 10

            output = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\n".format(
                guide.gene_id, guide.guide_exons, guide.pamseq, guide.chromosome, guide.start, guide.end, guide.strand,
                str(snp_score), str(domain_score), str(micrhomology_score), str(prox_CDS_penalty),
                str(prox_splicesite_penalty), str(transcript_covered_score), str(exon_ranking_score),
                str(total_score))

            return output

        except Exception as e:
            logging.exception('failed to score %s %s, %s', guide.gene_id, guide, e)
            return None


def find_guide_files(params):
    base_path = os.path.join('output', params['ensembl_release'], params['species_name'])
    guides_info_dir = os.path.join(base_path, guides_info.OUTPUT_FOLDER)
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

    guide_files = find_guide_files(params)

    if not guide_files:
        logging.warning('no input files')
        exit(0)

    runner = Runner(params)
    start = timeit.default_timer()
    for input_file_path in guide_files:
        runner.calculate_scores(input_file_path)

    stop = timeit.default_timer()
    logging.info('time to compute scores %dsec', stop - start)
