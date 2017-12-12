#!/usr/bin/python

import os
import sys
import json
import logging
import timeit
from Bio.Seq import Seq

from guides_info.exon_dict import ExonsInfo
import guides_info.chromosome_sequence_dict as chromosome_sequence_dict
import guides_info.guide_utils as guide_utils
import download.main as downloads


class Region(object):
    def __init__(self, chromosome, start, end):
        start = int(start)
        end = int(end)
        if start > end:
            raise Exception('invalid location: %d-%d' % (start, end))
        self.start = start
        self.end = end
        self.chromosome = chromosome

    @classmethod
    def parse(cls, str):
        chr = str[:str.index(':')]
        # start, stop = list(map(int, str[str.index(':') + 1:].split('-')))
        start, end = str[str.index(':') + 1:].split('-')
        return cls(chr, start, end)

    def __str__(self):
        return '%s:%d-%d' % (self.chromosome, self.start, self.end)

    def __repr__(self):
        return str(self)


class Guide(object):
    _seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

    def __init__(self, tsv):
        self._load_from_tsv(tsv)

    def _load_from_tsv(self, line):
        l = line.strip().split('\t')
        self.chromosome = l[0]
        self.start = int(l[1])
        self.end = int(l[2])
        self.strand = l[3]
        self.seq = l[4]
        self.offtarget_profile = l[10]
        self.uniqueness = self.offtarget_profile.split(",")[1]
        self.cutsite = self._compute_cutsite()

    def __str__(self):
        return '%s %s:%d-$d' % (self.strand, self.seq, self.start, self.end)

    def is_on_forward_strand(self):
        return self.strand == "+" or self.strand == "1"

    def _compute_cutsite(self):
        if self.is_on_forward_strand():
            return self.end - 2
        return self.start + 2

    # def set_pam(self, sequence):
    #     self.ngg = sequence
    #
    # def get_pam(self):
    #     return self.ngg
    #
    # pam = property(get_pam, set_pam)

    def is_cutsite_within(self, start, stop):
        return start <= self.cutsite <= stop

    def _reverse_complement(self, seq):
        # return "".join([Guide._seq_dict[base] for base in reversed(seq)])
        return str(Seq(seq).reverse_complement())

    def seq_with_pam(self, chromosome_sequence_dict):
        chr_seq = chromosome_sequence_dict[self.chromosome]

        if self.is_on_forward_strand():
            return chr_seq[self.start - 1: self.end + 3]
        seq = chr_seq[self.start - 3 - 1:self.end]
        seq = self._reverse_complement(seq)

        ##### test
        # seq_ngg_regex = "[ATGC]{20}[ATGC]{1}GG"
        # search_object_ngg = re.search(seq_ngg_regex, seq)
        # if search_object_ngg:
        #     ngg_outcome = seq
        # else:
        #     ngg_outcome = "nan"

        return seq

    def seq_for_microhomology_scoring(self, chromosome_sequence_dict):
        sequence_start = self.start - 13
        sequence_end = self.end + 27

        chr_seq = chromosome_sequence_dict[self.chromosome]
        if self.is_on_forward_strand():
            return chr_seq[sequence_start - 1:sequence_end]
        return self._reverse_complement(chr_seq[sequence_start - 1:sequence_end])


def make_guide_file_with_info(guide_file, output_file_path, sequence_dict, exon_dict):
    logging.debug('processing %s', guide_file)

    gene_name = os.path.basename(guide_file).split(".")[0]
    logging.debug('gene: %s', gene_name)

    if gene_name not in exon_dict.gene_exons:  # NMD genes for example
        logging.warning('not present in the parsed gtf %s', gene_name)

    buffer = []
    with open(guide_file) as guide_input_file:
        for line in guide_input_file:
            g = Guide(line)

            if g.chromosome not in exon_dict.chromosomes:
                # this happened for human, chromosome CHR_HSCHR6_MHC_QBL_CTG1 (in fasta but not in gtf)
                logging.warning("unknown chromosome %s, for %s", g.chromosome, line)
                continue

            # Searching for exons that this guide cuts
            # it's possible to have overlapping genes, and also overlapping transcripts where the same exon is
            # coding in one but non-coding in the other one!
            guide_exons = []
            if gene_name in exon_dict.gene_exons:
                for exon in exon_dict.gene_exons[gene_name]:  # a linear search is ok, only a handful of exons
                    exon_start = exon_dict.exon_gene_info[exon]['start']
                    exon_stop = exon_dict.exon_gene_info[exon]['stop']
                    # if exon_start <= exon_stop:
                    if g.is_cutsite_within(exon_start, exon_stop):
                        guide_exons.append(exon)

            guide_exons = sorted(guide_exons)
            logging.debug('gene %s, guide %s, guide_exon_list %s', gene_name, g.seq, guide_exons)

            try:
                # Extracting feature
                cds_start_stop_cut_site = guide_utils.calculate_cutsite18_dist_fom_exon_cds_start_stop_for_exon_list_modified(
                    guide_exons, str(g.cutsite), exon_dict)
                exon_genomic_features = guide_utils.extract_exon_features_from_gtf_for_exonlist_modified(
                    guide_exons, exon_dict)

                if guide_exons:
                    microhomology_sequence = g.seq_for_microhomology_scoring(sequence_dict)
                    exon_biotype_list = [exon_dict.exon_id_biotype[x]["biotype"] for x in guide_exons]
                else:
                    microhomology_sequence = "nan"
                    exon_biotype_list = "nan"

                output = "\t".join(
                    [gene_name, (g.seq_with_pam(sequence_dict)), g.seq, g.chromosome, str(g.start), str(g.end),
                     g.strand, g.uniqueness, str(guide_exons), str(exon_biotype_list), str(g.cutsite),
                     "\t".join(cds_start_stop_cut_site), "\t".join(exon_genomic_features), microhomology_sequence])
                buffer.append(output)

            except Exception as e:
                logging.exception('failed to collect info on %s - %s, error: %s', gene_name, line, e)

    if buffer:
        gene_output_file_name = gene_name + "_guides_info.txt"
        gene_output_file_name_path = os.path.join(output_file_path, gene_output_file_name)
        with open(gene_output_file_name_path, "w") as gene_output_file_handle:
            header = ["Gene_id", "Guide_seq_with_ngg", "Guide_seq", "Guide_chr_no", "Guide_start", "Guide_stop",
                      "Guide_strand",
                      "Guide_uniqness",
                      "Guide_exon_list", "Exon_biotype_list",
                      "cutsite18", "dist_cds_start_cutsite", "dist_cds_stop_cutsite",
                      "Total_transcript_count", "Transcript_list", "Total_exons_in_transcript_list",
                      "Exon_rank_in_transcript_list",
                      "microhomology_seq"]
            header_output = "\t".join(header) + "\n"
            gene_output_file_handle.write(header_output)
            gene_output_file_handle.write('\n'.join(buffer))
    else:
        logging.info('no guides for this gene %s', gene_name)


def prepareOutputDirectory(params):
    ensembl_relase = params['ensembl_release']
    species = params['species_name']
    base_path = os.path.join('output', ensembl_relase, species)
    dir_to_be_created = "Guide_files_with_information"
    guide_file_info_directory = os.path.join(base_path, dir_to_be_created)
    try:
        if not os.path.exists(guide_file_info_directory):
            logging.info('creating output directory %s', guide_file_info_directory)
            os.makedirs(guide_file_info_directory)
    except Exception as e:
        logging.exception('Cannot create output directories %s', e)
        exit("\t ... Cannot create output directories %s" % guide_file_info_directory)

    return guide_file_info_directory


def generate_guides(params, sequence_dict, output_directory, exon_dict):
    logging.info('generating guides info')
    guides_directory = os.path.join('input', params['ensembl_release'], params['species_name'], 'guides')
    logging.info('reading guides from %s', guides_directory)

    start = timeit.default_timer()
    progress = ProgressLogger(1000)
    # weird cases for testing:
    # for guide_file_path in ['ENSDARG00000100456.guides.txt']: # has single nucleotide exon
    # for guide_file_path in ['ENSG00000026036.guides.txt']: # NMD gene
    for guide_file_path in os.listdir(guides_directory):
        make_guide_file_with_info(os.path.join(guides_directory, guide_file_path), output_directory, sequence_dict,
                                  exon_dict)
        progress.log()
    stop = timeit.default_timer()
    logging.info('time to generate guides %dsec', stop - start)


def load_exons_info(params):
    parsed_file_path = downloads.get_parsed_gtf_filepath(params)
    logging.info('loading exons from %s', parsed_file_path)
    return ExonsInfo(parsed_file_path)


def load_sequence_dict(params):
    # chromosome_data_path = "output/v85_2017/danio_rerio/v_85/Raw_data_files/Danio_rerio.GRCz10.dna.toplevel.fa.gz"
    compressed_fasta_file_name = os.path.basename(params['DNA_top_level_fa'])
    fasta_file_name = os.path.splitext(compressed_fasta_file_name)[0]
    fasta_file_path = os.path.join('output', params['ensembl_release'], params['species_name'], 'Raw_data_files',
                                   fasta_file_name)
    logging.info('loading fasta from %s', fasta_file_name)
    sequence_dict = chromosome_sequence_dict.load_chromosome_sequence_dict_from_fasta(fasta_file_path)
    # fasta_file_path = '../ga-pipeline/input/human/fasta/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz'
    # logging.info('loading fasta from %s', fasta_file_name)
    # sequence_dict = chromosome_sequence_dict.load_chromosome_sequence_dict_list(fasta_file_path)
    return sequence_dict


class ProgressLogger(object):
    def __init__(self, report):
        self.counter = 0
        self.start = timeit.default_timer()
        self.report = report

    def log(self):
        self.counter += 1
        if self.counter % self.report == 0:
            stop = timeit.default_timer()
            logging.info('%d processed in %dsec', self.counter, stop - self.start)

if __name__ == "__main__":
    if len(sys.argv) < 2 or not os.path.exists(sys.argv[1]):
        exit('missing config file!')

    logging.basicConfig(filename=None, level=getattr(logging, 'DEBUG', None),
                        format='%(asctime)s %(levelname)s %(funcName)s %(message)s')
    params = json.load(open(sys.argv[1]))
    logging.info("computing info files, parameters: %s", params)

    output_directory = prepareOutputDirectory(params)
    exon_dict = load_exons_info(params)
    sequence_dict = load_sequence_dict(params)

    generate_guides(params, sequence_dict, output_directory, exon_dict)
