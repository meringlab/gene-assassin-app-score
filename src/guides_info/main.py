#!/usr/bin/python

import os
import sys
import json
import logging
import timeit

from domain.guide import Guide
from domain.progress import ProgressLogger
from guides_info.exon_dict import ExonsInfo
import guides_info.chromosome_sequence_dict as chromosome_sequence_dict
import guides_info.guide_utils as guide_utils
import download.main as downloads


def find_target_exons(guide, gene, exon_dict):
    '''
      Searching for exons that this guide cuts.

      It's possible to have overlapping genes, and also overlapping transcripts where the same exon is
      coding in one but non-coding in the other one!
    :param guide:
    :param gene:
    :param exon_dict:
    :return:
    '''
    guide_exons = []
    if gene in exon_dict.gene_exons:
        for exon in exon_dict.gene_exons[gene]:  # a linear search is ok, only a handful of exons
            exon_start = exon_dict.exon_gene_info[exon]['start']
            exon_stop = exon_dict.exon_gene_info[exon]['stop']
            # if exon_start <= exon_stop:
            if guide.is_cutsite_within(exon_start, exon_stop):
                guide_exons.append(exon)

    return sorted(guide_exons)

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
                # this happened for human (and mouse), chromosome CHR_HSCHR6_MHC_QBL_CTG1 (in fasta but not in gtf)
                logging.warning("unknown chromosome %s, for %s", g.chromosome, line)
                continue

            try:
                guide_exons = find_target_exons(g, gene_name, exon_dict)
                logging.debug('gene %s, guide %s, guide_exon_list %s', gene_name, g.seq, guide_exons)

                cds_start_stop_cut_site = guide_utils.cutsite_distance_fom_exon_cds(guide_exons, g.cutsite, exon_dict)
                exon_genomic_features = guide_utils.exon_features_from_gtf(guide_exons, exon_dict)
                microhomology_sequence = "nan"
                exon_biotype_list = "nan"

                if guide_exons:
                    microhomology_sequence = g.seq_for_microhomology_scoring(sequence_dict)
                    exon_biotype_list = [exon_dict.exon_id_biotype[x]["biotype"] for x in guide_exons]

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
    if not os.path.exists(guide_file_info_directory):
        logging.info('creating output directory %s', guide_file_info_directory)
        os.makedirs(guide_file_info_directory)

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
