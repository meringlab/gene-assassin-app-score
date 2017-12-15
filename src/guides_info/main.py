#!/usr/bin/python

import os
import sys
import json
import logging
import timeit

from helper.guide import Guide
from helper.progress import ProgressLogger
from guides_info.exon_dict import ExonsInfo
import guides_info.chromosome_sequence_dict as chromosome_sequence_dict
import download.main as downloads

OUTPUT_FOLDER = "infofiles"


def make_guide_file_with_info(guide_file, output_file_path, sequence_dict, exon_dict):
    logging.debug('processing %s', guide_file)

    gene_name = os.path.basename(guide_file).split(".")[0]
    logging.debug('gene: %s', gene_name)

    if gene_name not in exon_dict.gene_exons:  # NMD genes for example
        logging.warning('not present in the parsed GTF %s', gene_name)

    with open(guide_file) as f:
        records = f.readlines()

    guides = [Guide.load_guide(gene_name, record, exon_dict, sequence_dict) for record in records]
    guides = list(filter(None, guides))
    if len(records) != len(guides):
        logging.warning('not all guides loaded from %s', guide_file)

    if guides:
        gene_output_file_path = os.path.join(output_file_path, gene_name + "_guides_info.txt")
        with open(gene_output_file_path, "w") as f:
            f.write(Guide.full_tsv_header() + "\n")
            f.write('\n'.join([g.to_tsv() for g in guides]))
    else:
        logging.info('no guides for this gene %s', gene_name)


def get_output_directory(params):
    ensembl_relase = params['ensembl_release']
    species = params['species_name']
    base_path = os.path.join('output', ensembl_relase, species)
    guide_file_info_directory = os.path.join(base_path, OUTPUT_FOLDER)
    if not os.path.exists(guide_file_info_directory):
        logging.info('creating output directory %s', guide_file_info_directory)
        os.makedirs(guide_file_info_directory)

    return guide_file_info_directory


def generate_guides(params, sequence_dict, output_directory, exon_dict):
    guides_directory = os.path.join('input', params['ensembl_release'], params['species_name'], 'guides')
    logging.info('generating guides info for %s', guides_directory)

    progress = ProgressLogger(1000)
    # weird cases for testing:
    # for guide_file_path in ['ENSDARG00000100456.guides.txt']: # has single nucleotide exon
    # for guide_file_path in ['ENSG00000026036.guides.txt']: # NMD gene
    for guide_file_path in sorted(os.listdir(guides_directory)):
        inputfile = os.path.join(guides_directory, guide_file_path)
        make_guide_file_with_info(inputfile, output_directory, sequence_dict, exon_dict)
        progress.log()


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
    return sequence_dict


if __name__ == "__main__":
    if len(sys.argv) < 2 or not os.path.exists(sys.argv[1]):
        exit('missing config file!')

    logging.basicConfig(filename=None, level=getattr(logging, 'DEBUG', None),
                        format='%(asctime)s %(levelname)s %(funcName)s %(message)s')
    params = json.load(open(sys.argv[1]))
    logging.info("computing info files, parameters: %s", params)

    output_directory = get_output_directory(params)
    exon_dict = load_exons_info(params)
    sequence_dict = load_sequence_dict(params)

    generate_guides(params, sequence_dict, output_directory, exon_dict)

##################### Notes

###### Exon list empty means no protein-coding transcript
####### Transcript list empty means two things : No protein-coding transcript and also that there is no coding exons so it is empty
###### Dist cds start-stop empty: no protein-coding transcript; Dist cds start-stop "nan" : no coding exons

####### Total score becomes "-20": when Exon list empty and transcript list is empty so in cases where there no coding transcript and also there is only a single non-coding exon for a transcript.
#
# import ast
# def test_for_overlapping_exon_guide(exon_list, exon_dict):
#     transcript_list_all = []
#
#     for exon in exon_list:
#         exon_feautre = extract_exon_features_from_gtf(exon, exon_dict)
#         transcript_list_all = transcript_list_all + ast.literal_eval(exon_feautre[1])
#
#     # find common transcripts in transcript_list_all
#     len_transcript_list_all = len(transcript_list_all)
#     len_unique_transcript_list_all = len(list(set(transcript_list_all)))
#
#     if len_unique_transcript_list_all == len_unique_transcript_list_all:
#         outcome = "no_exon_overlap"
#     else:
#         outcome = "exon_overlap"
#
#     return (outcome)
