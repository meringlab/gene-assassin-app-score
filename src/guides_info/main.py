#!/usr/bin/python

import os
import sys
import json
import logging
import timeit

from guides_info.exon_dict import ExonsInfo
import guides_info.chromosome_sequence_dict as chromosome_sequence_dict
import guides_info.guide_utils as guide_utils
import download.main as downloads


def making_guide_file_with_info(guide_file, output_file_path, sequence_dict, exon_dict):
    gene_name = os.path.basename(guide_file).split(".")[0]
    logging.debug('gene: %s', gene_name)

    buffer = []
    with open(guide_file) as guide_input_file:
        for line in guide_input_file:

            l = line.strip().split('\t')
            guide_chr = l[0]
            guide_start = l[1]
            guide_stop = l[2]
            guide_strand = l[3]
            guide_seq = l[4]
            guide_offtarget_profile = l[10]

            if guide_chr not in exon_dict.chromosomes:
                # this happened for human, chromosome CHR_HSCHR6_MHC_QBL_CTG1 (in fasta but not in gtf)
                logging.warning("unknown chromosome %s, for %s", guide_chr, line)
                continue

            guide_cutsite18 = guide_utils.calculate_cutsite18_guide(guide_strand, guide_start, guide_stop)
            cutsite = int(guide_cutsite18)

            # Searching for exons that this guide cuts
            # it's possible to have overlapping genes, and also overlapping transcripts where the same exon is
            # coding in one but non-coding in the other one!

            guide_exons = []
            if gene_name not in exon_dict.gene_exons:  # NMD genes for example
                logging.warning('not present in the parsed gtf %s', gene_name)
            else:
                for exon in exon_dict.gene_exons[gene_name]:  # a linear search is ok, only a handful of exons
                    exon_start = exon_dict.exon_gene_info[exon]['start']
                    exon_stop = exon_dict.exon_gene_info[exon]['stop']
                    # if exon_start <= exon_stop:
                    if exon_start <= cutsite <= exon_stop:
                        guide_exons.append(exon)

            guide_exons = sorted(guide_exons)
            logging.debug('gene %s, guide %s, guide_exon_list %s', gene_name, guide_seq, guide_exons)

            try:
                # Guide_uniqueness
                offtarget_profile_list = guide_offtarget_profile.split(",")
                guide_uniq = offtarget_profile_list[1]

                guide_seq_with_ngg = guide_utils.get_guide_seq_with_ngg(guide_start, guide_stop, guide_strand,
                                                                         guide_chr,
                                                                         sequence_dict)

                # Extracting feature
                cds_start_stop_cut_site = guide_utils.calculate_cutsite18_dist_fom_exon_cds_start_stop_for_exon_list_modified(
                    guide_exons, guide_cutsite18, exon_dict)
                exon_genomic_features = guide_utils.extract_exon_features_from_gtf_for_exonlist_modified(
                    guide_exons, exon_dict)

                if len(guide_exons) != 0:
                    microhomology_sequence = guide_utils.generate_seq_for_microhomology_scoring(guide_start,
                                                                                                guide_stop,
                                                                                                guide_strand,
                                                                                                guide_chr,
                                                                                                sequence_dict)
                    exon_biotype_list = [exon_dict.exon_id_biotype[x]["biotype"] for x in
                                         (guide_exons)]
                else:
                    microhomology_sequence = "nan"
                    exon_biotype_list = "nan"

                output = "\t".join(
                    [gene_name, guide_seq_with_ngg, guide_seq, guide_chr, guide_start, guide_stop, guide_strand,
                     guide_uniq, str(guide_exons), str(exon_biotype_list), str(guide_cutsite18),
                     "\t".join(cds_start_stop_cut_site), "\t".join(exon_genomic_features), microhomology_sequence])
                buffer.append(output)
                # gene_output_file_handle.write("\n")

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
    start = timeit.default_timer()
    num_processed = 0
    guides_directory = os.path.join('input', params['ensembl_release'], params['species_name'], 'guides')
    logging.info('reading guides from %s', guides_directory)
    # weird cases for testing:
    # for guide_file_path in ['ENSDARG00000100456.guides.txt']: # has single nucleotide exon
    # for guide_file_path in ['ENSG00000026036.guides.txt']: # NMD gene
    for guide_file_path in os.listdir(guides_directory):
        logging.debug('processing %s', guide_file_path)
        making_guide_file_with_info(os.path.join(guides_directory, guide_file_path), output_directory, sequence_dict,
                                    exon_dict)
        num_processed += 1
        if num_processed % 1000 == 0:
            stop = timeit.default_timer()
            logging.info('%d processed in %dsec', num_processed, stop - start)
    stop = timeit.default_timer()
    logging.info('time to generate guides %dsec', stop - start)


if __name__ == "__main__":
    if len(sys.argv) < 2 or not os.path.exists(sys.argv[1]):
        exit('missing config file!')

    logging.basicConfig(filename=None, level=getattr(logging, 'DEBUG', None),
                        format='%(asctime)s %(levelname)s %(funcName)s %(message)s')
    params = json.load(open(sys.argv[1]))
    logging.info("computing info files, parameters: %s", params)

    output_directory = prepareOutputDirectory(params)

    parsed_file_path = downloads.get_parsed_gtf_filepath(params)
    logging.info('loading exons from %s', parsed_file_path)
    exon_dict = ExonsInfo(parsed_file_path)

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

    logging.info('generating guides info')
    generate_guides(params, sequence_dict, output_directory, exon_dict)
