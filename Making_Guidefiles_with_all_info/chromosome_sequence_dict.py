#!/usr/bin/env python
import sys
import gzip
import timeit

def load_chromosome_sequence_dict(fasta_filepath):
    '''
    TODO replace this with biopython:

    from Bio import SeqIO
    with gzip.open(fasta_file, "rt") as handle:  # only "r" doesnt work in python 3:
        chromosome_sequence_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    :param fasta_filepath:
    :return:
    '''
    # chromosome_data_path = "output/v85_2017/danio_rerio/v_85/Raw_data_files/Danio_rerio.GRCz10.dna.toplevel.fa.gz"

    # if not os.path.exists(fasta_filepath):
    #     exit("\t ... fasta file does not exist %s." % fasta_filepath)
    start = timeit.default_timer()

    try:
        input_file_handle = gzip.open(fasta_filepath, "rb")
    except:
        print "failed to read fasta file: ", fasta_filepath, sys.exc_info()[0]
        raise

    seq= ""
    chromosome_sequence_dict = {}

    for line in input_file_handle:
        if line.startswith(">"):
            if seq!= "":
                #print l ,chromosome_name, len(seq)  ####### Good check to see if everything is rightly stored
                if chromosome_name not in chromosome_sequence_dict:
                    chromosome_sequence_dict[chromosome_name] = seq

            seq = ""

            l = line.strip("\n").strip(">").split(" ")
            chromosome_name = l[0]

        else:
            seq = seq + line.strip("\n")

    chromosome_sequence_dict[chromosome_name] = seq

    stop = timeit.default_timer()
    print('time to build sequence index %dsec' % (stop - start))

    return chromosome_sequence_dict

