#!/usr/bin/env python
import sys
import gzip
import timeit
from Bio import SeqIO
from cStringIO import StringIO
import io


def load_chromosome_sequence_dict_bio(fasta_filepath):
    start = timeit.default_timer()
    with gzip.open(fasta_filepath, "rt") as handle:  # only "r" doesnt work in python 3:
        chromosome_sequence_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    stop = timeit.default_timer()
    print('time to build sequence index using SeqIO %dsec' % (stop - start))
    return chromosome_sequence_dict


def load_chromosome_sequence_dict(fasta_filepath):
    '''
    :param fasta_filepath: "output/v_85/danio_rerio/Raw_data_files/Danio_rerio.GRCz10.dna.toplevel.fa.gz"
    :return: dict { chromosome -> dna sequence }
    '''

    # if not os.path.exists(fasta_filepath):
    #     exit("\t ... fasta file does not exist %s." % fasta_filepath)
    start = timeit.default_timer()

    try:
        gz = gzip.open(fasta_filepath, "rb")
    except:
        print "failed to read fasta file: ", fasta_filepath, sys.exc_info()[0]
        raise

    seq= ""
    chromosome_sequence_dict = {}
    input_file_handle = io.BufferedReader(gz)
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

def load_chromosome_sequence_dict_from_fasta(fasta_uncompressed_filepath):
    '''
    this doesnt measure any faster than plain string concatenation

    :param fasta_uncompressed_filepath: "output/v_85/danio_rerio/Raw_data_files/Danio_rerio.GRCz10.dna.toplevel.fa"
    :return: dict { chromosome -> dna sequence }
    '''


    start = timeit.default_timer()
    buffer= StringIO()
    chromosome_sequence_dict = {}

    with open(fasta_uncompressed_filepath, "r") as input_file_handle:
        for line in input_file_handle:
            line = line.strip()
            if line[0] == '>':
                seq = buffer.getvalue()
                if seq!= '':
                    #print(l, chromosome_name, len(seq))  # test if everything is rightly stored
                    if chromosome_name not in chromosome_sequence_dict:
                        chromosome_sequence_dict[chromosome_name] = seq

                buffer = StringIO() # or just .seek(0) .truncate() ?

                l = line.strip(">").split(" ")
                chromosome_name = l[0]
            else:
                buffer.write(line)

    seq = buffer.getvalue()
    chromosome_sequence_dict[chromosome_name] = seq

    stop = timeit.default_timer()
    print('time to build sequence index using StringIO %dsec' % (stop - start))

    return chromosome_sequence_dict

def load_chromosome_sequence_dict_list(fasta_filepath):
    '''
    this doesnt measure any faster than plain string concatenation

    :param fasta_filepath: "output/v_85/danio_rerio/Raw_data_files/Danio_rerio.GRCz10.dna.toplevel.fa.gz"
    :return: dict { chromosome -> dna sequence }
    '''

    # if not os.path.exists(fasta_filepath):
    #     exit("\t ... fasta file does not exist %s." % fasta_filepath)
    start = timeit.default_timer()

    try:
        gz = gzip.open(fasta_filepath, "rb")
    except:
        print "failed to read fasta file: ", fasta_filepath, sys.exc_info()[0]
        raise

    buffer = []
    chromosome_sequence_dict = {}

    input_file_handle = io.BufferedReader(gz)
    for line in input_file_handle:
        if line.startswith(">"):
            seq = ''.join(buffer)
            if seq!= "":
                #print l ,chromosome_name, len(seq)  ####### Good check to see if everything is rightly stored
                if chromosome_name not in chromosome_sequence_dict:
                    chromosome_sequence_dict[chromosome_name] = seq

            # python 3 has .clear()
            # del buffer[:]
            buffer = []

            l = line.strip("\n").strip(">").split(" ")
            chromosome_name = l[0]

        else:
            buffer.append(line.strip())

    seq = ''.join(buffer)
    chromosome_sequence_dict[chromosome_name] = seq

    stop = timeit.default_timer()
    print('time to build sequence index using join on list %dsec' % (stop - start))

    return chromosome_sequence_dict

def test_read_gzip(fasta_filepath):
    start = timeit.default_timer()
    with gzip.open(fasta_filepath, "rb") as handle:
        for line in handle:
            if line[0] == '>':
                pass
    stop = timeit.default_timer()
    print('test time on %s %dsec' % (fasta_filepath, stop - start))

def test_read_gzip_buffer(fasta_filepath):
    c = 0
    start = timeit.default_timer()
    with gzip.open(fasta_filepath, "rb") as gz:
        f = io.BufferedReader(gz)
        for line in f:
            if line[0] == '>':
                c+=1
    stop = timeit.default_timer()
    print('test buffer read time on %s %dsec' % (fasta_filepath, stop - start))


if __name__ == '__main__':
    # test_read_gzip('../output/v_85/homo_sapiens/Raw_data_files/Homo_sapiens.GRCh38.dna.toplevel.fa.gz')
    # test_read_gzip('../output/v_85/danio_rerio/Raw_data_files/Danio_rerio.GRCz10.dna.toplevel.fa.gz')
    # output:
    # test time on ../output/v_85/homo_sapiens/Raw_data_files/Homo_sapiens.GRCh38.dna.toplevel.fa.gz 1451sec
    # test time on ../output/v_85/danio_rerio/Raw_data_files/Danio_rerio.GRCz10.dna.toplevel.fa.gz 50sec
    #
    # time { gunzip Homo_sapiens.GRCh38.dna.toplevel.fa.gz ; cat Homo_sapiens.GRCh38.dna.toplevel.fa | wc -l; }
    # 797852943
    # real	1m56.876s
    #
    # conclusion: READING GZIP DOMINATES THE PERFORMANCE

    #
    test_read_gzip('../output/v_85/danio_rerio/Raw_data_files/Danio_rerio.GRCz10.dna.toplevel.fa.gz')
    test_read_gzip_buffer('../output/v_85/danio_rerio/Raw_data_files/Danio_rerio.GRCz10.dna.toplevel.fa.gz')
    # test_read_gzip_buffer('Homo_sapiens.GRCh38.dna.toplevel.fa.gz')
    # test buffer read time on Homo_sapiens.GRCh38.dna.toplevel.fa.gz 576sec