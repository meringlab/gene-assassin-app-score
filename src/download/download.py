import os
import sys
import json
import urllib2
import shutil
from subprocess import call
import logging

OUTPUT_DIR = "Raw_data_files"

def download_from_url(url_link, dst_dir):
    file_name = os.path.basename(url_link)
    destination = os.path.join(dst_dir, file_name)
    logging.info('downloading %s to %s', url_link, destination)

    if not os.path.exists(dst_dir):
        os.makedirs(dst_dir)

    with urllib2.urlopen(url_link) as file_handle:
        with open(destination, "w") as f:
            shutil.copyfileobj(file_handle, f)
    logging.info('saved at %s', destination)

    return destination


def download(params):
    dst_dir = os.path.join("output" , params['ensembl_release'], params['species_name'], OUTPUT_DIR)

    download_from_url(params['base_url_gtf'], dst_dir)

    destination = download_from_url(params['DNA_top_level_fa'], dst_dir)
    # reading gzip is much much slower than uncompressing and reading uncompressed:
    call(["gunzip", destination])

    destination = download_from_url(params['GVF_file'], dst_dir)
    call(["gunzip", destination])


if __name__ == "__main__":
    if len(sys.argv) < 2 or not os.path.exists(sys.argv[1]):
        exit('missing config file!')

    params = json.load(open(sys.argv[1]))
    download(params)

# base_url_gtf_1 = "ftp://ftp.ensembl.org/pub/release-90/gtf/danio_rerio/Danio_rerio.GRCz10.90.gtf.gz"
# DNA_top_level_fa = "ftp://ftp.ensembl.org/pub/release-90/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz"
# GVF_file = "ftp://ftp.ensembl.org/pub/release-90/variation/gvf/danio_rerio/Danio_rerio.gvf.gz"
# ensembl_version = "v_90"
