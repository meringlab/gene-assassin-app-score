import os
import sys
import json
from subprocess import call
import logging

OUTPUT_DIR = "Raw_data_files"

def download_from_url(url_link, dst_dir):
    logging.info('downloading %s to %s', url_link, dst_dir)
    file_name = os.path.basename(url_link)
    destination = os.path.join(dst_dir, file_name)

    if not os.path.exists(dst_dir):
        os.makedirs(dst_dir)

    call(["rsync", "-av", url_link, dst_dir])

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
