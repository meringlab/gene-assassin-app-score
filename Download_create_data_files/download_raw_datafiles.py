import os
import sys
import json
import urllib2
import shutil

def download_from_url(url_link, downlading_path, species_name, file_category):
    file_name = url_link.split("/")[-1]
    full_path = os.path.join(downlading_path, species_name, file_category)
    destination = os.path.join(full_path, file_name)
    print('downloading %s to %s' % (url_link, destination))

    if not os.path.exists(full_path):
        try:
            os.makedirs(full_path)
        except Exception as e:
            print e
            exit("Cannot create directory %s" % downlading_path)

    try:
        file_link = url_link
        file_handle = urllib2.urlopen(file_link)
    except Exception as e:
        print e
        exit("\t ... URL link does not work %s." % url_link)

    with open(destination, "w") as f:
        shutil.copyfileobj(file_handle, f)

    file_handle.close()

    print('\t... saved at %s \n' % destination)


if __name__ == "__main__":
    if len(sys.argv) < 2 or not os.path.exists(sys.argv[1]):
        exit('missing config file!')
    params = json.load(open(sys.argv[1]))

    tempfile_path = "output/%s" % params['ensembl_release']
    # ensembl_version = "v_85"
    file_category = "Raw_data_files"

    download_from_url(params['base_url_gtf'], tempfile_path, params['species_name'], file_category)
    download_from_url(params['DNA_top_level_fa'], tempfile_path, params['species_name'],  file_category)
    download_from_url(params['GVF_file'], tempfile_path, params['species_name'], file_category)

#base_url_gtf_1 = "ftp://ftp.ensembl.org/pub/release-90/gtf/danio_rerio/Danio_rerio.GRCz10.90.gtf.gz"
#DNA_top_level_fa = "ftp://ftp.ensembl.org/pub/release-90/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz"
#GVF_file = "ftp://ftp.ensembl.org/pub/release-90/variation/gvf/danio_rerio/Danio_rerio.gvf.gz"
#ensembl_version = "v_90"
