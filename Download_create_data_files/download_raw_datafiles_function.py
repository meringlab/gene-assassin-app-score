import os, sys
import urllib
import urllib2
import re
import gzip
import shutil


def download_datafile_ensembl(url_link, downlading_path, species_name, release_version,file_category):
    
    try:
        file_link = url_link
        file_handle = urllib2.urlopen(file_link)

    except Exception as e:
        print e
        exit("\t ... URL link does not work %s."%url_link)
        
    ####### Cp the file

    file_name = url_link.split("/")[-1]
    
    downlading_path_species= os.path.join(downlading_path, species_name)
    downlading_path_version_name =os.path.join(downlading_path_species, str(release_version))
    downlading_path_file_category = os.path.join(downlading_path_version_name, file_category)
    
    
    
    
    ########## Making data files diretcory
    
    if not os.path.exists(downlading_path):
        try :
            os.makedirs(downlading_path)
        except Exception as e:
            print e 
            exit("Cannot create directory %s" %downlading_path)
            
    
    if not os.path.exists(downlading_path_species):
        try :
            os.makedirs(downlading_path_species)
        except Exception as e:
            print e 
            exit("Cannot create directory %s" %downlading_path_species)
            
    
    if not os.path.exists(downlading_path_version_name):
        try :
            os.makedirs(downlading_path_version_name)
        except Exception as e:
            print e 
            exit("Cannot create directory %s" %downlading_path_version_name)
            
    if not os.path.exists(downlading_path_file_category):
        try :
            os.makedirs(downlading_path_file_category)
        except Exception as e:
            print e 
            exit("Cannot create directory %s" %downlading_path_file_category)
            
    
    
    
    tempfile_name = os.path.join(downlading_path_file_category, file_name)
    tempfile_file_hande = open(tempfile_name, "w")
    
    shutil.copyfileobj(file_handle, tempfile_file_hande)
    
    tempfile_file_hande.close()
    file_handle.close()
    
    sys.stdout.write('\t... saved at %s \n' % tempfile_name)
    

#species_name = "danio_rerio"
#base_url_gtf_1 = "ftp://ftp.ensembl.org/pub/release-90/gtf/danio_rerio/Danio_rerio.GRCz10.90.gtf.gz"
#DNA_top_level_fa = "ftp://ftp.ensembl.org/pub/release-90/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz"
#GVF_file = "ftp://ftp.ensembl.org/pub/release-90/variation/gvf/danio_rerio/Danio_rerio.gvf.gz"
#tempfile_path = "/Users/neha/Desktop/Crispr_dummy/Raw_data_files"
#ensembl_version = "v_90"
#
#download_datafile_ensembl(base_url_gtf_1, tempfile_path, ensembl_version,species_name)
#download_datafile_ensembl(DNA_top_level_fa, tempfile_path, ensembl_version,species_name)  
#download_datafile_ensembl(GVF_file, tempfile_path, ensembl_version,species_name)  
#    
    
    
    