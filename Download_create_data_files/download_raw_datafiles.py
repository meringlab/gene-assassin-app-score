import download_raw_datafiles_function as downaload_fn


############### Danio_rerio

####### v 85

species_name = "danio_rerio"
base_url_gtf_85 = "ftp://ftp.ensembl.org/pub/release-85/gtf/danio_rerio/Danio_rerio.GRCz10.85.gtf.gz"
DNA_top_level_fa_85 = "ftp://ftp.ensembl.org/pub/release-85/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz"
GVF_file_85= "ftp://ftp.ensembl.org/pub/release-85/variation/gvf/danio_rerio/Danio_rerio.gvf.gz"
tempfile_path = "output/v85_2017"
ensembl_version = "v_85"
file_category = "Raw_data_files"

downaload_fn.download_datafile_ensembl(base_url_gtf_85, tempfile_path,species_name,ensembl_version,file_category)
downaload_fn.download_datafile_ensembl(DNA_top_level_fa_85, tempfile_path, species_name,ensembl_version,file_category)  
downaload_fn.download_datafile_ensembl(GVF_file_85, tempfile_path, species_name, ensembl_version,file_category)  

