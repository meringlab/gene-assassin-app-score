import os


# parsed_gtf_file_path = "output/v_85/danio_rerio/Processed_data_files/parsed_Danio_rerio.GRCz10.85.gtf.txt"

class ExonsInfo(object):
    def __init__(self, parsed_gtf_file_path):
        if not os.path.exists(parsed_gtf_file_path):
            raise Exception("parsed gtf file path does not exist %s." % parsed_gtf_file_path)

        self.exon_gene_info = {}
        self.gene_exons = {}
        self.chromosomes = set()
        self.exon_transcript_info = {}
        self.exon_id_biotype = {}

        self._load_(parsed_gtf_file_path)

    def _load_(self, parsed_gtf_file_path):
        input_file_handle = open(parsed_gtf_file_path)

        for line in input_file_handle:

            l = line.strip("\n").split("\t")

            gene_id = l[0];
            gene_strand = l[2];
            total_transcript = l[3]
            transcript_id = l[4];
            total_exons = l[5]
            exon_id = l[6];
            exon_location = l[7];
            exon_rank = l[-2];
            exon_cds_location = l[-1]

            if gene_id != "Gene_id":
                chr = exon_location[:exon_location.index(':')]
                self.chromosomes.add(chr)
                if gene_id not in self.gene_exons:
                    self.gene_exons[gene_id] = set()
                self.gene_exons[gene_id].add(exon_id)

                ######## Exon_gene_info :
                if exon_id not in self.exon_gene_info:
                    self.exon_gene_info[exon_id] = {}

                start, stop = exon_location[exon_location.index(':') + 1:].split('-')
                self.exon_gene_info[exon_id]["start"] = int(start)
                self.exon_gene_info[exon_id]["stop"] = int(stop)
                self.exon_gene_info[exon_id]["gene_strand"] = gene_strand
                self.exon_gene_info[exon_id]["total_transcript_per_gene"] = total_transcript

                ########## Exon_id_biotype
                if exon_id not in self.exon_id_biotype:
                    self.exon_id_biotype[exon_id] = {}

                if exon_cds_location == "nan":
                    self.exon_id_biotype[exon_id]["biotype"] = "non-coding"
                else:
                    self.exon_id_biotype[exon_id]["biotype"] = "coding"

                ######## Exon_transcipt_info
                if exon_id not in self.exon_transcript_info:
                    self.exon_transcript_info[exon_id] = {}

                    self.exon_transcript_info[exon_id]["transcript_list"] = []
                    self.exon_transcript_info[exon_id]["exon_rank_list"] = []
                    self.exon_transcript_info[exon_id]["exon_cds_list"] = []
                    self.exon_transcript_info[exon_id]["total_exons_in_transcript"] = []

                # if exon_cds_location != "nan":                             ####### only_protein_coding info
                if transcript_id not in self.exon_transcript_info[exon_id]["transcript_list"]:
                    self.exon_transcript_info[exon_id]["transcript_list"].append(transcript_id)
                    self.exon_transcript_info[exon_id]["exon_rank_list"].append(exon_rank)
                    self.exon_transcript_info[exon_id]["exon_cds_list"].append(exon_cds_location)
                    self.exon_transcript_info[exon_id]["total_exons_in_transcript"].append(total_exons)

# input_file_handle = open(parsed_gtf_file_path)
#
# location_exon_id = {}
# exon_gene_info = {}
# exon_transcript_info = {}
# exon_id_biotype = {}
#
#
# for line in input_file_handle:
#
#     l = line.strip("\n").split("\t")
#
#     gene_id = l[0]; gene_strand = l[2]; total_transcript = l[3]
#     transcript_id = l[4]; total_exons = l[5]
#     exon_id = l[6]; exon_location = l[7]; exon_rank = l[-2]; exon_cds_location = l[-1]
#
#
#     #print gene_id, gene_strand, total_transcript
#     #print transcript_id, total_exons
#     #print exon_id, exon_location, exon_rank, exon_cds_location
#
#     if gene_id!= "Gene_id":
#
#
#         ###### Location_exon_id
#
#
#
#         if exon_location not in location_exon_id:
#
#             location_exon_id[exon_location] = []
#
#         if exon_id not in location_exon_id[exon_location]:
#             location_exon_id[exon_location].append(exon_id)
#
#
#
#         ######## Exon_gene_info :
#
#
#         if exon_id not in exon_gene_info:
#             exon_gene_info[exon_id] = {}
#
#         exon_gene_info[exon_id]["gene_id"] = gene_id
#         exon_gene_info[exon_id]["gene_strand"] = gene_strand
#         exon_gene_info[exon_id]["total_transcript_per_gene"] = total_transcript
#
#
#         ########## Exon_id_biotype
#
#         if exon_id not in exon_id_biotype:
#             exon_id_biotype[exon_id] = {}
#
#         if exon_cds_location == "nan":
#             exon_id_biotype[exon_id]["biotype"] = "non-coding"
#         else:
#             exon_id_biotype[exon_id]["biotype"] = "coding"
#
#
#
#         ######## Exon_transcipt_info
#
#
#         if exon_id not in exon_transcript_info:
#             exon_transcript_info[exon_id] = {}
#
#             exon_transcript_info[exon_id]["transcript_list"] = []
#             exon_transcript_info[exon_id]["exon_rank_list"] = []
#             exon_transcript_info[exon_id]["exon_cds_list"] = []
#             exon_transcript_info[exon_id]["total_exons_in_transcript"] = []
#
#         #if exon_cds_location != "nan":                             ####### only_protein_coding info
#
#         if transcript_id not in exon_transcript_info[exon_id]["transcript_list"]:
#
#             exon_transcript_info[exon_id]["transcript_list"].append(transcript_id)
#             exon_transcript_info[exon_id]["exon_rank_list"].append(exon_rank)
#             exon_transcript_info[exon_id]["exon_cds_list"].append(exon_cds_location)
#             exon_transcript_info[exon_id]["total_exons_in_transcript"].append(total_exons)
#
