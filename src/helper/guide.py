from Bio.Seq import Seq
import logging
import ast


def extract_exon_features_from_gtf(exon_id, exon_dict):
    if exon_id in exon_dict.exon_transcript_info:
        total_transcript_count = exon_dict.exon_gene_info[exon_id]["total_transcript_per_gene"]  ###### gene_feature

        transcript_list = exon_dict.exon_transcript_info[exon_id]["transcript_list"]
        exon_rank_list = exon_dict.exon_transcript_info[exon_id]["exon_rank_list"]
        total_no_exon_transcript_list = exon_dict.exon_transcript_info[exon_id]["total_exons_in_transcript"]

        return (total_transcript_count, str(transcript_list), str(total_no_exon_transcript_list), str(exon_rank_list))
    else:
        return ("nan", "nan", "nan", "nan")


def exon_features_from_gtf(exon_list, exon_dict):
    if len(exon_list) == 0:  # no protein coding trsnacript
        total_transcipt_count = "nan"
        transcript_list_all = []
        exon_rank_list_all = []
    else:
        transcript_list_all = []
        exon_rank_list_all = []

        for exon in exon_list:
            exon_biotype = exon_dict.exon_id_biotype[exon]["biotype"]
            exon_feautre = extract_exon_features_from_gtf(exon, exon_dict)

            if exon_biotype == "coding":
                transcript_list_all = transcript_list_all + ast.literal_eval(exon_feautre[1])
                exon_rank_list_all = exon_rank_list_all + ast.literal_eval(exon_feautre[3])
            else:
                no_of_transcript_covered = len(ast.literal_eval(exon_feautre[1]))
                list_exon_rank_insert = [0 for x in range(no_of_transcript_covered)]
                exon_rank_list_all = exon_rank_list_all + list_exon_rank_insert

        total_transcipt_count = exon_feautre[0]  # since it is gene feature remains same for all the exons

    return (total_transcipt_count, str(transcript_list_all), str(exon_rank_list_all))


class Guide(object):
    _seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

    full_info_fields_order = ('gene_id',  # 1
                              'pamseq',  # 2
                              'seq',
                              'chromosome',
                              'start',  # 5
                              'end',
                              'strand',
                              'uniqueness',
                              'guide_exons',
                              'exon_biotype_list',  # 10
                              'cds_start_cutsite',
                              'cds_stop_cutsite',
                              'transcript_count',
                              'transcripts',
                              'exon_rank_list',  # 15
                              'microhomology_sequence'
                              )

    types = {
        'start': int,
        'end': int,
        'transcript_count': int,
    }

    def __init__(self):
        # for f in Guide.full_info_fields_order:
        #     self[f] = ''
        pass

    @classmethod
    def full_tsv_header(cls):
        return '\t'.join([f.capitalize() for f in Guide.full_info_fields_order])

    @classmethod
    def from_simple_tsv(cls, line):
        g = cls()
        g._load_from_tsv(line)
        return g

    @classmethod
    def from_full_tsv(cls, line):
        g = cls()
        g._load_from_full_info_tsv(line)
        return g

    @classmethod
    def load_guide(cls, gene_name, guide_tsv, exon_dict, sequence_dict):
        g = cls.from_simple_tsv(guide_tsv)

        if g.chromosome not in exon_dict.chromosomes:
            # this happened for human (and mouse), chromosome CHR_HSCHR6_MHC_QBL_CTG1 (in fasta but not in gtf)
            logging.warning("unknown chromosome %s, for %s", g.chromosome, guide_tsv)
            return None

        try:
            guide_exons = g.find_target_exons(gene_name, exon_dict)
            logging.debug('gene %s, guide %s, guide_exon_list %s', gene_name, g.seq, guide_exons)
            g.calculate_cutsite_distance_fom_exon_cds(guide_exons, exon_dict)
            g.set_microhomology_sequence(guide_exons, sequence_dict)
            g.set_exon_biotype_list(guide_exons, exon_dict)
            (transcript_count, transcripts, exon_ranks) = exon_features_from_gtf(guide_exons, exon_dict)
            g.set_info(gene_name, g.seq_with_pam(sequence_dict), guide_exons, transcript_count, transcripts, exon_ranks)
            return g
        except Exception as e:
            logging.exception('failed to collect info on %s - %s, error: %s', gene_name, guide_tsv, e)
        return None

    def find_target_exons(self, gene, exon_dict):
        '''
          Searching for exons that this guide cuts.

          It's possible to have overlapping genes, and also overlapping transcripts where the same exon is
          coding in one but non-coding in the other one!
        :param guide:
        :param gene:
        :param exon_dict:
        :return:
        '''
        guide_exons = []
        if gene in exon_dict.gene_exons:
            for exon in exon_dict.gene_exons[gene]:  # a linear search is ok, only a handful of exons
                exon_start = exon_dict.exon_gene_info[exon]['start']
                exon_stop = exon_dict.exon_gene_info[exon]['stop']
                # if exon_start <= exon_stop:
                if self.is_cutsite_within(exon_start, exon_stop):
                    guide_exons.append(exon)

        return sorted(guide_exons)

    def _load_from_full_info_tsv(self, line):
        r = line.strip().split('\t')
        if len(r) != len(Guide.full_info_fields_order):
            raise Exception('invalid record %s' % line)

        for i in range(len(Guide.full_info_fields_order)):
            field = Guide.full_info_fields_order[i]
            if field in Guide.types and r[i] != 'nan':
                r[i] = Guide.types[field](r[i])
            setattr(self, field, r[i])

    def to_tsv(self):
        return '\t'.join([str(getattr(self, f)) for f in Guide.full_info_fields_order])

    def _load_from_tsv(self, line):
        l = line.strip().split('\t')
        self.chromosome = l[0]
        self.start = int(l[1])
        self.end = int(l[2])
        self.strand = l[3]
        self.seq = l[4]
        self.offtarget_profile = l[10]
        self.uniqueness = self.offtarget_profile.split(",")[1]
        self.cutsite = self._compute_cutsite()

    def set_info(self, gene_id, seq_with_pam, guide_exons,
                 transcript_count, transcript_id_list, exon_rank_list):
        self.gene_id = gene_id
        self.pamseq = seq_with_pam
        self.transcript_count = transcript_count
        self.transcripts = transcript_id_list
        self.exon_rank_list = exon_rank_list
        self.guide_exons = guide_exons

    def set_microhomology_sequence(self, guide_exons, sequence_dict):
        self.microhomology_sequence = 'nan'
        if guide_exons:
            self.microhomology_sequence = self.seq_for_microhomology_scoring(sequence_dict)

    def set_exon_biotype_list(self, guide_exons, exon_dict):
        self.exon_biotype_list = 'nan'
        if guide_exons:
            self.exon_biotype_list = [exon_dict.exon_id_biotype[x]["biotype"] for x in guide_exons]

    def _calculate_cutsite18_dist_fom_exon_cds_start_stop_list(self, guide_exon, exon_dict):
        # if guide_exon in exon_dict.exon_transcript_info:
        #     return ((),())
        exon_cds_start_stop_list = exon_dict.exon_transcript_info[guide_exon]["exon_cds_list"]
        exon_strand = exon_dict.exon_gene_info[guide_exon]["gene_strand"]

        exon_cds_start_stop_list_coding = [x for x in exon_cds_start_stop_list if x != "nan"]

        # split cds for coding exons
        exon_cds_start_list = [float(x.split(":")[1].split("-")[0]) for x in exon_cds_start_stop_list_coding]
        exon_cds_stop_list = [float(x.split(":")[1].split("-")[-1]) for x in exon_cds_start_stop_list_coding]

        # add nan to the distance
        if "nan" in exon_cds_start_stop_list:
            no_transcript_covered = len(exon_cds_start_stop_list)

            dist_cds_start_cutsite_list_noncoding = ["nan" for x in range(no_transcript_covered)]
            dist_cds_stop_cutsite_list_noncoding = ["nan" for x in range(no_transcript_covered)]

            dist_cds_start_cutsite_list_coding = []
            dist_cds_stop_cutsite_list_coding = []
        else:
            dist_cds_start_cutsite_list_noncoding = []
            dist_cds_stop_cutsite_list_noncoding = []

            # Calculating distances
            if exon_strand == "-":  ####### real cds start is stop according to annotation
                dist_cds_start_cutsite_list_coding = [abs(self.cutsite - x) for x in exon_cds_stop_list]
                dist_cds_stop_cutsite_list_coding = [abs(self.cutsite - x) for x in exon_cds_start_list]

            if exon_strand == "+":
                dist_cds_start_cutsite_list_coding = [abs(self.cutsite - x) for x in exon_cds_start_list]
                dist_cds_stop_cutsite_list_coding = [abs(self.cutsite - x) for x in exon_cds_stop_list]

        dist_cds_start_cutsite_list = dist_cds_start_cutsite_list_coding + dist_cds_start_cutsite_list_noncoding
        dist_cds_stop_cutsite_list = dist_cds_stop_cutsite_list_coding + dist_cds_stop_cutsite_list_noncoding

        return (dist_cds_start_cutsite_list, dist_cds_stop_cutsite_list)

    def calculate_cutsite_distance_fom_exon_cds(self, exon_list, exon_dict):
        self.cds_start_cutsite = []
        self.cds_stop_cutsite = []
        for exon in exon_list:
            distance_feature = self._calculate_cutsite18_dist_fom_exon_cds_start_stop_list(exon, exon_dict)
            self.cds_start_cutsite.extend(distance_feature[0])
            self.cds_stop_cutsite.extend(distance_feature[1])

    def __str__(self):
        return '%s %s:%d-$d' % (self.strand, self.seq, self.start, self.end)

    def __repr__(self):
        return self.__str__()

    def is_on_forward_strand(self):
        return self.strand == "+" or self.strand == "1"

    def _compute_cutsite(self):
        if self.is_on_forward_strand():
            return self.end - 2
        return self.start + 2

    # def set_pam(self, sequence):
    #     self.ngg = sequence
    #
    # def get_pam(self):
    #     return self.ngg
    #
    # pam = property(get_pam, set_pam)

    def is_cutsite_within(self, start, stop):
        return start <= self.cutsite <= stop

    def _reverse_complement(self, seq):
        # return "".join([Guide._seq_dict[base] for base in reversed(seq)])
        return str(Seq(seq).reverse_complement())

    def seq_with_pam(self, chromosome_sequence_dict):
        chr_seq = chromosome_sequence_dict[self.chromosome]

        if self.is_on_forward_strand():
            return chr_seq[self.start - 1: self.end + 3]
        seq = chr_seq[self.start - 3 - 1:self.end]
        seq = self._reverse_complement(seq)

        ##### test
        # seq_ngg_regex = "[ATGC]{20}[ATGC]{1}GG"
        # search_object_ngg = re.search(seq_ngg_regex, seq)
        # if search_object_ngg:
        #     ngg_outcome = seq
        # else:
        #     ngg_outcome = "nan"

        return seq

    def seq_for_microhomology_scoring(self, chromosome_sequence_dict):
        sequence_start = self.start - 13
        sequence_end = self.end + 27

        chr_seq = chromosome_sequence_dict[self.chromosome]
        if self.is_on_forward_strand():
            return chr_seq[sequence_start - 1:sequence_end]
        return self._reverse_complement(chr_seq[sequence_start - 1:sequence_end])
