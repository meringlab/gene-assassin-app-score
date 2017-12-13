from Bio.Seq import Seq


class Guide(object):
    _seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

    def __init__(self, tsv):
        self._load_from_tsv(tsv)

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

    def __str__(self):
        return '%s %s:%d-$d' % (self.strand, self.seq, self.start, self.end)

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
