import unittest

from helper.guide import Guide


class GuideTest(unittest.TestCase):

    def test_parse_simple(self):
        g = Guide.from_simple_tsv(
            '4	15012700	15012720	1	AAAGATGTCGTCTCCAAGTC	0.3791	0.549	0.995	0.7826	110,145,0	2,1,0,0,0')
        self.assertGuideEqualBasic(g, ['4', 15012700, 15012720, '1', 'AAAGATGTCGTCTCCAAGTC', '1'])

        g = Guide.from_simple_tsv(
            '4	15052112	15052132	-1	TTCCAAACGCCCCCTTCATA	0.3057	0.3467	1.0	0.7305	138,117,0	1,1,0,0,0')
        self.assertGuideEqualBasic(g, ['4', 15052112, 15052132, '-1', 'TTCCAAACGCCCCCTTCATA', '1'])

    def test_parse_full(self):
        full_tsv = "ENSDARG00000000001	AGAAGAAGAGTCTTTCTGAGTGG	AGAAGAAGAGTCTTTCTGAG	9	34305562	34305582	-1	1	['ENSDARE00000000004']	['coding']	[68.0, 68.0]	[3.0, 3.0]	2	['ENSDART00000000004', 'ENSDART00000169788']	['7', '6']	AACCAAGGCTAACACAGAGAGTGAATCAGAAGAAGAGTCTTTCTGAGTGGGATTTGAACTC"
        g = Guide.from_full_tsv(full_tsv)
        self.assertGuideEqualBasic(g, ['9', 34305562, 34305582, '-1', 'AGAAGAAGAGTCTTTCTGAG', '1'])
        expected = ['ENSDARG00000000001', 'AGAAGAAGAGTCTTTCTGAGTGG', "['ENSDARE00000000004']", "['coding']",
                    '[68.0, 68.0]',
                    '[3.0, 3.0]', 2, "['ENSDART00000000004', 'ENSDART00000169788']", "['7', '6']",
                    "AACCAAGGCTAACACAGAGAGTGAATCAGAAGAAGAGTCTTTCTGAGTGGGATTTGAACTC"]
        self.assertGuideEqualExtra(g, expected)
        self.assertEquals(full_tsv, g.to_tsv())
        self.assertGuideEqualExtra(Guide.from_full_tsv(g.to_tsv()), expected)

    def test_pam(self):
        g = Guide()
        g.chromosome, g.start, g.end, g.strand = ('1', 28, 47, '1')
        chr = {'1': 'AACCAAGGCTAACACAGAGAGTGAATCAGAAGAAGAGTCTTTCTGAGTGGGATTTGAACTC'}
        g.seq = 'AGAAGAAGAGTCTTTCTGAG'
        self.assertEqual('AGAAGAAGAGTCTTTCTGAGTGG', g.seq_with_pam(chr))

        g.chromosome, g.start, g.end, g.strand = ('1', 6, 25, '-1')
        chr = {'1': 'AACCAAGGCTAACACAGAGAGTGAATCAGAAGAAGAGTCTTTCTGAGTGGGATTTGAACTC'}
        g.seq = 'TTCACTCTCTGTGTTAGCCT'

        self.assertEqual('TTCACTCTCTGTGTTAGCCTTGG', g.seq_with_pam(chr))

        # if g.is_on_forward_strand():
        #     kmer_start = exon_start + position - 37
        #     kmer_end = kmer_start + kmer_length
        #     kmer_d = str(dna_dict[chromosome].seq[kmer_start:kmer_end])
        #     pam_d = str(dna_dict[chromosome].seq[kmer_end:kmer_end + pam_length])
        # else:
        #     kmer_start = exon_end - position + 16
        #     kmer_end = kmer_start + kmer_length
        #     kmer_d = str(dna_dict[chromosome].seq[kmer_start:kmer_end].reverse_complement())
        #     pam_d = str(dna_dict[chromosome].seq[kmer_start - pam_length:kmer_start].reverse_complement())

        ## test
        # seq_ngg_regex = "[ATGC]{20}[ATGC]{1}GG"
        # search_object_ngg = re.search(seq_ngg_regex, seq)
        # if search_object_ngg:
        #     ngg_outcome = seq
        # else:
        #     ngg_outcome = "nan"


    def assertGuideEqualExtra(self, guide, values):
        self.assertEqual(values[0], guide.gene_id)
        self.assertEqual(values[1], guide.pamseq)
        self.assertEqual(values[2], guide.guide_exons)
        self.assertEqual(values[3], guide.exon_biotype_list)
        self.assertEqual(values[4], guide.cds_start_cutsite)
        self.assertEqual(values[5], guide.cds_stop_cutsite)
        self.assertEqual(values[6], guide.transcript_count)
        self.assertEqual(values[7], guide.transcripts)
        self.assertEqual(values[8], guide.exon_rank_list)
        self.assertEqual(values[9], guide.microhomology_sequence)

    def assertGuideEqualBasic(self, guide, values):
        self.assertEqual(values[0], guide.chromosome)
        self.assertEqual(values[1], guide.start)
        self.assertEqual(values[2], guide.end)
        self.assertEqual(values[3], guide.strand)
        self.assertEqual(values[4], guide.seq)
        self.assertEqual(values[5], guide.uniqueness)
