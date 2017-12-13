import unittest

from helper.guide import Guide


class GuideTest(unittest.TestCase):

    def test_parse(self):
        g = Guide(
            '4	15012700	15012720	1	AAAGATGTCGTCTCCAAGTC	0.3791	0.549	0.995	0.7826	110,145,0	2,1,0,0,0')
        self.assertGuideEqual(g, ['4', 15012700, 15012720, '1', 'AAAGATGTCGTCTCCAAGTC', '1'])

        g = Guide(
            '4	15052112	15052132	-1	TTCCAAACGCCCCCTTCATA	0.3057	0.3467	1.0	0.7305	138,117,0	1,1,0,0,0')
        self.assertGuideEqual(g, ['4', 15052112, 15052132, '-1', 'TTCCAAACGCCCCCTTCATA', '1'])

    def assertGuideEqual(self, guide, values):
        self.assertEqual(values[0], guide.chromosome)
        self.assertEqual(values[1], guide.start)
        self.assertEqual(values[2], guide.end)
        self.assertEqual(values[3], guide.strand)
        self.assertEqual(values[4], guide.seq)
        self.assertEqual(values[5], guide.uniqueness)
