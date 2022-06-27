import unittest
from generate_transcript_data.cdot_json import _cigar_to_gap_and_length


class UTAConversionTestCase(unittest.TestCase):
    def test_cigar_to_gap_and_length(self):
        cigar = '194=1D60=1D184='
        expected_gap = 'M194 I1 M60 I1 M184'

        gap, exon_length = _cigar_to_gap_and_length(cigar)
        self.assertEqual(gap, expected_gap)

    def test_cigar_full_match(self):
        """ Should return None as perfect match """
        cigar = '194='
        expected_gap = None

        gap, exon_length = _cigar_to_gap_and_length(cigar)
        self.assertEqual(gap, expected_gap)

    def test_cigar_merged_matches(self):
        cigar = '194=100='
        expected_gap = None

        gap, exon_length = _cigar_to_gap_and_length(cigar)
        self.assertEqual(gap, expected_gap)

    def test_cigar_mismatch(self):
        cigar = '195=1X1D430='  # X will become match and should merge w/first
        expected_gap = "M196 I1 M430"

        gap, exon_length = _cigar_to_gap_and_length(cigar)
        self.assertEqual(gap, expected_gap)

    def test_cigar_deletion_exon_length(self):
        cigar = '100=50D100='  # 100 match, 50 deletion, 100 match = 200 exon length

        _, exon_length = _cigar_to_gap_and_length(cigar)
        self.assertEqual(exon_length, 200)


if __name__ == '__main__':
    unittest.main()
