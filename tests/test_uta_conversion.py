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

    def test_cigar_mismatch(self):
        cigar = '195=1X430='
        expected_gap = "M195 M1 M430"

        gap, exon_length = _cigar_to_gap_and_length(cigar)
        self.assertEqual(gap, expected_gap)


if __name__ == '__main__':
    unittest.main()
