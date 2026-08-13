import unittest

from generate_transcript_data.cdot_json import _cigar_to_gap_and_length, _add_cds_start_end


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
        cigar = '100=50I100='  # 100 match, 50 insertion (in ref, del in transcript), 100 match = 200 exon length

        _, exon_length = _cigar_to_gap_and_length(cigar)
        self.assertEqual(exon_length, 200)


class AddCdsStartEndTestCase(unittest.TestCase):
    """ UTA data has no genomic cds_start/cds_end so they are derived from the codon
        transcript positions. Expected values are 0-based half-open, matching the GFF parser """

    @staticmethod
    def _transcript_data(strand):
        # Single exon, genomic [1000, 1100), cDNA 1..100, CDS at transcript 0-based [10, 90)
        return {
            "id": "NM_FAKE.1",
            "gene_name": "FAKE",
            "start_codon": 10,
            "stop_codon": 90,
            "genome_builds": {
                "GRCh37": {
                    "contig": "NC_000001.10",
                    "strand": strand,
                    "exons": [[1000, 1100, 0, 1, 100, None]],
                }
            },
        }

    def test_add_cds_start_end_positive_strand(self):
        transcript_data = self._transcript_data("+")
        _add_cds_start_end("GRCh37", transcript_data)
        build_coordinates = transcript_data["genome_builds"]["GRCh37"]
        self.assertEqual(build_coordinates["cds_start"], 1010)
        self.assertEqual(build_coordinates["cds_end"], 1090)

    def test_add_cds_start_end_negative_strand(self):
        """ Both bounds were 1 too high on the minus strand - #95 """
        transcript_data = self._transcript_data("-")
        _add_cds_start_end("GRCh37", transcript_data)
        build_coordinates = transcript_data["genome_builds"]["GRCh37"]
        self.assertEqual(build_coordinates["cds_start"], 1010)
        self.assertEqual(build_coordinates["cds_end"], 1090)

    def test_add_cds_start_end_alignment_gap(self):
        """ A 'D' op (genomic bases missing from the transcript) shifts genomic positions
            past it, so CDS bounds after the gap move by its length """
        transcript_data = self._transcript_data("+")
        build_coordinates = transcript_data["genome_builds"]["GRCh37"]
        # Genomic [1000, 1102), cDNA 1..100, 2 extra genomic bases after the 50th
        build_coordinates["exons"] = [[1000, 1102, 0, 1, 100, "M50 D2 M50"]]
        _add_cds_start_end("GRCh37", transcript_data)
        # start_codon=10 is before the gap (unshifted), stop_codon=90 is past it (+2)
        self.assertEqual(build_coordinates["cds_start"], 1010)
        self.assertEqual(build_coordinates["cds_end"], 1092)


if __name__ == '__main__':
    unittest.main()
