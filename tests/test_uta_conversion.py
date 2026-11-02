import gzip
import json
import os
import shutil
import tempfile
import unittest
from argparse import Namespace
from unittest.mock import patch

from generate_transcript_data.cdot_json import _cigar_to_gap_and_length, _add_cds_start_end, _convert_uta_exons, \
    uta_to_json


class UTAConversionTestCase(unittest.TestCase):
    def test_cigar_to_gap_and_length(self):
        cigar = '194=1D60=1D184='
        expected_gap = 'M194 I1 M60 I1 M184'

        gap, exon_length = _cigar_to_gap_and_length(cigar)
        self.assertEqual(gap, expected_gap)
        # D (deleted from the genome) is transcript-only sequence, so counts towards the exon's cDNA length
        self.assertEqual(exon_length, 440)

    def test_cigar_full_match(self):
        """ Should return None as perfect match """
        cigar = '194='
        expected_gap = None

        gap, exon_length = _cigar_to_gap_and_length(cigar)
        self.assertEqual(gap, expected_gap)
        self.assertEqual(exon_length, 194)

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


class ConvertUTAExonsTestCase(unittest.TestCase):
    def test_transcript_coordinates_follow_cigar(self):
        """ UTA lists exons in transcript order with 0-based genomic coordinates. Exon cDNA ranges are
            1-based and run on from the previous exon's cDNA length (from the CIGAR, not the genomic span) """
        exons = _convert_uta_exons("1000,2000,3000", "1100,2100,3200", "100=,50=2D48=,100=1I99=")
        self.assertEqual(exons, [
            (1000, 1100, 0, 1, 100, None),
            (2000, 2100, 1, 101, 200, "M50 I2 M48"),  # 98 genomic bases, 100 cDNA bases
            (3000, 3200, 2, 201, 399, "M100 D1 M99"),  # 200 genomic bases, 199 cDNA bases
        ])

    def test_minus_strand_sorted_into_genomic_order(self):
        """ UTA minus strand exons come in transcript order (descending genomic), cdot exons are in genomic order """
        exons = _convert_uta_exons("3000,1000", "3100,1100", "100=,100=")
        self.assertEqual(exons, [
            (1000, 1100, 1, 101, 200, None),
            (3000, 3100, 0, 1, 100, None),
        ])


class UTAToJSONTestCase(unittest.TestCase):
    """ uta_to_json on a CSV in the layout the uta_grch38.template query produces """
    URL = "postgresql://uta.biocommons.org/uta_20241220"
    HEADER = "ac,hgnc,origin_url,contig,strand,cds_start_i,cds_end_i,exon_starts,exon_ends,cigars,protein"

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def _run(self, rows) -> dict:
        csv_filename = os.path.join(self.tmp_dir, "uta.csv")
        with open(csv_filename, "w") as f:
            f.write(self.HEADER + "\n")
            for row in rows:
                f.write(",".join(f'"{v}"' if "," in v else v for v in row) + "\n")
        output_filename = os.path.join(self.tmp_dir, "uta.json.gz")
        args = Namespace(uta_csv_filename=csv_filename, url=self.URL, genome_build="GRCh38", output=output_filename)
        with patch("builtins.print") as mock_print:
            uta_to_json(args)
        self.printed = [str(c.args[0]) for c in mock_print.call_args_list]
        with gzip.open(output_filename, "rt") as f:
            return json.load(f)

    def test_uta_to_json(self):
        rows = [
            # + strand, 2 exons, CDS from 2nd base of exon 2 (transcript 0-based 101) to 160
            ["NM_000059.4", "BRCA2", "http://www.ncbi.nlm.nih.gov/refseq/", "NC_000013.11", "1", "101", "160",
             "1000,2000", "1100,2100", "100=,100=", "NP_000050.3,NP_000050.3"],
            # - strand, exons given in transcript order, alignment gap in the first (transcript) exon
            ["NM_001754.5", "RUNX1", "http://www.ncbi.nlm.nih.gov/refseq/", "NC_000021.9", "-1", "10", "160",
             "3000,1000", "3100,1100", "50=1D49=,100=", "NP_001745.2"],
            # Non-coding, no protein
            ["NR_000001.1", "BRCA2", "http://www.ncbi.nlm.nih.gov/refseq/", "NC_000013.11", "1", "", "",
             "5000", "5100", "100=", ""],
            # Aligned to both X and Y (PAR) - skipped
            ["NM_000451.4", "SHOX", "http://www.ncbi.nlm.nih.gov/refseq/", "NC_000023.11,NC_000024.10", "1", "0", "90",
             "600000,700000", "600100,700100", "100=,100=", "NP_000442.1"],
            # Odd accession - skipped
            ["NM_001202404.1/111..1620", "GENE", "http://www.ncbi.nlm.nih.gov/refseq/", "NC_000001.11", "1", "0", "90",
             "100", "200", "100=", ""],
            # CDS end beyond the exons, so the genomic CDS coordinates can't be worked out - skipped
            ["NM_000002.1", "BADCDS", "http://www.ncbi.nlm.nih.gov/refseq/", "NC_000001.11", "1", "10", "500",
             "100", "200", "100=", "NP_000002.1"],
        ]
        data = self._run(rows)

        self.assertEqual(sorted(data["transcripts"]), ["NM_000059.4", "NM_001754.5", "NR_000001.1"])
        self.assertEqual(data["genome_builds"], ["GRCh38"])
        self.assertEqual(data["cdot_version"], data["metadata"]["cdot_version"])
        self.assertEqual(data["metadata"]["method"], "Conversion of Universal Transcript Archive dump")
        self.assertEqual(data["metadata"]["url_counts"], {self.URL: 3})

        brca2 = data["transcripts"]["NM_000059.4"]
        self.assertEqual(brca2["gene_name"], "BRCA2")
        self.assertEqual(brca2["gene_version"], "_BRCA2")
        self.assertEqual(brca2["protein"], "NP_000050.3")
        self.assertEqual((brca2["start_codon"], brca2["stop_codon"]), (101, 160))
        build_data = brca2["genome_builds"]["GRCh38"]
        self.assertEqual(build_data["url"], self.URL)
        self.assertEqual(build_data["contig"], "NC_000013.11")
        self.assertEqual(build_data["strand"], "+")
        self.assertEqual(build_data["exons"], [[1000, 1100, 0, 1, 100, None], [2000, 2100, 1, 101, 200, None]])
        self.assertEqual((build_data["cds_start"], build_data["cds_end"]), (2001, 2060))

        runx1 = data["transcripts"]["NM_001754.5"]
        build_data = runx1["genome_builds"]["GRCh38"]
        self.assertEqual(build_data["strand"], "-")
        self.assertEqual(build_data["exons"], [[1000, 1100, 1, 101, 200, None], [3000, 3100, 0, 1, 100, "M50 I1 M49"]])
        # start_codon 10 is 10 bases into the first transcript exon from its genomic end (3099 - 10 = 3089),
        # the last CDS base (transcript 159) is 59 bases into the second from its end (1099 - 59 = 1040)
        self.assertEqual((build_data["cds_start"], build_data["cds_end"]), (1040, 3090))

        ncrna = data["transcripts"]["NR_000001.1"]
        for key in ["start_codon", "stop_codon", "protein"]:
            self.assertNotIn(key, ncrna)
        self.assertNotIn("cds_start", ncrna["genome_builds"]["GRCh38"])

        # Fake genes (not real accessions, so merge_historical swaps them for real ones) list their transcripts
        self.assertEqual(data["genes"], {
            "_BRCA2": {"gene_symbol": "BRCA2", "transcripts": ["NM_000059.4", "NR_000001.1"], "url": self.URL},
            "_RUNX1": {"gene_symbol": "RUNX1", "transcripts": ["NM_001754.5"], "url": self.URL},
            "_BADCDS": {"gene_symbol": "BADCDS", "transcripts": ["NM_000002.1"], "url": self.URL},
        })
        self.assertTrue(any(p.startswith("Warning: Skipped NM_000002.1") for p in self.printed), self.printed)


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
