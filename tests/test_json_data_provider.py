import unittest

import hgvs
from hgvs.assemblymapper import AssemblyMapper

from seedot.hgvs.dataproviders.json_data_provider import JSONDataProvider


class TestJSONDataProvider(unittest.TestCase):
    # TODO: Test loading from file, filename

    def setUp(self):
        GRCH37_TRANSCRIPTS = {
            "NM_001637.3":
                {'chrom': 'NC_000007.13',
                 'start': 36552548,
                 'end': 36764154,
                 'strand': '-',
                 'cds_start': 36552857,
                 'cds_end': 36763753,
                 'exons': [[36763626, 36764154],
                           [36726303, 36726399],
                           [36713547, 36713614],
                           [36698770, 36698870],
                           [36677456, 36677516],
                           [36671641, 36671712],
                           [36662795, 36662856],
                           [36661315, 36661386],
                           [36660386, 36660435],
                           [36657902, 36657951],
                           [36655985, 36656080],
                           [36633944, 36634036],
                           [36616179, 36616262],
                           [36589044, 36589081],
                           [36588217, 36588292],
                           [36579924, 36580097],
                           [36571891, 36571950],
                           [36571752, 36571812],
                           [36570023, 36570120],
                           [36561644, 36561721],
                           [36552548, 36552986]],
                 'start_codon_transcript_pos': 401,
                 'stop_codon_transcript_pos': 2129,
                 'cdna_match': [[36763626, 36764154, 1, 528, None],
                                [36726303, 36726399, 529, 624, None],
                                [36713547, 36713614, 625, 691, None],
                                [36698770, 36698870, 692, 791, None],
                                [36677456, 36677516, 792, 851, None],
                                [36671641, 36671712, 852, 922, None],
                                [36662795, 36662856, 923, 983, None],
                                [36661315, 36661386, 984, 1054, None],
                                [36660386, 36660435, 1055, 1103, None],
                                [36657902, 36657951, 1104, 1152, None],
                                [36655985, 36656080, 1153, 1247, None],
                                [36633944, 36634036, 1248, 1339, None],
                                [36616179, 36616262, 1340, 1422, None],
                                [36589044, 36589081, 1423, 1459, None],
                                [36588217, 36588292, 1460, 1534, None],
                                [36579924, 36580097, 1535, 1707, None],
                                [36571891, 36571950, 1708, 1766, None],
                                [36571752, 36571812, 1767, 1826, None],
                                [36570023, 36570120, 1827, 1923, None],
                                [36561644, 36561721, 1924, 2000, None],
                                [36552548, 36552986, 2001, 2440, 'M196 I1 M61 I1 M181']],
                 'biotype': 'protein_coding',
                 'id': 'NM_001637.3',
                 'gene_version': '313',
                 'gene_name': 'AOAH',
                 'url': 'http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/ref_GRCh37.p13_top_level.gff3.gz'}
        }
        self.json_data_provider = JSONDataProvider(json_data=GRCH37_TRANSCRIPTS)

    def test_hgvs_conversion(self):
        am = AssemblyMapper(self.json_data_provider,
                            assembly_name='GRCh37', alt_aln_method='splign', replace_reference=True)
        HGVS_C_TO_G = [
            ('NM_001637.3:c.1582G>A', 'NC_000007.13:g.36561662C>T'),
        ]

        hp = hgvs.parser.Parser()
        for hgvs_c, expected_hgvs_g in HGVS_C_TO_G:
            var_c = hp.parse_hgvs_variant(hgvs_c)
            var_g = am.c_to_g(var_c)
            self.assertEqual(str(var_g), expected_hgvs_g)


if __name__ == '__main__':
    unittest.main()
