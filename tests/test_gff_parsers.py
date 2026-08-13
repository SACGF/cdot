
import gzip
import os
import tempfile
from inspect import getsourcefile
import unittest
from generate_transcript_data.cdot_json import add_gencode_hgnc
from generate_transcript_data.gff_parser import GTFParser, GFF3Parser


class Test(unittest.TestCase):
    this_file_dir = os.path.dirname(os.path.abspath(getsourcefile(lambda: 0)))
    test_data_dir = os.path.join(this_file_dir, "test_data")
    ENSEMBL_104_GTF_FILENAME = os.path.join(test_data_dir, "ensembl_test.GRCh38.104.gtf")
    ENSEMBL_111_GTF_FILENAME = os.path.join(test_data_dir, "ensembl_test.GRCh38.111.gtf")
    # Older RefSeq, before Genbank => GenBank changed
    REFSEQ_GFF3_FILENAME_2021 = os.path.join(test_data_dir, "refseq_test.GRCh38.p13_genomic.109.20210514.gff")
    # Newer RefSeq, before Genbank => GenBank changed
    REFSEQ_GFF3_FILENAME_2023 = os.path.join(test_data_dir, "refseq_test.GRCh38.p14_genomic.RS_2023_03.gff")
    # Synthetic - alignment leaves a hole in the transcript coordinates (issue #123)
    REFSEQ_GFF3_FILENAME_COORDINATE_HOLE = os.path.join(test_data_dir,
                                                        "refseq_test.transcript_coordinate_hole.gff")
    # NCBI historical alignments (issue #51) - annotation + alignments files concatenated (annotation first).
    # Contains ACADM NM_000016.2/.3, gapped NM_000028.1 and NM_000066.1, partial-start NM_002521.1,
    # plus an alignment-only NM_001005277.1 (annotation deliberately omitted) to exercise skip_missing_parents
    REFSEQ_GFF3_FILENAME_HISTORICAL = os.path.join(test_data_dir, "refseq_test.historical_RS_2024_08.gff")
    REFSEQ_GFF3_FILENAME_GRCH37_MT = os.path.join(test_data_dir, "refseq_grch37_mt.gff")
    REFSEQ_GFF3_FILENAME_GRCH38_MT = os.path.join(test_data_dir, "refseq_grch38.p14_mt.gff")
    UCSC_GTF_FILENAME = os.path.join(test_data_dir, "hg19_chrY_300kb_genes.gtf")
    FAKE_URL = "http://fake.url"

    FAKE_MT_TRANSCRIPTS = [
        "fake-rna-ATP6", "fake-rna-ATP8", "fake-rna-COX1", "fake-rna-COX2", "fake-rna-COX3", "fake-rna-CYTB",
        "fake-rna-ND1", "fake-rna-ND2", "fake-rna-ND3", "fake-rna-ND4", "fake-rna-ND4L", "fake-rna-ND5", "fake-rna-ND6"
    ]

    def _test_exon_length(self, transcripts, genome_build, transcript_id, expected_length):
        transcript = transcripts[transcript_id]
        exons = transcript["genome_builds"][genome_build]["exons"]
        length = sum([exon[1] - exon[0] for exon in exons])
        self.assertEqual(expected_length, length, "%s exons sum" % transcript_id)

    def test_ucsc_gtf(self):
        genome_build = "GRCh37"
        parser = GTFParser(self.UCSC_GTF_FILENAME, genome_build, self.FAKE_URL)
        _, transcripts = parser.get_genes_and_transcripts()
        self._test_exon_length(transcripts, genome_build, "NM_013239", 2426)

    def test_ensembl_gtf(self):
        genome_build = "GRCh38"
        parser = GTFParser(self.ENSEMBL_104_GTF_FILENAME, genome_build, self.FAKE_URL)
        genes, transcripts = parser.get_genes_and_transcripts()
        self._test_exon_length(transcripts, genome_build, "ENST00000357654.9", 7088)

        # Ensure that geneID was inserted with a version
        expected_gene_version = "ENSG00000012048.23"

        transcript = transcripts["ENST00000357654.9"]
        transcript_gene_version = transcript["gene_version"]
        self.assertEqual(expected_gene_version, transcript_gene_version, "Transcript gene has version")

        self.assertTrue(expected_gene_version in genes, f"{expected_gene_version=} in genes")

        protein = transcript.get("protein")
        self.assertEqual(protein, "ENSP00000350283.3")

    def test_refseq_gff3_2021(self):
        genome_build = "GRCh38"
        parser = GFF3Parser(self.REFSEQ_GFF3_FILENAME_2021, genome_build, self.FAKE_URL)
        _, transcripts = parser.get_genes_and_transcripts()
        self._test_exon_length(transcripts, genome_build, "NM_007294.4", 7088)

        transcript = transcripts["NM_015120.4"]
        protein = transcript.get("protein")
        self.assertEqual(protein, "NP_055935.4")

    def test_refseq_gff3_2023(self):
        genome_build = "GRCh38"
        parser = GFF3Parser(self.REFSEQ_GFF3_FILENAME_2023, genome_build, self.FAKE_URL)
        _, transcripts = parser.get_genes_and_transcripts()
        self._test_exon_length(transcripts, genome_build, "NM_007294.4", 7088)

        transcript = transcripts["NM_015120.4"]
        protein = transcript.get("protein")
        self.assertEqual(protein, "NP_055935.4")

    def test_refseq_gff3_historical(self):
        """ NCBI historical transcript alignments, issue #51 """
        genome_build = "GRCh38"
        parser = GFF3Parser(self.REFSEQ_GFF3_FILENAME_HISTORICAL, genome_build, self.FAKE_URL,
                            skip_missing_parents=True)
        _, transcripts = parser.get_genes_and_transcripts()

        # The alignment-only transcript is skipped, everything else is kept
        self.assertEqual(sorted(transcripts),
                         ["NM_000016.2", "NM_000016.3", "NM_000028.1", "NM_000066.1", "NM_002521.1"])
        self.assertEqual(parser.skipped_features_no_parents["cDNA_match"], 1)

        # Two historical versions of the same transcript, each with its own alignment
        # (in 2023 a data-generation bug meant every exon appeared twice - guard against that)
        self._test_exon_length(transcripts, genome_build, "NM_000016.2", 2192)
        self._test_exon_length(transcripts, genome_build, "NM_000016.3", 2423)
        for accession in ["NM_000016.2", "NM_000016.3"]:
            exons = transcripts[accession]["genome_builds"][genome_build]["exons"]
            exon_coords = [(e[0], e[1]) for e in exons]
            self.assertEqual(len(exon_coords), len(set(exon_coords)), f"{accession} has duplicated exons")

        # Alignment gaps come through from the cDNA_match Gap attribute
        agl = transcripts["NM_000028.1"]["genome_builds"][genome_build]
        gaps = [e[5] for e in agl["exons"] if e[5]]
        self.assertEqual(gaps, ["M148 D1 M35 D2 M371 I2 M1135 D1 M794"])

        # CDS from the annotation file combined with alignment coordinates
        acadm = transcripts["NM_000016.3"]
        self.assertEqual(acadm["gene_name"], "ACADM")
        self.assertEqual(acadm["start_codon"], 430)
        self.assertEqual(acadm["stop_codon"], 1696)

        # Partial alignment: the transcript's first base doesn't align, cDNA coordinates start at 2
        nppb = transcripts["NM_002521.1"]["genome_builds"][genome_build]
        last_exon = nppb["exons"][-1]  # - strand, so first transcript exon is last in genomic order
        self.assertEqual(last_exon[3], 2, "NM_002521.1 alignment starts at base 2 of the transcript")

        # Without skip_missing_parents the alignment-only transcript is an error
        parser = GFF3Parser(self.REFSEQ_GFF3_FILENAME_HISTORICAL, genome_build, self.FAKE_URL)
        with self.assertRaises(ValueError):
            parser.get_genes_and_transcripts()

    def test_exons_in_genomic_order(self):
        genome_build = "GRCh38"
        parser = GTFParser(self.ENSEMBL_104_GTF_FILENAME, genome_build, self.FAKE_URL)
        _, transcripts = parser.get_genes_and_transcripts()
        transcript = transcripts["ENST00000357654.9"]
        exons = transcript["genome_builds"][genome_build]["exons"]
        first_exon = exons[0]
        last_exon = exons[-1]
        self.assertGreater(last_exon[0], first_exon[0])

        parser = GFF3Parser(self.REFSEQ_GFF3_FILENAME_2021, genome_build, self.FAKE_URL)
        _, transcripts = parser.get_genes_and_transcripts()
        transcript = transcripts["NM_007294.4"]
        self.assertEqual(transcript.get("hgnc"), "1100", f"{transcript} has HGNC:1100")
        exons = transcript["genome_builds"][genome_build]["exons"]
        first_exon = exons[0]
        last_exon = exons[-1]
        self.assertGreater(last_exon[0], first_exon[0])

        parser = GFF3Parser(self.REFSEQ_GFF3_FILENAME_2023, genome_build, self.FAKE_URL)
        _, transcripts = parser.get_genes_and_transcripts()
        transcript = transcripts["NM_007294.4"]
        self.assertEqual(transcript.get("hgnc"), "1100", f"{transcript} has HGNC:1100")
        exons = transcript["genome_builds"][genome_build]["exons"]
        first_exon = exons[0]
        last_exon = exons[-1]
        self.assertGreater(last_exon[0], first_exon[0])

    def test_ensembl_gtf_tags(self):
        genome_build = "GRCh38"
        parser = GTFParser(self.ENSEMBL_111_GTF_FILENAME, genome_build, self.FAKE_URL)
        genes, transcripts = parser.get_genes_and_transcripts()
        transcript = transcripts["ENST00000641515.2"]
        tag = transcript["genome_builds"][genome_build].get("tag")
        self.assertIn("MANE_Select", tag)

    def test_chrom_contig_conversion(self):
        genome_build = "GRCh38"
        parser = GTFParser(self.ENSEMBL_111_GTF_FILENAME, genome_build, self.FAKE_URL)
        _, transcripts = parser.get_genes_and_transcripts()
        transcript = transcripts["ENST00000641515.2"]
        contig = transcript["genome_builds"][genome_build].get("contig")
        self.assertEqual(contig, "NC_000001.11")

    def test_ncrna_gene(self):
        """ We were incorrectly missing ncRNA gene info @see https://github.com/SACGF/cdot/issues/72 """
        genome_build = "GRCh38"
        parser = GTFParser(self.ENSEMBL_111_GTF_FILENAME, genome_build, self.FAKE_URL)
        genes, transcripts = parser.get_genes_and_transcripts()
        gene = genes["ENSG00000210156"]
        gene_symbol = gene["gene_symbol"]
        self.assertEqual(gene_symbol, "MT-TK")

    def _test_mito(self, filename, genome_build):
        parser = GFF3Parser(filename, genome_build, self.FAKE_URL)
        genes, transcripts = parser.get_genes_and_transcripts()

        for transcript_accession in self.FAKE_MT_TRANSCRIPTS:
            self.assertIn(transcript_accession, transcripts)

        transcript = transcripts["fake-rna-ATP6"]
        exons = transcript["genome_builds"][genome_build]["exons"]
        first_exon = exons[0]
        self.assertEqual(first_exon[0], 8526)
        self.assertEqual(first_exon[1], 9207)

    def test_mito_mrna(self):
        """ Need to make fake MT transcripts for RefSeq @see https://github.com/SACGF/cdot/issues/72 """
        self._test_mito(self.REFSEQ_GFF3_FILENAME_GRCH38_MT, "GRCh38")

    def test_mito_no_mrna(self):
        """ Need to make fake MT transcripts for RefSeq @see https://github.com/SACGF/cdot/issues/72 """
        self._test_mito(self.REFSEQ_GFF3_FILENAME_GRCH37_MT, "GRCh37")

    def test_transcript_position_across_coordinate_hole(self):
        """ A few RefSeq alignments leave a run of transcript bases unaligned between two exons, so
            the exon transcript coordinates have a hole in them. Codon positions must stay in whole
            transcript coordinates rather than collapsing the hole out.
            @see https://github.com/SACGF/cdot/issues/123 """
        # (alt_start, alt_end, exon_id, cds_start, cds_end, gap) - stranded order.
        # Exon 1 ends at transcript 200, exon 2 picks up at 231, so 30 bases align nowhere.
        exons = [
            (1_000, 1_100, 0, 1, 100, None),
            (2_000, 2_100, 1, 101, 200, None),
            (3_000, 3_100, 2, 231, 330, None),
        ]
        # 50 bases into the exon that follows the hole
        self.assertEqual(GFF3Parser._get_transcript_position(True, exons, 3_050), 280)
        # Exons before the hole are unaffected
        self.assertEqual(GFF3Parser._get_transcript_position(True, exons, 2_050), 150)

        # Same again on the minus strand (exons in stranded order, so genomic order is reversed)
        rev_exons = [
            (3_000, 3_100, 0, 1, 100, None),
            (2_000, 2_100, 1, 101, 200, None),
            (1_000, 1_100, 2, 231, 330, None),
        ]
        self.assertEqual(GFF3Parser._get_transcript_position(False, rev_exons, 1_050), 280)
        self.assertEqual(GFF3Parser._get_transcript_position(False, rev_exons, 2_050), 150)

    def test_codon_positions_across_coordinate_hole(self):
        """ End to end version of the above: the codon positions the parser writes out must be in
            the same coordinate system as the exon cds_start/cds_end.
            @see https://github.com/SACGF/cdot/issues/123 """
        genome_build = "GRCh38"
        parser = GFF3Parser(self.REFSEQ_GFF3_FILENAME_COORDINATE_HOLE, genome_build, self.FAKE_URL)
        _, transcripts = parser.get_genes_and_transcripts()
        transcript = transcripts["NM_000123.1"]

        exons = transcript["genome_builds"][genome_build]["exons"]
        self.assertEqual([tuple(e[3:5]) for e in exons], [(1, 100), (101, 200), (231, 330)])

        # CDS runs from genomic 1010 (10 bases into exon 1) to 3050 (50 bases into exon 3)
        self.assertEqual(transcript["start_codon"], 10)
        self.assertEqual(transcript["stop_codon"], 280)

    def test_ensembl_hgnc_injection(self):
        """ Test that GENCODE HGNC metadata is properly injected into Ensembl GTF transcripts and genes.
            @see https://github.com/SACGF/cdot/issues/97 """
        genome_build = "GRCh38"
        parser = GTFParser(self.ENSEMBL_111_GTF_FILENAME, genome_build, self.FAKE_URL)
        genes, transcripts = parser.get_genes_and_transcripts()

        # Verify no HGNC present before injection
        self.assertIsNone(transcripts["ENST00000641515.2"].get("hgnc"))
        self.assertIsNone(transcripts["ENST00000387421.1"].get("hgnc"))

        # Create a minimal GENCODE HGNC metadata file matching the two transcripts in the test GTF
        # OR4F5 = HGNC:14825, MT-TK = HGNC:7489 (confirmed from gencode.v45.metadata.HGNC.gz)
        hgnc_content = (
            "ENST00000641515.2\tOR4F5\tHGNC:14825\n"
            "ENST00000387421.1\tMT-TK\tHGNC:7489\n"
        )
        with tempfile.NamedTemporaryFile(suffix=".gz", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            with gzip.open(tmp_path, "wt") as f:
                f.write(hgnc_content)
            add_gencode_hgnc(tmp_path, genes, transcripts)
        finally:
            os.unlink(tmp_path)

        # HGNC injected into transcripts
        self.assertEqual(transcripts["ENST00000641515.2"]["hgnc"], "14825")
        self.assertEqual(transcripts["ENST00000387421.1"]["hgnc"], "7489")

        # HGNC injected into the gene entries (looked up via each transcript's gene_version key)
        gene_version_or4f5 = transcripts["ENST00000641515.2"]["gene_version"]
        gene_version_mttk = transcripts["ENST00000387421.1"]["gene_version"]
        self.assertEqual(genes[gene_version_or4f5]["hgnc"], "14825")
        self.assertEqual(genes[gene_version_mttk]["hgnc"], "7489")
