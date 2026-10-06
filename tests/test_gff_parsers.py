
import argparse
import gzip
import json
import os
import tempfile
from inspect import getsourcefile
import unittest
from generate_transcript_data.cdot_json import add_gencode_hgnc, combine_builds, write_cdot_json
from generate_transcript_data.annotation_consortium import RefSeq
from generate_transcript_data.gff_parser import GTFParser, GFF3Parser
from generate_transcript_data.transcript_builder import TranscriptBuilder
from generate_transcript_data.transcript_coordinates import get_transcript_position


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
    # SELENOM (- strand, selenocysteine codon in exon 2) from RefSeq 110 and Ensembl 108
    REFSEQ_GFF3_FILENAME_SELENOPROTEIN = os.path.join(test_data_dir, "refseq_test.selenoprotein.gff")
    ENSEMBL_108_GTF_FILENAME_SELENOPROTEIN = os.path.join(test_data_dir,
                                                          "ensembl_test.GRCh38.108.selenoprotein.gtf")
    # GPX4 (+ strand) and SELENOP (- strand, 10 selenocysteines) from RefSeq RS_2024_08
    REFSEQ_GFF3_FILENAME_SELENOPROTEIN_RS_2024_08 = os.path.join(test_data_dir,
                                                                 "refseq_test.RS_2024_08.selenoprotein.gff")
    # GPX4 (+ strand) and SELENOH ENST00000528798 (cds_start_NF, CDS out of frame) from Ensembl 115
    ENSEMBL_115_GTF_FILENAME_SELENOPROTEIN = os.path.join(test_data_dir,
                                                          "ensembl_test.GRCh38.115.selenoprotein.gtf")
    # PEG10 (-1 frameshift), OAZ1 (+1) and OAZ2 (+1, - strand) from RefSeq RS_2025_08
    REFSEQ_GFF3_FILENAME_RIBOSOMAL_SLIPPAGE = os.path.join(test_data_dir,
                                                           "refseq_test.RS_2025_08.ribosomal_slippage.gff")
    # ACTN3: the GRCh37 reference has the R577X stop allele, which RefSeq corrects on GRCh37 only
    REFSEQ_GFF3_FILENAME_GRCH37_ACTN3 = os.path.join(test_data_dir, "refseq_test.GRCh37.105.20220307.ACTN3.gff")
    REFSEQ_GFF3_FILENAME_GRCH38_ACTN3 = os.path.join(test_data_dir, "refseq_test.RS_2025_08.ACTN3.gff")
    UCSC_GTF_FILENAME = os.path.join(test_data_dir, "hg19_chrY_300kb_genes.gtf")
    # Ensembl GFF3 (chrMT up to the first CDS). CDS rows have the protein version from release 114
    ENSEMBL_115_GFF3_FILENAME = os.path.join(test_data_dir, "ensembl_test.GRCh38.115.MT.gff3")
    # MT-ND1 (stop codon completed by the poly(A) tail), MT-TI (tRNA), MT-ATP8 and MT-ND6 (- strand)
    ENSEMBL_115_GTF_FILENAME_MT = os.path.join(test_data_dir, "ensembl_test.GRCh38.115.MT.gtf")
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

    def test_ensembl_gff3(self):
        """ Ensembl GFF3 CDS rows have the protein version in 'version' from release 114
            @see https://github.com/SACGF/cdot/issues/135 """
        parser = GFF3Parser(self.ENSEMBL_115_GFF3_FILENAME, "GRCh38", self.FAKE_URL)
        genes, transcripts = parser.get_genes_and_transcripts()
        self.assertEqual(transcripts["ENST00000361390.2"].get("protein"), "ENSP00000354687.2")
        # Gene version is read off ncRNA_gene rows too, not just 'gene'
        self.assertEqual(transcripts["ENST00000387365.1"]["gene_version"], "ENSG00000210100.1")
        self.assertIn("ENSG00000210100.1", genes)

    def test_ensembl_gff3_before_114_not_supported(self):
        """ Ensembl GFF3 CDS rows before release 114 have protein_id but no version.
            The error must say so rather than just 'missing version' @see https://github.com/SACGF/cdot/issues/101 """
        with open(self.ENSEMBL_115_GFF3_FILENAME) as f:
            gff3 = f.read().replace("protein_id=ENSP00000354687;version=2", "protein_id=ENSP00000354687")
        with tempfile.NamedTemporaryFile("w", suffix=".gff3") as f:
            f.write(gff3)
            f.flush()
            parser = GFF3Parser(f.name, "GRCh38", self.FAKE_URL)
            with self.assertRaisesRegex(ValueError, "Ensembl GFF3 files before release 114"):
                parser.get_genes_and_transcripts()

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
        gene = genes["ENSG00000210156.1"]
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

    def test_mito_transl_except_and_transl_table(self):
        """ RefSeq MT CDS rows carry transl_table=2, and transl_except=(...,aa:TERM) where the stop
            codon is completed by the poly(A) tail """
        genome_build = "GRCh38"
        parser = GFF3Parser(self.REFSEQ_GFF3_FILENAME_GRCH38_MT, genome_build, self.FAKE_URL)
        _, transcripts = parser.get_genes_and_transcripts()
        for transcript_accession in self.FAKE_MT_TRANSCRIPTS:
            self.assertEqual(transcripts[transcript_accession]["translation"]["transl_table"], 2, transcript_accession)

        # ND1: CDS of 956 bases, the last 2 bases (TA) of the stop codon are in the genome
        nd1 = transcripts["fake-rna-ND1"]
        self.assertEqual(nd1["stop_codon"] - nd1["start_codon"], 956)
        self.assertEqual(nd1["translation"]["transl_except"], {"TERM": [319]})
        # ATP8 has a complete stop codon
        self.assertNotIn("transl_except", transcripts["fake-rna-ATP8"]["translation"])

    def test_ensembl_gtf_mito_transl_table(self):
        """ Ensembl GTFs don't name the genetic code, so coding transcripts on MT get the vertebrate
            mitochondrial code, as RefSeq writes on its MT CDS rows """
        genome_build = "GRCh38"
        parser = GTFParser(self.ENSEMBL_115_GTF_FILENAME_MT, genome_build, self.FAKE_URL)
        _, transcripts = parser.get_genes_and_transcripts()
        for transcript_accession in ["ENST00000361390.2", "ENST00000361851.1", "ENST00000361681.2"]:
            transcript = transcripts[transcript_accession]
            self.assertEqual(transcript["genome_builds"][genome_build]["contig"], "NC_012920.1")
            self.assertEqual(transcript["translation"], {"transl_table": 2}, transcript_accession)
        self.assertNotIn("translation", transcripts["ENST00000387365.1"])  # MT-TI, a tRNA

    def test_refseq_gff3_selenocysteine(self):
        """ RefSeq names the selenocysteine codon in transl_except on every CDS row """
        genome_build = "GRCh38"
        parser = GFF3Parser(self.REFSEQ_GFF3_FILENAME_SELENOPROTEIN, genome_build, self.FAKE_URL)
        _, transcripts = parser.get_genes_and_transcripts()
        transcript = transcripts["NM_080430.4"]
        self.assertEqual(transcript["translation"], {"transl_except": {"Sec": [48]}})  # no transl_table on nuclear CDS

    def test_ensembl_gtf_selenocysteine(self):
        """ Ensembl GTF writes the selenocysteine codon as a 'Selenocysteine' row """
        genome_build = "GRCh38"
        parser = GTFParser(self.ENSEMBL_108_GTF_FILENAME_SELENOPROTEIN, genome_build, self.FAKE_URL)
        _, transcripts = parser.get_genes_and_transcripts()
        transcript = transcripts["ENST00000400299.6"]
        self.assertEqual(transcript["translation"], {"transl_except": {"Sec": [48]}})
        # The Selenocysteine row must not change the transcript, exons or CDS
        exons = transcript["genome_builds"][genome_build]["exons"]
        self.assertEqual(len(exons), 5)
        self.assertEqual(transcript["stop_codon"] - transcript["start_codon"], 438)
        self.assertNotIn("warnings", transcript["genome_builds"][genome_build])

    def test_refseq_gff3_selenocysteine_plus_strand_and_multiple(self):
        genome_build = "GRCh38"
        parser = GFF3Parser(self.REFSEQ_GFF3_FILENAME_SELENOPROTEIN_RS_2024_08, genome_build, self.FAKE_URL)
        _, transcripts = parser.get_genes_and_transcripts()
        # GPX4 U73
        gpx4 = transcripts["NM_002085.5"]
        self.assertEqual(gpx4["genome_builds"][genome_build]["strand"], "+")
        self.assertEqual(gpx4["translation"]["transl_except"], {"Sec": [73]})
        # SELENOP, as in UniProt P49908
        selenop = transcripts["NM_005410.4"]
        self.assertEqual(selenop["translation"]["transl_except"], {"Sec": [59, 300, 318, 330, 345, 352, 367, 369, 376, 378]})
        for transcript in (gpx4, selenop):
            self.assertNotIn("warnings", transcript["genome_builds"][genome_build])

    def test_ensembl_gtf_selenocysteine_plus_strand(self):
        genome_build = "GRCh38"
        parser = GTFParser(self.ENSEMBL_115_GTF_FILENAME_SELENOPROTEIN, genome_build, self.FAKE_URL)
        _, transcripts = parser.get_genes_and_transcripts()
        gpx4 = transcripts["ENST00000354171.13"]
        self.assertEqual(gpx4["translation"]["transl_except"], {"Sec": [73]})
        self.assertNotIn("warnings", gpx4["genome_builds"][genome_build])

    def test_ensembl_gtf_selenocysteine_unplaced(self):
        """ SELENOH ENST00000528798 is cds_start_NF, and its CDS starts out of frame, so the
            selenocysteine can't be given a codon number. It must be flagged, not silently dropped """
        genome_build = "GRCh38"
        parser = GTFParser(self.ENSEMBL_115_GTF_FILENAME_SELENOPROTEIN, genome_build, self.FAKE_URL)
        with self.assertLogs(level="WARNING"):
            _, transcripts = parser.get_genes_and_transcripts()
        transcript = transcripts["ENST00000528798.1"]
        self.assertNotIn("translation", transcript)
        self.assertEqual(transcript["genome_builds"][genome_build]["warnings"], {"transl_except_unplaced": ["Sec"]})

    def test_ensembl_gtf_selenocysteine_split_across_exons(self):
        """ A codon split across exons would be a Selenocysteine row per exon (not seen in real data yet).
            GPX4 codon 60 spans exons 2/3, the 2nd piece is mid-codon but must not be flagged unplaced """
        with open(self.ENSEMBL_115_GTF_FILENAME_SELENOPROTEIN) as f:
            lines = f.readlines()
        sec_line = next(line for line in lines if "\tSelenocysteine\t1105403\t1105405\t" in line)
        split_lines = [sec_line.replace("\t1105403\t1105405\t", "\t1105279\t1105280\t"),
                       sec_line.replace("\t1105403\t1105405\t", "\t1105366\t1105366\t")]
        lines = [line for line in lines if line is not sec_line] + split_lines

        genome_build = "GRCh38"
        with tempfile.NamedTemporaryFile("w", suffix=".gtf") as f:
            f.writelines(lines)
            f.flush()
            parser = GTFParser(f.name, genome_build, self.FAKE_URL)
            _, transcripts = parser.get_genes_and_transcripts()
        gpx4 = transcripts["ENST00000354171.13"]
        self.assertEqual(gpx4["translation"]["transl_except"], {"Sec": [60]})
        self.assertNotIn("warnings", gpx4["genome_builds"][genome_build])

    def test_refseq_gff3_ribosomal_slippage(self):
        """ RefSeq splits the CDS at a ribosomal frameshift: rows overlapping by a base is -1 (base read
            twice), a 1 base gap is +1 (base skipped) @see https://github.com/SACGF/cdot/issues/76 """
        genome_build = "GRCh38"
        parser = GFF3Parser(self.REFSEQ_GFF3_FILENAME_RIBOSOMAL_SLIPPAGE, genome_build, self.FAKE_URL)
        _, transcripts = parser.get_genes_and_transcripts()
        # (cds_position, shift, protein length incl stop codon)
        expected = {
            "NM_015068.3": (957, -1, 709),  # PEG10
            "NM_004152.3": (205, 1, 229),  # OAZ1, skips the U of UGA at codon 69
            "NM_002537.3": (97, 1, 190),  # OAZ2, - strand
        }
        for transcript_accession, (cds_position, shift, num_codons) in expected.items():
            transcript = transcripts[transcript_accession]
            self.assertEqual(transcript["translation"]["ribosomal_slippage"],
                             [{"cds_position": cds_position, "shift": shift}],
                             transcript_accession)
            self.assertEqual(transcript["translation"]["exceptions"], ["ribosomal slippage"])
            self.assertNotIn("warnings", transcript["genome_builds"][genome_build])
            # Applying the shift gives a whole number of codons
            translated_length = transcript["stop_codon"] - transcript["start_codon"] - shift
            self.assertEqual(translated_length, num_codons * 3, transcript_accession)

    def test_refseq_gff3_ribosomal_slippage_unplaced(self):
        """ A slippage exception where the CDS rows don't show where must be flagged """
        with open(self.REFSEQ_GFF3_FILENAME_RIBOSOMAL_SLIPPAGE) as f:
            lines = f.readlines()
        # Move PEG10's 2nd CDS row off the overlap
        lines = [line.replace("\t94664513\t94665682\t", "\t94664520\t94665682\t") for line in lines]

        genome_build = "GRCh38"
        with tempfile.NamedTemporaryFile("w", suffix=".gff") as f:
            f.writelines(lines)
            f.flush()
            parser = GFF3Parser(f.name, genome_build, self.FAKE_URL)
            with self.assertLogs(level="WARNING"):
                _, transcripts = parser.get_genes_and_transcripts()
        transcript = transcripts["NM_015068.3"]
        self.assertNotIn("ribosomal_slippage", transcript["translation"])
        self.assertEqual(transcript["genome_builds"][genome_build]["warnings"], {"ribosomal_slippage_unplaced": True})

    def test_refseq_genome_mismatch_per_build(self):
        """ RefSeq transl_except that correct a genome codon depend on the build, so go in genome_builds,
            and must stay with their build when builds are combined """
        transcript_accession = "NM_001104.4"
        builds = {"GRCh37": self.REFSEQ_GFF3_FILENAME_GRCH37_ACTN3, "GRCh38": self.REFSEQ_GFF3_FILENAME_GRCH38_ACTN3}
        with tempfile.TemporaryDirectory() as temp_dir:
            args = {"output": os.path.join(temp_dir, "combined.json.gz")}
            for genome_build, filename in builds.items():
                _, transcripts = GFF3Parser(filename, genome_build, self.FAKE_URL).get_genes_and_transcripts()
                transcript = transcripts[transcript_accession]
                self.assertNotIn("translation", transcript)
                build_data = transcript["genome_builds"][genome_build]
                if genome_build == "GRCh37":
                    self.assertEqual(build_data["genome_mismatch"],
                                     {"transl_except": {"Arg": [577]},
                                      "exceptions": ["annotated by transcript or proteomic data"],
                                      "transcript": {"substitutions": 1},
                                      "protein": {"substitutions": 1}})
                else:
                    self.assertNotIn("genome_mismatch", build_data)
                args[genome_build.lower()] = os.path.join(temp_dir, f"{genome_build}.json.gz")
                write_cdot_json(args[genome_build.lower()], "test", [], {}, transcripts, [genome_build])
            args["t2t_chm13v2"] = os.path.join(temp_dir, "T2T.json.gz")
            write_cdot_json(args["t2t_chm13v2"], "test", [], {}, {}, ["T2T-CHM13v2.0"])

            combine_builds(argparse.Namespace(**args))
            with gzip.open(args["output"]) as f:
                combined = json.load(f)["transcripts"][transcript_accession]["genome_builds"]
        self.assertEqual(combined["GRCh37"]["genome_mismatch"]["transl_except"], {"Arg": [577]})
        self.assertNotIn("genome_mismatch", combined["GRCh38"])

    def test_refseq_genome_mismatch_note(self):
        """ #115 - RefSeq Note counts of how the transcript/protein differ from the genome """
        parse = RefSeq._parse_genome_mismatch_note
        self.assertEqual(parse("The RefSeq transcript has 1 substitution%2C 1 non-frameshifting indel "
                               "compared to this genomic sequence"),
                         {"transcript": {"substitutions": 1, "non_frameshifting_indels": 1}})
        self.assertEqual(parse("The RefSeq transcript has 5 substitutions, 5 non-frameshifting indels "
                               "compared to this genomic sequence"),
                         {"transcript": {"substitutions": 5, "non_frameshifting_indels": 5}})
        self.assertEqual(parse("The RefSeq transcript has 2 substitutions, 1 frameshift and aligns at 99% "
                               "coverage compared to this genomic sequence"),
                         {"transcript": {"substitutions": 2, "frameshifts": 1, "pct_coverage": 99}})
        # Other sentences around it, and the protein
        self.assertEqual(parse("isoform 1 is encoded by transcript variant 1%3B The RefSeq protein has "
                               "1 non-frameshifting indel compared to this genomic sequence"),
                         {"protein": {"non_frameshifting_indels": 1}})
        self.assertEqual(parse("isoform 1 is encoded by transcript variant 1"), {})
        with self.assertLogs(level="WARNING"):
            self.assertEqual(parse("The RefSeq transcript has 1 substitution, something odd compared to "
                                   "this genomic sequence"),
                             {"transcript": {"substitutions": 1}})

    def test_refseq_genome_mismatch_note_and_alignment(self):
        """ #115 - Note counts and the cDNA_match alignment stats go in the build's genome_mismatch """
        parser = GFF3Parser(self.REFSEQ_GFF3_FILENAME_HISTORICAL, "GRCh38", self.FAKE_URL, skip_missing_parents=True)
        _, transcripts = parser.get_genes_and_transcripts()
        genome_mismatch = transcripts["NM_000066.1"]["genome_builds"]["GRCh38"]["genome_mismatch"]
        self.assertEqual(genome_mismatch["transcript"], {"substitutions": 1, "non_frameshifting_indels": 1})
        self.assertEqual(genome_mismatch["alignment"], {"num_mismatch": 1, "gap_count": 1,
                                                        "pct_identity_gap": 99.8998, "pct_coverage": 100})

    def test_refseq_genome_mismatch_alignment_only_when_imperfect(self):
        """ #115 - cDNA_match rows of a perfect alignment don't add genome_mismatch """
        class Feature:
            def __init__(self, **attr):
                self.attr = attr

        builder = TranscriptBuilder("GRCh38", self.FAKE_URL, {})
        perfect = {"num_mismatch": "0", "gap_count": "0", "pct_identity_gap": "100", "pct_coverage": "100"}
        RefSeq()._add_alignment_stats("NM_000059.4", Feature(**perfect), builder)
        self.assertNotIn("NM_000059.4", builder.transcript_genome_mismatch)

        substitution = {**perfect, "num_mismatch": "1", "pct_identity_gap": "99.9"}
        RefSeq()._add_alignment_stats("NM_001754.5", Feature(**substitution), builder)
        self.assertEqual(builder.transcript_genome_mismatch["NM_001754.5"]["alignment"],
                         {"num_mismatch": 1, "gap_count": 0, "pct_identity_gap": 99.9, "pct_coverage": 100})

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
        self.assertEqual(get_transcript_position(True, exons, 3_050), 280)
        # Exons before the hole are unaffected
        self.assertEqual(get_transcript_position(True, exons, 2_050), 150)

        # Same again on the minus strand (exons in stranded order, so genomic order is reversed)
        rev_exons = [
            (3_000, 3_100, 0, 1, 100, None),
            (2_000, 2_100, 1, 101, 200, None),
            (1_000, 1_100, 2, 231, 330, None),
        ]
        self.assertEqual(get_transcript_position(False, rev_exons, 1_050), 280)
        self.assertEqual(get_transcript_position(False, rev_exons, 2_050), 150)

    def test_transcript_position_after_gap_deletion(self):
        """ A Gap 'D' is genomic bases the transcript doesn't have, so they count towards the position
            in the exon. Otherwise the base after a D can't be placed, and gaps after the coordinate are
            applied too, eg the start codon of historical NM_005561.2 (CDS 191..1441) came out as 187 """
        # transcript 1-10 = genomic 1000-1009, genomic 1010-1011 aren't in the transcript,
        # transcript 11-15 = genomic 1012-1016, transcript 16-18 aren't in the genome, 19-28 = genomic 1017-1026
        exons = [(1_000, 1_027, 0, 1, 28, "M10 D2 M5 I3 M10")]
        self.assertEqual(get_transcript_position(True, exons, 1_012), 10)  # 1st base after the D
        self.assertEqual(get_transcript_position(True, exons, 1_015), 13)  # before the I
        self.assertEqual(get_transcript_position(True, exons, 1_017), 18)  # after the I
        # Minus strand: the Gap is in transcript order, the coordinate is the base's end
        self.assertEqual(get_transcript_position(False, exons, 1_015), 10)  # 1st base after the D
        self.assertEqual(get_transcript_position(False, exons, 1_012), 13)  # before the I
        self.assertEqual(get_transcript_position(False, exons, 1_010), 18)  # after the I

    def test_transcript_position_end_before_gap_deletion(self):
        """ A range (eg the stop codon) can end right before a D. Its end is then at the D's first
            base, which is only an error for a base. Eg historical NM_014040.1 (CDS 148..504) had
            stop_codon 503 rather than 504 """
        # Same exon as above: genomic 1010-1011 aren't in the transcript
        exons = [(1_000, 1_027, 0, 1, 28, "M10 D2 M5 I3 M10")]
        self.assertEqual(get_transcript_position(True, exons, 1_010, end=True), 10)
        self.assertEqual(get_transcript_position(False, exons, 1_017, end=True), 10)  # D is genomic 1015-1016
        with self.assertRaises(ValueError):
            get_transcript_position(True, exons, 1_010)

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
