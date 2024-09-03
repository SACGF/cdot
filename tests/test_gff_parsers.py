
import os
from inspect import getsourcefile
import unittest
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
    UCSC_GTF_FILENAME = os.path.join(test_data_dir, "hg19_chrY_300kb_genes.gtf")
    FAKE_URL = "http://fake.url"

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
        self.assertEqual(tag, "MANE_Select")

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
