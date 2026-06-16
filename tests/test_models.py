import json
import os
import unittest
from inspect import getsourcefile
from os.path import abspath

# Optional typed layer (issue #37) - skip the whole module if msgspec isn't installed
try:
    from cdot import models
    HAVE_MSGSPEC = True
except ImportError:
    HAVE_MSGSPEC = False

THIS_DIR = os.path.dirname(abspath(getsourcefile(lambda: 0)))
REFSEQ_JSON = os.path.join(THIS_DIR, "test_data/cdot.refseq.grch37.json")
ENSEMBL_JSON = os.path.join(THIS_DIR, "test_data/cdot.ensembl.grch38.json")


@unittest.skipUnless(HAVE_MSGSPEC, "msgspec not installed")
class TestModels(unittest.TestCase):
    def test_load_refseq(self):
        data = models.load(REFSEQ_JSON)
        self.assertEqual(data.cdot_version, "0.2.10")
        self.assertEqual(data.genome_builds, ["GRCh37"])

        tx = data.transcripts["NM_001637.3"]
        self.assertEqual(tx.id, "NM_001637.3")
        self.assertEqual(tx.gene_name, "AOAH")
        self.assertEqual(tx.protein, "NP_001628.1")
        self.assertEqual(tx.start_codon, 401)
        self.assertEqual(tx.stop_codon, 2129)
        self.assertEqual(tx.biotype, ["protein_coding"])

        build = tx.genome_builds["GRCh37"]
        self.assertEqual(build.contig, "NC_000007.13")
        self.assertEqual(build.strand, "-")
        self.assertEqual(build.cds_start, 36552857)
        self.assertEqual(build.cds_end, 36763753)

    def test_exon_array_is_named(self):
        """The positional exon array is exposed with real field names."""
        data = models.load(REFSEQ_JSON)
        build = data.transcripts["NM_001637.3"].genome_builds["GRCh37"]

        first = build.exons[0]
        self.assertEqual(first.alt_start, 36552548)
        self.assertEqual(first.alt_end, 36552986)
        self.assertEqual(first.exon_id, 20)
        self.assertEqual(first.cds_start, 2001)
        self.assertEqual(first.cds_end, 2440)
        self.assertEqual(first.gap, "M196 I1 M61 I1 M181")

        # An exon that aligns cleanly has gap == None
        self.assertIsNone(build.exons[1].gap)

    def test_non_coding_transcript_optionals(self):
        """A non-coding transcript has no protein / codon fields."""
        data = models.load(REFSEQ_JSON)
        # The second refseq transcript in the test data is non-coding
        non_coding = [t for t in data.transcripts.values() if t.protein is None]
        self.assertTrue(non_coding, "expected a non-coding transcript in test data")
        tx = non_coding[0]
        self.assertIsNone(tx.protein)
        self.assertIsNone(tx.start_codon)
        self.assertIsNone(tx.stop_codon)

    def test_gene_info(self):
        data = models.load(REFSEQ_JSON)
        gene = data.genes["GATA2"]
        self.assertEqual(gene.gene_symbol, "GATA2")
        self.assertEqual(gene.map_location, "3q21.3")
        self.assertEqual(gene.hgnc, "4171")
        self.assertIn("MONOMAC", gene.aliases)
        self.assertEqual(gene.biotype, "protein_coding")

    def test_load_ensembl_build_fields(self):
        """Ensembl builds carry a 'tag' and (in this data) no start/stop."""
        data = models.load(ENSEMBL_JSON)
        tx = next(iter(data.transcripts.values()))
        build = next(iter(tx.genome_builds.values()))
        self.assertIn("MANE_Select", build.tag)
        self.assertIsNone(build.start)
        self.assertIsNone(build.stop)
        self.assertEqual(tx.biotype, ["mRNA", "protein_coding"])

    def test_transcript_from_dict(self):
        """On-demand path: build a typed Transcript from a plain dict."""
        with open(REFSEQ_JSON) as f:
            raw = json.load(f)
        d = raw["transcripts"]["NM_001637.3"]
        tx = models.transcript_from_dict(d)
        self.assertEqual(tx.gene_name, "AOAH")
        self.assertEqual(tx.genome_builds["GRCh37"].exons[0].exon_id, 20)

    def test_typed_matches_raw_for_all_transcripts(self):
        """Data-driven: the typed view faithfully represents the raw JSON."""
        for path in (REFSEQ_JSON, ENSEMBL_JSON):
            with open(path) as f:
                raw = json.load(f)
            data = models.load(path)
            self.assertEqual(set(data.transcripts), set(raw["transcripts"]))
            for tx_id, raw_tx in raw["transcripts"].items():
                tx = data.transcripts[tx_id]
                self.assertEqual(tx.gene_name, raw_tx.get("gene_name"))
                self.assertEqual(tx.protein, raw_tx.get("protein"))
                self.assertEqual(tx.start_codon, raw_tx.get("start_codon"))
                self.assertEqual(set(tx.genome_builds), set(raw_tx["genome_builds"]))
                for build_name, raw_build in raw_tx["genome_builds"].items():
                    build = tx.genome_builds[build_name]
                    self.assertEqual(build.contig, raw_build["contig"])
                    self.assertEqual(build.strand, raw_build["strand"])
                    self.assertEqual(len(build.exons), len(raw_build["exons"]))
                    for exon, raw_exon in zip(build.exons, raw_build["exons"]):
                        self.assertEqual(
                            [exon.alt_start, exon.alt_end, exon.exon_id,
                             exon.cds_start, exon.cds_end, exon.gap],
                            raw_exon,
                        )


if __name__ == "__main__":
    unittest.main()
