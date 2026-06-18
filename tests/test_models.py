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
# Minimal fake files differing only in gene biotype representation (issue #111):
# the comma-separated str form (<= 0.2.19) and the list form (>= 0.2.20).
GENE_BIOTYPE_STR_JSON = os.path.join(THIS_DIR, "test_data/biotype/cdot.gene_biotype_str.json")
GENE_BIOTYPE_LIST_JSON = os.path.join(THIS_DIR, "test_data/biotype/cdot.gene_biotype_list.json")


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
        # Legacy (<= 0.2.19) str biotype is normalised to list[str] on load - see issue #111
        self.assertEqual(gene.biotype, ["protein_coding"])

    def test_load_ensembl_build_fields(self):
        """Ensembl builds carry a 'tag' and (in this data) no start/stop."""
        data = models.load(ENSEMBL_JSON)
        tx = next(iter(data.transcripts.values()))
        build = next(iter(tx.genome_builds.values()))
        self.assertIn("MANE_Select", build.tag)
        self.assertIsNone(build.start)
        self.assertIsNone(build.stop)
        self.assertEqual(tx.biotype, ["mRNA", "protein_coding"])

    def test_gene_biotype_backwards_compatible(self):
        """Old (str) and new (list) gene biotype both normalise to list[str] - issue #111.

        Gene biotype switched from a comma-separated str (<= 0.2.19) to a list (>= 0.2.20)
        without a detectable schema bump, so both forms appear in the wild.
        """
        expected = ["mRNA", "protein_coding"]

        str_data = models.load(GENE_BIOTYPE_STR_JSON)
        self.assertEqual(str_data.cdot_version, "0.2.19")
        self.assertEqual(str_data.genes["FAKE1"].biotype, expected)

        list_data = models.load(GENE_BIOTYPE_LIST_JSON)
        self.assertEqual(list_data.cdot_version, "0.2.20")
        self.assertEqual(list_data.genes["FAKE1"].biotype, expected)

    def test_gene_from_dict_normalises_biotype(self):
        """The on-demand gene_from_dict path also normalises legacy str biotype."""
        self.assertEqual(models.gene_from_dict({"biotype": "protein_coding"}).biotype, ["protein_coding"])
        self.assertEqual(models.gene_from_dict({"biotype": ["mRNA"]}).biotype, ["mRNA"])
        # Empty / absent biotype stays falsy rather than becoming ['']
        self.assertEqual(models.gene_from_dict({"biotype": ""}).biotype, [])
        self.assertIsNone(models.gene_from_dict({}).biotype)

    def test_genome_build_keeps_current_schema_fields(self):
        """source/ccds/transcript_support_level (data schema >= 0.2.32/0.2.33) must be
        retained, not silently dropped. The typed struct is meant to be a drop-in for
        the plain build dicts, so attribute, item and .get() access must all work."""
        d = {
            "id": "NM_FAKE.1",
            "gene_name": "FAKE1",
            "genome_builds": {
                "GRCh38": {
                    "contig": "NC_000001.11",
                    "strand": "+",
                    "exons": [[100, 200, 1, 0, 100, None]],
                    "source": "BestRefSeq",
                    "ccds": "CCDS123.1",
                    "transcript_support_level": "1",
                }
            },
        }
        build = models.transcript_from_dict(d).genome_builds["GRCh38"]
        self.assertEqual(build.source, "BestRefSeq")
        self.assertEqual(build.ccds, "CCDS123.1")
        self.assertEqual(build.transcript_support_level, "1")
        # dict-compat access (used by the data providers / external subclasses)
        self.assertEqual(build["ccds"], "CCDS123.1")
        self.assertEqual(build.get("source"), "BestRefSeq")

    def test_loads_accepts_genome_build_source_str_or_list(self):
        """Regression: genome_builds[...].source is a str in early 0.2.32 data but a
        list (e.g. ['BestRefSeq']) from 0.2.33 on. The strict msgspec loads() path
        - the one JSONDataProvider uses - must accept both. Before the fix, the list
        form failed to load with msgspec 'Expected `str | null`, got `array`', so
        every current production JSON.gz file was unloadable."""
        def _payload(source):
            return json.dumps({
                "cdot_version": "0.2.33",
                "genome_builds": ["GRCh38"],
                "transcripts": {
                    "NM_FAKE.1": {
                        "id": "NM_FAKE.1",
                        "genome_builds": {
                            "GRCh38": {
                                "contig": "NC_000001.11",
                                "strand": "+",
                                "exons": [[100, 200, 1, 0, 100, None]],
                                "source": source,
                            }
                        },
                    }
                },
            })

        # 0.2.33+ list form - this is exactly what regressed
        list_data = models.loads(_payload(["BestRefSeq"]))
        self.assertEqual(
            list_data.transcripts["NM_FAKE.1"].genome_builds["GRCh38"].source,
            ["BestRefSeq"])

        # early 0.2.32 str form must still load
        str_data = models.loads(_payload("BestRefSeq"))
        self.assertEqual(
            str_data.transcripts["NM_FAKE.1"].genome_builds["GRCh38"].source,
            "BestRefSeq")

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
