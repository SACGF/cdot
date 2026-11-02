""" cdot_json.add_gene_info / add_canonical_transcripts - enrich parsed genes and transcripts from the
    gene info JSON (cdot_gene_info.py) and the Ensembl GRCh37 canonical transcripts CSV """
import gzip
import json
import os
import shutil
import tempfile
import unittest

from generate_transcript_data.cdot_json import add_gene_info, add_canonical_transcripts


class TestAddGeneInfo(unittest.TestCase):
    API_RETRIEVAL_DATE = "2024-01-02T03:04:05"

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.gene_info_filename = os.path.join(self.tmp_dir, "gene_info.json.gz")
        with gzip.open(self.gene_info_filename, "wt") as f:
            json.dump({
                "api_retrieval_date": self.API_RETRIEVAL_DATE,
                "gene_info": {
                    "672": {"gene_symbol": "BRCA1", "map_location": "17q21.31", "summary": "BRCA1 summary",
                            "description": "BRCA1 DNA repair associated"},
                    "861": {"gene_symbol": "RUNX1", "map_location": "21q22.12", "summary": "RUNX1 summary"},
                },
            }, f)

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def test_refseq_by_gene_id(self):
        genes = {"672": {"gene_symbol": "BRCA1", "description": "from the GFF"}, "99999": {"gene_symbol": "UNKNOWN"}}
        date = add_gene_info(self.gene_info_filename, genes)
        self.assertEqual(date, self.API_RETRIEVAL_DATE)
        self.assertEqual(genes["672"], {"gene_symbol": "BRCA1", "map_location": "17q21.31", "summary": "BRCA1 summary",
                                        "description": "BRCA1 DNA repair associated"})  # gene info wins
        self.assertEqual(genes["99999"], {"gene_symbol": "UNKNOWN"})

    def test_ensembl_by_symbol(self):
        """ Ensembl gene accessions aren't in the (Entrez ID keyed) gene info, so look up the symbol """
        genes = {"ENSG00000159216.22": {"gene_symbol": "RUNX1"}, "ENSG00000000001.1": {"gene_symbol": None}}
        add_gene_info(self.gene_info_filename, genes)
        self.assertEqual(genes["ENSG00000159216.22"],
                         {"gene_symbol": "RUNX1", "map_location": "21q22.12", "summary": "RUNX1 summary"})
        self.assertEqual(genes["ENSG00000000001.1"], {"gene_symbol": None})

    def test_no_genes(self):
        """ Nothing to add, so the file isn't even read """
        self.assertIsNone(add_gene_info(os.path.join(self.tmp_dir, "missing.json.gz"), {}))


class TestAddCanonicalTranscripts(unittest.TestCase):
    GENOME_BUILD = "GRCh37"

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.csv_filename = os.path.join(self.tmp_dir, "canonical.csv")

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def _write(self, content):
        with open(self.csv_filename, "w") as f:
            f.write(content)

    def _transcript(self, gene_version, tag=None):
        build_data = {"url": "http://fake.url", "contig": "NC_000021.8", "strand": "-", "exons": []}
        if tag:
            build_data["tag"] = tag
        return {"gene_version": gene_version, "genome_builds": {self.GENOME_BUILD: build_data}}

    def test_canonical_tag(self):
        self._write("gene_id,canonical_transcript\n"
                    "ENSG00000159216.14,ENST00000300305.3\n"
                    "ENSG00000012048.15,ENST00000357654.3\n")
        transcripts = {
            "ENST00000300305.3": self._transcript("ENSG00000159216.14", tag="basic,CCDS"),
            "ENST00000344691.4": self._transcript("ENSG00000159216.14", tag="basic"),
            "ENST00000357654.3": self._transcript("ENSG00000012048.15"),
            "ENST00000000001.1": self._transcript("ENSG00000000001.1"),  # Gene not in the CSV
        }
        add_canonical_transcripts(self.csv_filename, self.GENOME_BUILD, transcripts)
        tags = {accession: t["genome_builds"][self.GENOME_BUILD].get("tag") for accession, t in transcripts.items()}
        self.assertEqual(tags, {
            "ENST00000300305.3": "basic,CCDS,Ensembl_canonical",
            "ENST00000344691.4": "basic",
            "ENST00000357654.3": "Ensembl_canonical",
            "ENST00000000001.1": None,
        })

    def test_bad_header(self):
        self._write("gene,transcript\nENSG00000159216.14,ENST00000300305.3\n")
        with self.assertRaisesRegex(ValueError, "unexpected header line: 'gene,transcript'"):
            add_canonical_transcripts(self.csv_filename, self.GENOME_BUILD, {"ENST00000300305.3": {}})

    def test_no_csv(self):
        """ Only Ensembl GRCh37 passes a CSV """
        transcripts = {"NM_000059.4": self._transcript("675")}
        add_canonical_transcripts(None, self.GENOME_BUILD, transcripts)
        self.assertNotIn("tag", transcripts["NM_000059.4"]["genome_builds"][self.GENOME_BUILD])
