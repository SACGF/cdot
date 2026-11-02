""" cdot_json.py gtf_to_json / gff3_to_json / merge_historical / combine_builds, run as the data release
    Snakefile does

    merge_historical / combine_builds stream their input with ijson, which by default reads non-integer
    numbers as Decimal, which json can't write back out. RefSeq alignment stats (eg pct_identity_gap, #115)
    are floats.
"""
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from inspect import getsourcefile

THIS_FILE_DIR = os.path.dirname(os.path.abspath(getsourcefile(lambda: 0)))
REPO_DIR = os.path.dirname(THIS_FILE_DIR)
CDOT_JSON = os.path.join(REPO_DIR, "generate_transcript_data", "cdot_json.py")
TEST_DATA_DIR = os.path.join(THIS_FILE_DIR, "test_data")
GOLDEN_DIR = os.path.join(TEST_DATA_DIR, "gff_parser_golden")
GOLDEN_JSON = os.path.join(GOLDEN_DIR, "refseq_test.historical_RS_2024_08.gff.json")


class TestCdotJsonMerge(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def _run(self, *args):
        env = dict(os.environ, PYTHONPATH=REPO_DIR)
        return subprocess.run([sys.executable, CDOT_JSON, *args], env=env, check=True,
                              capture_output=True, text=True).stdout

    def _path(self, filename):
        return os.path.join(self.tmp_dir, filename)

    @staticmethod
    def _load(filename):
        with gzip.open(filename, "rt") as f:
            return json.load(f)

    @classmethod
    def _floats(cls, obj, path=()):
        """ yields (path, value) for every float in nested dicts """
        if isinstance(obj, float):
            yield path, obj
        elif isinstance(obj, dict):
            for key, value in obj.items():
                yield from cls._floats(value, path + (key,))

    def _gene_info(self, gene_info) -> str:
        return self._write("gene_info.json.gz", {"api_retrieval_date": "2024-01-02T03:04:05", "gene_info": gene_info})

    @staticmethod
    def _golden(fixture) -> dict:
        with open(os.path.join(GOLDEN_DIR, fixture + ".json")) as f:
            return json.load(f)

    def test_gff3_to_json(self):
        """ The RefSeq parse (as the goldens check it), with gene info added by GeneID """
        fixture = "refseq_test.transcript_coordinate_hole.gff"
        gene_info_filename = self._gene_info({"99999": {"gene_symbol": "TEST", "summary": "A test gene"}})
        output_filename = self._path("refseq.json.gz")
        url = "http://example.com/refseq.gff"
        stdout = self._run("gff3_to_json", os.path.join(TEST_DATA_DIR, fixture), "--url", url,
                           "--genome-build=GRCh38", "--gene-info-json", gene_info_filename, "--output", output_filename)
        self.assertIn("Using refseq annotation conventions", stdout)
        data = self._load(output_filename)

        golden = self._golden(fixture)
        expected_transcripts = json.loads(json.dumps(golden["transcripts"]).replace("http://fake.url", url))
        self.assertEqual(data["transcripts"], expected_transcripts)
        self.assertEqual(data["genes"]["99999"]["summary"], "A test gene")
        self.assertEqual(data["genes"]["99999"]["gene_symbol"], "TEST")
        self.assertEqual(data["refseq_gene_summary_api_retrieval_date"], "2024-01-02T03:04:05")
        self.assertEqual(data["genome_builds"], ["GRCh38"])
        self.assertEqual(data["metadata"]["method"], "Single GFF3 file")
        self.assertEqual(data["metadata"]["url_counts"], {url: 1})

    def test_gtf_to_json(self):
        """ The Ensembl GRCh37 parse with GENCODE HGNC and the canonical transcripts CSV, and
            --keep-contigs-with-underscores / --no-contig-conversion passed through to the parser """
        fixture = "ensembl_test.GRCh37.87.par.gtf"
        gene_info_filename = self._gene_info({"90001": {"gene_symbol": "PARX", "summary": "A PAR gene"}})
        gencode_filename = self._path("gencode.metadata.HGNC.gz")
        with gzip.open(gencode_filename, "wt") as f:
            f.write("ENST00000900001.1\tPARX\tHGNC:90001\n")
        canonical_filename = self._path("canonical.csv")
        with open(canonical_filename, "w") as f:
            f.write("gene_id,canonical_transcript\nENSG00000900001.1,ENST00000900001.1\n")
        output_filename = self._path("ensembl.json.gz")
        url = "http://example.com/ensembl.gtf"
        stdout = self._run("gtf_to_json", os.path.join(TEST_DATA_DIR, fixture), "--url", url,
                           "--genome-build=GRCh37", "--gene-info-json", gene_info_filename,
                           "--gencode-hgnc-metadata", gencode_filename,
                           "--gene-canonical-transcripts-csv", canonical_filename,
                           "--keep-contigs-with-underscores", "--no-contig-conversion",
                           "--output", output_filename)
        self.assertIn("Using ensembl annotation conventions", stdout)
        data = self._load(output_filename)

        transcript = data["transcripts"]["ENST00000900001.1"]
        build_data = transcript["genome_builds"]["GRCh37"]
        golden_build_data = self._golden(fixture)["transcripts"]["ENST00000900001.1"]["genome_builds"]["GRCh37"]
        self.assertEqual(build_data["exons"], golden_build_data["exons"])
        self.assertEqual(build_data["contig"], "X")  # --no-contig-conversion
        self.assertEqual(build_data["tag"], "basic,Ensembl_canonical")
        self.assertEqual(transcript["hgnc"], "90001")
        gene = data["genes"]["ENSG00000900001.1"]
        self.assertEqual(gene["hgnc"], "90001")
        self.assertEqual(gene["summary"], "A PAR gene")  # By symbol, Ensembl accessions aren't in gene info
        self.assertEqual(data["metadata"]["method"], "Single GTF file")

    def test_merge_and_combine_keep_floats(self):
        with open(GOLDEN_JSON) as f:
            golden = json.load(f)
        build_filenames = {}
        for build, arg in [("GRCh37", "--grch37"), ("GRCh38", "--grch38"), ("T2T-CHM13v2.0", "--t2t_chm13v2")]:
            # Golden data is GRCh38, relabel it so there's an input for each build
            source = json.loads(json.dumps(golden))
            source["genome_builds"] = [build]
            for transcript in source["transcripts"].values():
                transcript["genome_builds"] = {build: transcript["genome_builds"]["GRCh38"]}
            source_filename = self._path(f"source_{build}.json.gz")
            with gzip.open(source_filename, "wt") as f:
                json.dump(source, f)
            build_filenames[arg] = self._path(f"{build}.json.gz")
            self._run("merge_historical", source_filename, f"--genome-build={build}",
                      "--output", build_filenames[arg])

        merged = self._load(build_filenames["--grch38"])
        self.assertEqual(golden["transcripts"].keys(), merged["transcripts"].keys())

        combined_filename = self._path("combined.json.gz")
        build_args = [a for arg, filename in build_filenames.items() for a in (arg, filename)]
        self._run("combine_builds", "--output", combined_filename, *build_args)
        combined = self._load(combined_filename)

        golden_floats = list(self._floats(golden["transcripts"]))
        self.assertTrue(golden_floats, "Test data should contain float values")
        for data in [merged, combined]:
            for path, value in golden_floats:
                result = data["transcripts"]
                for key in path:
                    result = result[key]
                self.assertIsInstance(result, float)
                self.assertEqual(value, result)

        self.assertIn("Method: Combine multiple genome builds",
                      self._run("release_notes", combined_filename, "--show-urls"))

    def _run_error(self, *args) -> str:
        """ Runs cdot_json.py expecting it to fail, returns stderr """
        env = dict(os.environ, PYTHONPATH=REPO_DIR)
        result = subprocess.run([sys.executable, CDOT_JSON, *args], env=env,
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        return result.stderr

    def _write(self, filename, data) -> str:
        path = self._path(filename)
        with gzip.open(path, "wt") as f:
            json.dump(data, f)
        return path

    def test_release_notes_no_metadata(self):
        """ #159 - missing key gave a bare StopIteration rather than the ValueError """
        filename = self._write("no_metadata.json.gz", {"transcripts": {}, "genes": {}})
        stderr = self._run_error("release_notes", filename)
        self.assertIn("ValueError: No metadata in JSON", stderr)

    def test_combine_builds_no_genome_builds(self):
        """ #159 - input without "genome_builds" should give the build mismatch ValueError """
        filename = self._write("no_builds.json.gz", {"transcripts": {}, "genes": {}})
        stderr = self._run_error("combine_builds", "--output", self._path("combined.json.gz"),
                                 "--grch37", filename, "--grch38", filename, "--t2t_chm13v2", filename)
        self.assertIn("JSON file provided for GRCh37 needs to have only GRCh37 data (has None)", stderr)

    def _build_input(self, genome_build):
        """ A single transcript file for one genome build """
        url = "http://example.com/refseq.gff3"
        transcript = self._transcript("NM_007294.3", "BRCA1", "672", url)
        transcript["genome_builds"] = {genome_build: transcript["genome_builds"]["GRCh38"]}
        return self._write(f"{genome_build}.json.gz", {
            "genome_builds": [genome_build],
            "genes": {"672": {"gene_symbol": "BRCA1", "url": url}},
            "transcripts": {"NM_007294.3": transcript},
        })

    def test_combine_builds_subset_of_builds(self):
        """ #174 - any one or more builds can be combined (the Snakefile passes the builds in its YAML) """
        grch38_filename = self._build_input("GRCh38")
        t2t_filename = self._build_input("T2T-CHM13v2.0")

        combined_filename = self._path("combined.json.gz")
        self._run("combine_builds", "--grch38", grch38_filename, "--t2t_chm13v2", t2t_filename,
                  "--output", combined_filename)
        combined = self._load(combined_filename)
        self.assertEqual(["GRCh38", "T2T-CHM13v2.0"], combined["genome_builds"])
        self.assertEqual({"GRCh38", "T2T-CHM13v2.0"}, set(combined["transcripts"]["NM_007294.3"]["genome_builds"]))

        single_filename = self._path("single.json.gz")
        self._run("combine_builds", "--grch38", grch38_filename, "--output", single_filename)
        single = self._load(single_filename)
        self.assertEqual(["GRCh38"], single["genome_builds"])
        self.assertEqual(["GRCh38"], list(single["transcripts"]["NM_007294.3"]["genome_builds"]))

    def test_combine_builds_no_builds(self):
        """ #174 - at least one build argument is needed """
        stderr = self._run_error("combine_builds", "--output", self._path("combined.json.gz"))
        self.assertIn("combine_builds needs at least one genome build JSON file (--grch37, --grch38, --t2t_chm13v2)",
                      stderr)

    @staticmethod
    def _transcript(transcript_accession, gene_symbol, gene_version, url):
        return {
            "id": transcript_accession,
            "gene_name": gene_symbol,
            "gene_version": gene_version,
            "genome_builds": {
                "GRCh38": {"url": url, "contig": "NC_000017.11", "strand": "-",
                           "exons": [[100, 200, 1, 0, 100, None]]},
            },
        }

    def _merge_inputs(self):
        """ An older RefSeq style file with real gene accessions, then a newer UTA style file with
            fake "_" + symbol genes (as written by uta_to_json) """
        refseq_url = "http://example.com/refseq.gff3"
        uta_url = "http://example.com/uta.csv"
        refseq = {
            "genome_builds": ["GRCh38"],
            "genes": {
                "672": {"gene_symbol": "BRCA1", "url": refseq_url},
                # 2 real genes with the same symbol
                "99": {"gene_symbol": "DUPE", "url": refseq_url},
                "100": {"gene_symbol": "DUPE", "url": refseq_url},
            },
            "transcripts": {
                "NM_007294.3": self._transcript("NM_007294.3", "BRCA1", "672", refseq_url),
            },
        }
        uta = {
            "genome_builds": ["GRCh38"],
            "genes": {
                "_BRCA1": {"gene_symbol": "BRCA1", "url": uta_url},
                "_DUPE": {"gene_symbol": "DUPE", "url": uta_url},
                "_NOREAL": {"gene_symbol": "NOREAL", "url": uta_url},
            },
            "transcripts": {
                "NM_007294.3": self._transcript("NM_007294.3", "BRCA1", "_BRCA1", uta_url),
                "NM_007300.4": self._transcript("NM_007300.4", "BRCA1", "_BRCA1", uta_url),
                "NM_000001.1": self._transcript("NM_000001.1", "DUPE", "_DUPE", uta_url),
                "NR_000001.1": self._transcript("NR_000001.1", "NOREAL", "_NOREAL", uta_url),
            },
        }
        return [self._write("refseq.json.gz", refseq), self._write("uta.json.gz", uta)]

    def test_merge_historical_fake_uta_genes(self):
        """ #160 - fake UTA genes are swapped for real ones, deterministically """
        input_filenames = self._merge_inputs()
        outputs = []
        for seed in ["0", "1", "2", "3"]:
            output_filename = self._path(f"merged_{seed}.json.gz")
            env = dict(os.environ, PYTHONPATH=REPO_DIR, PYTHONHASHSEED=seed)
            subprocess.run([sys.executable, CDOT_JSON, "merge_historical", *input_filenames,
                            "--genome-build=GRCh38", "--output", output_filename],
                           env=env, check=True, capture_output=True)
            outputs.append(self._load(output_filename))

        merged = outputs[0]
        gene_versions = {ta: t["gene_version"] for ta, t in merged["transcripts"].items()}
        self.assertEqual({
            "NM_007294.3": "672",  # From the same transcript in the older file
            "NM_007300.4": "672",  # From the gene symbol
            "NM_000001.1": "100",  # Symbol has 2 genes, lowest accession used
            "NR_000001.1": "_NOREAL",  # No real gene, kept fake
        }, gene_versions)
        self.assertEqual({"672", "100", "_NOREAL"}, set(merged["genes"]))
        # Latest transcript data (UTA) is kept
        self.assertEqual("http://example.com/uta.csv",
                         merged["transcripts"]["NM_007294.3"]["genome_builds"]["GRCh38"]["url"])

        for other in outputs[1:]:
            self.assertEqual(merged["transcripts"], other["transcripts"])
            self.assertEqual(merged["genes"], other["genes"])

    def test_merge_historical_missing_gene(self):
        """ #160 - input written with --no-genes gives a clear error """
        filename = self._write("no_genes.json.gz", {
            "genome_builds": ["GRCh38"],
            "transcripts": {
                "NM_007294.3": self._transcript("NM_007294.3", "BRCA1", "672", "http://example.com/refseq.gff3"),
            },
        })
        stderr = self._run_error("merge_historical", filename, "--genome-build=GRCh38",
                                 "--output", self._path("merged.json.gz"))
        self.assertIn("transcript 'NM_007294.3' has gene_version '672' which is not in 'genes'", stderr)
        self.assertIn("--no-genes", stderr)

    def test_merge_historical_wrong_genome_build(self):
        """ #160 - input genome_builds must match --genome-build (as combine_builds checks) """
        input_filenames = self._merge_inputs()
        stderr = self._run_error("merge_historical", *input_filenames, "--genome-build=GRCh37",
                                 "--output", self._path("merged.json.gz"))
        self.assertIn("needs to have only GRCh37 data (has ['GRCh38'])", stderr)

    def test_merge_historical_later_file_overrides(self):
        """ Transcripts in several files take the data from the latest, and only genes still used are kept """
        old_url = "http://example.com/old.gff3"
        new_url = "http://example.com/new.gff3"
        old = {
            "genome_builds": ["GRCh38"],
            "genes": {
                "672": {"gene_symbol": "BRCA1", "url": old_url, "summary": "old summary"},
                "861": {"gene_symbol": "RUNX1", "url": old_url},
            },
            "transcripts": {
                "NM_007294.3": self._transcript("NM_007294.3", "BRCA1", "672", old_url),
                "NM_001754.4": self._transcript("NM_001754.4", "RUNX1", "861", old_url),
            },
        }
        new = {
            "genome_builds": ["GRCh38"],
            "genes": {
                "672": {"gene_symbol": "BRCA1", "url": new_url, "summary": "new summary"},
                "999": {"gene_symbol": "UNUSED", "url": new_url},
            },
            "transcripts": {
                "NM_007294.3": self._transcript("NM_007294.3", "BRCA1", "672", new_url),
                "NM_007294.4": self._transcript("NM_007294.4", "BRCA1", "672", new_url),
            },
        }
        new["transcripts"]["NM_007294.3"]["genome_builds"]["GRCh38"]["exons"] = [[100, 300, 1, 0, 200, None]]
        input_filenames = [self._write("old.json.gz", old), self._write("new.json.gz", new)]
        output_filename = self._path("merged.json.gz")
        self._run("merge_historical", *input_filenames, "--genome-build=GRCh38", "--output", output_filename)
        merged = self._load(output_filename)

        self.assertEqual(sorted(merged["transcripts"]), ["NM_001754.4", "NM_007294.3", "NM_007294.4"])
        brca1_old_version = merged["transcripts"]["NM_007294.3"]["genome_builds"]["GRCh38"]
        self.assertEqual(brca1_old_version["url"], new_url)
        self.assertEqual(brca1_old_version["exons"], [[100, 300, 1, 0, 200, None]])
        self.assertEqual(merged["transcripts"]["NM_001754.4"]["genome_builds"]["GRCh38"]["url"], old_url)
        # Gene data from the latest file, genes without transcripts dropped
        self.assertEqual(sorted(merged["genes"]), ["672", "861"])
        self.assertEqual(merged["genes"]["672"]["summary"], "new summary")
        self.assertEqual(merged["genome_builds"], ["GRCh38"])
        self.assertEqual(merged["metadata"]["method"], "Merge historical")
        self.assertEqual(merged["metadata"]["input_files"], input_filenames)
        self.assertEqual(merged["metadata"]["url_counts"], {new_url: 2, old_url: 1})

    def test_merge_historical_no_genes(self):
        input_filenames = self._merge_inputs()
        output_filename = self._path("merged.json.gz")
        self._run("merge_historical", *input_filenames, "--genome-build=GRCh38", "--no-genes",
                  "--output", output_filename)
        merged = self._load(output_filename)
        self.assertNotIn("genes", merged)
        self.assertEqual(len(merged["transcripts"]), 4)
        # Transcripts still point at (real, where possible) gene accessions
        self.assertEqual(merged["transcripts"]["NM_007300.4"]["gene_version"], "672")

    def test_combine_builds_codon_mismatch(self):
        """ When a transcript's codons differ between builds, the earlier build's coordinates were against
            a different version of the data so are dropped and only the latest build kept """
        grch37_url = "http://example.com/grch37.gff3"
        grch38_url = "http://example.com/grch38.gff3"
        inputs = {}
        for genome_build, arg, url, mismatch_start_codon in [("GRCh37", "--grch37", grch37_url, 10),
                                                              ("GRCh38", "--grch38", grch38_url, 13)]:
            transcripts = {}
            for accession, start_codon in [("NM_007294.3", mismatch_start_codon), ("NM_001754.4", 20)]:
                transcript = self._transcript(accession, "GENE", "1", url)
                transcript["genome_builds"] = {genome_build: transcript["genome_builds"]["GRCh38"]}
                transcript["start_codon"] = start_codon
                transcript["stop_codon"] = 90
                transcripts[accession] = transcript
            inputs[arg] = self._write(f"{genome_build}.json.gz", {
                "genome_builds": [genome_build],
                "genes": {"1": {"gene_symbol": "GENE", "url": url}},
                "transcripts": transcripts,
            })

        combined_filename = self._path("combined.json.gz")
        stdout = self._run("combine_builds", "--grch37", inputs["--grch37"], "--grch38", inputs["--grch38"],
                           "--output", combined_filename)
        combined = self._load(combined_filename)
        mismatch = combined["transcripts"]["NM_007294.3"]
        self.assertEqual(mismatch["start_codon"], 13)
        self.assertEqual(list(mismatch["genome_builds"]), ["GRCh38"])
        same = combined["transcripts"]["NM_001754.4"]
        self.assertEqual(sorted(same["genome_builds"]), ["GRCh37", "GRCh38"])
        self.assertEqual(combined["metadata"]["url_counts"], {grch38_url: 2, grch37_url: 1})
        self.assertIn("Some transcripts were removed as they had different coding coordinates from latest", stdout)
        self.assertIn(f"{grch37_url}: NM_007294.3", stdout)

    def test_release_notes_input_files(self):
        grch38_filename = self._build_input("GRCh38")
        combined_filename = self._path("combined.json.gz")
        self._run("combine_builds", "--grch38", grch38_filename, "--output", combined_filename)
        stdout = self._run("release_notes", combined_filename)
        self.assertEqual(stdout.splitlines(), [
            "### combined.json.gz",
            "Method: Combine multiple genome builds",
            "Input files:",
            f"- {grch38_filename}",
        ])
