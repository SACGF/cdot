""" cdot_json.py merge_historical / combine_builds, run as the data release Snakefile does

    These stream their input with ijson, which by default reads non-integer numbers as Decimal, which
    json can't write back out. RefSeq alignment stats (eg pct_identity_gap, #115) are floats.
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
GOLDEN_JSON = os.path.join(THIS_FILE_DIR, "test_data", "gff_parser_golden",
                           "refseq_test.historical_RS_2024_08.gff.json")


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

    def test_merge_and_combine_keep_floats(self):
        with open(GOLDEN_JSON) as f:
            golden = json.load(f)
        build_filenames = {}
        for build, arg in [("GRCh37", "--grch37"), ("GRCh38", "--grch38"), ("T2T-CHM13v2.0", "--t2t_chm13v2")]:
            # Golden data is GRCh38, relabel it so there's an input for each build
            source = json.loads(json.dumps(golden))
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
