""" Check the generated JSON format docs match cdot/models.py (issue #178)

    docs/json_data_format.md and docs/cdot-json-schema.json are produced by
    generate_transcript_data/generate_json_docs.py from the typed models. This renders both into a
    temporary directory and compares them, byte for byte, with the committed files so they cannot
    drift from the models (or from each other) unnoticed.

    If this fails after an intentional change to cdot/models.py or the generator, regenerate with:

        python generate_transcript_data/generate_json_docs.py

    and commit the updated docs/ files.
"""
import difflib
import os
import tempfile
import unittest

from generate_transcript_data import generate_json_docs

REGENERATE_COMMAND = "python generate_transcript_data/generate_json_docs.py"


class GeneratedDocsInSyncTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp_dir = tempfile.TemporaryDirectory()
        cls.generated_paths = generate_json_docs.generate_docs(cls.tmp_dir.name)

    @classmethod
    def tearDownClass(cls):
        cls.tmp_dir.cleanup()

    def test_generates_both_files(self):
        self.assertEqual(sorted(os.path.basename(p) for p in self.generated_paths),
                         sorted([generate_json_docs.SCHEMA_FILENAME, generate_json_docs.MARKDOWN_FILENAME]))

    def test_committed_docs_match_generated(self):
        for generated_path in self.generated_paths:
            filename = os.path.basename(generated_path)
            committed_path = os.path.join(generate_json_docs.DOCS_DIR, filename)
            with self.subTest(filename=filename):
                self.assertTrue(os.path.exists(committed_path),
                                f"{committed_path} is missing. Run: {REGENERATE_COMMAND}")
                with open(generated_path) as f:
                    generated = f.read()
                with open(committed_path) as f:
                    committed = f.read()
                if generated != committed:
                    diff = "".join(difflib.unified_diff(
                        committed.splitlines(keepends=True), generated.splitlines(keepends=True),
                        fromfile=f"docs/{filename} (committed)", tofile=f"docs/{filename} (generated)"))
                    self.fail(f"docs/{filename} is out of date with cdot/models.py / generate_json_docs.py.\n"
                              f"Regenerate it and commit the result:\n\n    {REGENERATE_COMMAND}\n\n{diff}")


if __name__ == "__main__":
    unittest.main()
