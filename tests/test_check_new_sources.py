""" check_new_sources.py, offline: fetch is patched to fail so nothing talks to the FTP sites """
import io
import re
import sys
import unittest
from contextlib import redirect_stdout
from unittest.mock import patch

import pytest

pytest.importorskip("yaml")  # check_new_sources reads cdot_transcripts.yaml, pyyaml is not a test dependency

from generate_transcript_data import check_new_sources  # noqa: E402
from generate_transcript_data.check_new_sources import NoSourcesConfigured, latest_configured  # noqa: E402

BUILDS = ["GRCh37", "GRCh38", "T2T-CHM13v2.0"]
EMPTY_CONFIG = {consortium: {build: {} for build in BUILDS} for consortium in ["ensembl", "refseq"]}


class TestCheckNewSources(unittest.TestCase):
    def setUp(self):
        patcher = patch.object(check_new_sources, "fetch", side_effect=AssertionError("no network in tests"))
        self.mock_fetch = patcher.start()
        self.addCleanup(patcher.stop)

    def test_latest_configured(self):
        self.assertEqual(latest_configured([81, 116, 100]), 116)
        self.assertEqual(latest_configured(["RS_2024_08", "RS_2025_08"]), "RS_2025_08")
        with self.assertRaises(NoSourcesConfigured):
            latest_configured([])

    def test_checks_without_sources(self):
        """ Each check reports the configuration problem before looking upstream """
        checks = [check_new_sources.check_ensembl_grch38, check_new_sources.check_ensembl_grch37,
                  check_new_sources.check_ensembl_t2t, check_new_sources.check_refseq_historical,
                  check_new_sources.check_uta]
        for check in checks:
            with self.subTest(check=check.__name__), self.assertRaises(NoSourcesConfigured):
                check(EMPTY_CONFIG)
        self.mock_fetch.assert_not_called()

    def test_main_reports_no_sources(self):
        with patch.object(check_new_sources, "load_config", return_value=EMPTY_CONFIG), \
                patch.object(sys, "argv", ["check_new_sources.py"]), \
                redirect_stdout(io.StringIO()) as stdout:
            with self.assertRaises(SystemExit) as cm:
                check_new_sources.main()
        self.assertEqual(cm.exception.code, 1)
        output = stdout.getvalue()
        for check in ["check_ensembl_grch38", "check_ensembl_grch37", "check_ensembl_t2t", "check_refseq_historical",
                      "check_uta"]:
            self.assertRegex(output, rf"(?m)^{check}\s+no sources configured in cdot_transcripts\.yaml$")
        # A failed upstream listing (here, no network) is a different problem and still reported as such
        self.assertRegex(output, r"(?m)^check_refseq\s+ERROR .*\(check the listing by hand\)$")

    def test_yaml_configures_every_source(self):
        """ The checked-in cdot_transcripts.yaml has the layout the checks expect """
        config = check_new_sources.load_config()
        self.assertEqual(sorted(config), ["ensembl", "refseq"])
        for consortium in config:
            for build in BUILDS:
                with self.subTest(consortium=consortium, build=build):
                    self.assertTrue(check_new_sources.source_urls(config, consortium, build))
        for check in [check_new_sources.check_ensembl_grch38, check_new_sources.check_ensembl_grch37,
                      check_new_sources.check_ensembl_t2t, check_new_sources.check_refseq_historical,
                      check_new_sources.check_uta]:
            with self.subTest(check=check.__name__), self.assertRaises(AssertionError) as cm:
                check(config)  # Gets as far as the (patched) network
            self.assertEqual(str(cm.exception), "no network in tests")

    def test_report(self):
        with redirect_stdout(io.StringIO()) as stdout:
            self.assertFalse(check_new_sources.report("Ensembl GRCh38", "116", "116", []))
            self.assertTrue(check_new_sources.report("Ensembl GRCh38", "116", "118", ["117", "118"], ["a note"]))
        lines = stdout.getvalue().splitlines()
        self.assertEqual(re.sub(r"\s+", " ", lines[0]), "Ensembl GRCh38 yaml has 116 upstream latest 116 up to date")
        self.assertEqual(re.sub(r"\s+", " ", lines[1]), "Ensembl GRCh38 yaml has 116 upstream latest 118 NEW: 117, 118")
        self.assertEqual(lines[2].strip(), "! a note")
