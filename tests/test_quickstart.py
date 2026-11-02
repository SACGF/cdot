"""
Runs examples/quickstart.py's code path against the test data, so the quickstart (and the README
example it mirrors) is pinned by the suite without downloading a data release.
"""
import importlib.util
import os
from unittest.mock import patch

import pytest
from hgvs.dataproviders.seqfetcher import SeqFetcher

from cdot.hgvs.dataproviders import ChainedSeqFetcher, JSONDataProvider
from tests.mock_seqfetcher import MockSeqFetcher

TESTS_DIR = os.path.dirname(__file__)
TEST_DATA_DIR = os.path.join(TESTS_DIR, "test_data")
QUICKSTART_PY = os.path.join(os.path.dirname(TESTS_DIR), "examples", "quickstart.py")


def _load_quickstart():
    spec = importlib.util.spec_from_file_location("quickstart", QUICKSTART_PY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def quickstart():
    return _load_quickstart()


@pytest.fixture(scope="module")
def hdp():
    # Same setup as test_json_data_provider_refseq: transcript sequences from the mock, genome from SeqFetcher
    mock_seqfetcher = MockSeqFetcher(os.path.join(TEST_DATA_DIR, "transcript_sequences.json"))
    seqfetcher = ChainedSeqFetcher(mock_seqfetcher, SeqFetcher())
    return JSONDataProvider([os.path.join(TEST_DATA_DIR, "cdot.refseq.grch37.json")], seqfetcher=seqfetcher)


def test_resolve_fixes_then_maps_to_genome(quickstart, hdp):
    examples = [
        "AOAH(NM_001637.3):c.1582G>A p.(Arg528His)",  # gene/transcript swapped, trailing p.
        "nm_001637.3 c.1582G>A",  # lowercase prefix, missing colon
    ]
    results = list(quickstart.resolve(hdp, "GRCh37", examples))
    assert [(original, fixed, str(var_g)) for original, fixed, _fixes, var_g in results] == [
        ("AOAH(NM_001637.3):c.1582G>A p.(Arg528His)", "NM_001637.3(AOAH):c.1582G>A", "NC_000007.13:g.36561662C>T"),
        ("nm_001637.3 c.1582G>A", "NM_001637.3:c.1582G>A", "NC_000007.13:g.36561662C>T"),
    ]
    for _original, _fixed, fixes, _var_g in results:
        assert fixes  # every example needed at least one fix


def test_build_provider_default_seqfetcher(quickstart):
    hdp = quickstart.build_provider(os.path.join(TEST_DATA_DIR, "cdot.refseq.grch37.json"))
    assert hdp.get_tx_identity_info("NM_001637.3")["hgnc"] == "AOAH"


def test_download_latest_data_skips_existing_file(quickstart, tmp_path):
    existing = tmp_path / "cdot-0.2.36.refseq.GRCh38.json.gz"
    existing.write_bytes(b"")
    url = f"https://github.com/SACGF/cdot/releases/download/data_v0.2.36/{existing.name}"
    with patch.object(quickstart, "get_latest_combo_file_urls", return_value=[url]), \
            patch.object(quickstart.requests, "get", side_effect=AssertionError("must not download")):
        assert quickstart.download_latest_data(str(tmp_path)) == str(existing)
