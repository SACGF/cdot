"""
Loading cdot JSON files closes the file handles it opens, and leaves caller-supplied file objects open (#169)
"""
import gc
import gzip
import io
import os
import shutil
import warnings
from inspect import getsourcefile
from os.path import abspath

import pytest

from cdot import models
from cdot.hgvs.dataproviders import JSONDataProvider
from cdot.pyhgvs.pyhgvs_transcript import JSONPyHGVSTranscriptFactory

THIS_DIR = os.path.dirname(abspath(getsourcefile(lambda: 0)))
REFSEQ_JSON = os.path.join(THIS_DIR, "test_data/cdot.refseq.grch37.json")
TRANSCRIPT_ID = "NM_001637.3"


@pytest.fixture(scope="module")
def refseq_json_gz(tmp_path_factory):
    gz_filename = str(tmp_path_factory.mktemp("json_loading") / "cdot.refseq.grch37.json.gz")
    with open(REFSEQ_JSON, "rb") as f_in, gzip.open(gz_filename, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    return gz_filename


def _assert_no_resource_warning(construct):
    # An unclosed handle emits its ResourceWarning from __del__ when it is collected, which can't
    # propagate as an exception, so record the warnings and check none were emitted
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always", ResourceWarning)
        obj = construct()
        del obj
        gc.collect()
    resource_warnings = [str(w.message) for w in record if issubclass(w.category, ResourceWarning)]
    assert not resource_warnings


@pytest.mark.parametrize("filename", [REFSEQ_JSON, "gz"])
def test_json_data_provider_closes_files(filename, refseq_json_gz):
    if filename == "gz":
        filename = refseq_json_gz
    _assert_no_resource_warning(lambda: JSONDataProvider([filename]))


@pytest.mark.parametrize("filename", [REFSEQ_JSON, "gz"])
def test_pyhgvs_factory_closes_files(filename, refseq_json_gz):
    if filename == "gz":
        filename = refseq_json_gz
    _assert_no_resource_warning(lambda: JSONPyHGVSTranscriptFactory([filename]))


def test_json_data_provider_gz_matches_plain(refseq_json_gz):
    plain = JSONDataProvider([REFSEQ_JSON])
    gz = JSONDataProvider([refseq_json_gz])
    assert gz.get_tx_identity_info(TRANSCRIPT_ID) == plain.get_tx_identity_info(TRANSCRIPT_ID)


def test_pyhgvs_factory_gz_matches_plain(refseq_json_gz):
    plain = JSONPyHGVSTranscriptFactory([REFSEQ_JSON])
    gz = JSONPyHGVSTranscriptFactory([refseq_json_gz])
    assert gz._get_transcript(TRANSCRIPT_ID) == plain._get_transcript(TRANSCRIPT_ID)


def test_file_objects_are_left_open():
    """ A caller-supplied file object belongs to the caller, so loading must not close it """
    with open(REFSEQ_JSON, "rb") as f:
        hdp = JSONDataProvider([f])
        assert not f.closed
        f.seek(0)
        factory = JSONPyHGVSTranscriptFactory([f])
        assert not f.closed
    assert TRANSCRIPT_ID in hdp.transcripts
    assert factory._get_transcript(TRANSCRIPT_ID)["id"] == TRANSCRIPT_ID


def test_open_json_closes_filename_but_not_file_object():
    with models.open_json(REFSEQ_JSON) as f:
        assert f.read(1) == b"{"
    assert f.closed

    bio = io.BytesIO(b"{}")
    with models.open_json(bio) as f:
        assert f is bio
    assert not bio.closed
