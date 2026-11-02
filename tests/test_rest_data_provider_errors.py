"""
Tests for RESTDataProvider HTTP error handling, timeouts and supported assemblies.

Issues SACGF/cdot#154 (transient errors cached as missing, no timeouts, missing Content-Type)
and SACGF/cdot#155 (T2T-CHM13v2.0 contigs rejected).

No network: the provider's requests.Session is replaced by a mock.
"""
from unittest import mock

import pytest
import requests

from cdot.hgvs.dataproviders.json_data_provider import RESTDataProvider


def _response(status_code, json_data=None, content_type="application/json"):
    response = requests.Response()
    response.status_code = status_code
    response.url = "https://example.org/mock"
    if content_type is not None:
        response.headers["Content-Type"] = content_type
    if json_data is not None:
        response.json = lambda: json_data
    return response


@pytest.fixture
def dp():
    provider = RESTDataProvider(url="https://example.org", session=mock.Mock(spec=requests.Session))
    return provider


def test_transient_error_raises_and_is_not_cached(dp):
    dp.session.get.return_value = _response(503)
    with pytest.raises(requests.HTTPError):
        dp._get_transcript("NM_000059.4")
    assert "NM_000059.4" not in dp.transcripts

    # Server recovers - the next call goes back to the server
    dp.session.get.return_value = _response(200, {"id": "NM_000059.4"})
    assert dp._get_transcript("NM_000059.4") == {"id": "NM_000059.4"}
    assert dp.session.get.call_count == 2


def test_rate_limited_gene_raises_and_is_not_cached(dp):
    dp.session.get.return_value = _response(429)
    with pytest.raises(requests.HTTPError):
        dp._get_gene("BRCA2")
    assert "BRCA2" not in dp.genes


def test_404_is_cached_as_missing(dp):
    dp.session.get.return_value = _response(404)
    assert dp._get_transcript("NM_999999.9") is None
    assert dp.transcripts["NM_999999.9"] is None
    assert dp._get_transcript("NM_999999.9") is None
    assert dp.session.get.call_count == 1  # second lookup is a cache hit


def test_missing_content_type_raises_value_error(dp):
    dp.session.get.return_value = _response(200, content_type=None)
    with pytest.raises(ValueError, match="Non-json response"):
        dp._get_transcript("NM_000059.4")


def test_post_missing_content_type_raises_value_error(dp):
    dp.session.post.return_value = _response(200, content_type=None)
    with pytest.raises(ValueError, match="Non-json response"):
        dp._post_to_url("https://example.org/transcripts", {"ids": ["NM_000059.4"]})


def test_default_timeout_passed_through(dp):
    dp.session.get.return_value = _response(200, {"id": "NM_000059.4"})
    dp._get_transcript("NM_000059.4")
    assert dp.session.get.call_args.kwargs["timeout"] == RESTDataProvider.DEFAULT_TIMEOUT

    dp.session.post.return_value = _response(200, {"NM_007294.4": {"id": "NM_007294.4"}})
    dp.prefetch(["NM_007294.4"])
    assert dp.session.post.call_args.kwargs["timeout"] == RESTDataProvider.DEFAULT_TIMEOUT


def test_custom_timeout_passed_through():
    session = mock.Mock(spec=requests.Session)
    session.get.return_value = _response(200, {"id": "NM_000059.4"})
    dp = RESTDataProvider(url="https://example.org", timeout=5, session=session)
    dp._get_transcript("NM_000059.4")
    assert session.get.call_args.kwargs["timeout"] == 5


def test_uses_a_session_by_default():
    dp = RESTDataProvider(url="https://example.org")
    assert isinstance(dp.session, requests.Session)


def test_prefetch_413_still_shrinks_and_retries(dp):
    def fake_post(url, json=None, timeout=None):
        if len(json["ids"]) > 2:
            return _response(413)
        return _response(200, {ac: {"id": ac} for ac in json["ids"]})

    dp.session.post.side_effect = fake_post
    accs = ["NM_000059.4", "NM_001754.5", "NM_007294.4", "NM_000546.6"]
    assert dp.prefetch(accs, batch_size=4) == 4
    assert all(dp.transcripts[ac] == {"id": ac} for ac in accs)


def test_prefetch_404_endpoint_falls_back_to_singles(dp):
    dp.session.post.return_value = _response(404)

    def fake_get(url, timeout=None):
        ac = url.rsplit("/", 1)[1]
        if ac == "NM_999999.9":
            return _response(404)
        if ac == "NM_001754.5":
            return _response(503)
        return _response(200, {"id": ac})

    dp.session.get.side_effect = fake_get
    n = dp.prefetch(["NM_000059.4", "NM_999999.9", "NM_001754.5"])
    assert n == 2
    assert dp.transcripts["NM_000059.4"] == {"id": "NM_000059.4"}
    assert dp.transcripts["NM_999999.9"] is None  # 404 cached as missing
    assert "NM_001754.5" not in dp.transcripts     # 503 not cached


def test_t2t_contig_supported():
    dp = RESTDataProvider(url="https://example.org")
    assert "T2T-CHM13v2.0" in dp.assembly_maps
    assert dp.assembly_by_contig["NC_060925.1"] == "T2T-CHM13v2.0"  # T2T chr1
    assert dp.assembly_by_contig["NC_000001.11"] == "GRCh38"
    assert dp.assembly_by_contig["NC_000001.10"] == "GRCh37"


def test_t2t_transcript_coordinates():
    dp = RESTDataProvider(url="https://example.org")
    build_coords = {"contig": "NC_060925.1", "strand": "+", "exons": [[100, 200, 0, 1, 101, None]]}
    dp.transcripts["NM_000059.4"] = {"id": "NM_000059.4", "genome_builds": {"T2T-CHM13v2.0": build_coords}}
    assert dp._get_transcript_coordinates_for_contig(dp.transcripts["NM_000059.4"], "NC_060925.1") == build_coords


def test_assemblies_argument():
    dp = RESTDataProvider(url="https://example.org", assemblies=["GRCh38"])
    assert list(dp.assembly_maps) == ["GRCh38"]


def test_timeout_propagates_and_is_not_cached(dp):
    dp.session.get.side_effect = requests.Timeout("read timed out")
    with pytest.raises(requests.Timeout):
        dp._get_transcript("NM_000059.4")
    assert "NM_000059.4" not in dp.transcripts

    # The next call goes back to the server
    dp.session.get.side_effect = None
    dp.session.get.return_value = _response(200, {"id": "NM_000059.4"})
    assert dp._get_transcript("NM_000059.4") == {"id": "NM_000059.4"}


def test_connection_error_on_gene_propagates_and_is_not_cached(dp):
    dp.session.get.side_effect = requests.ConnectionError("refused")
    with pytest.raises(requests.ConnectionError):
        dp._get_gene("BRCA2")
    assert "BRCA2" not in dp.genes


def test_server_error_on_get_tx_versions_raises(dp):
    dp.session.get.return_value = _response(500)
    with pytest.raises(requests.HTTPError):
        dp.get_tx_versions("NM_000059")
    assert dp.transcripts == {}


@pytest.mark.parametrize("content_type", ["text/html; charset=utf-8", "text/plain", ""])
def test_non_json_content_type_raises_value_error(dp, content_type):
    dp.session.get.return_value = _response(200, {"id": "NM_000059.4"}, content_type=content_type)
    with pytest.raises(ValueError, match="Non-json response"):
        dp._get_transcript("NM_000059.4")
    assert "NM_000059.4" not in dp.transcripts


def test_post_non_json_content_type_raises_value_error(dp):
    dp.session.post.return_value = _response(200, {}, content_type="text/html")
    with pytest.raises(ValueError, match="Non-json response"):
        dp._post_to_url("https://example.org/transcripts", {"ids": ["NM_000059.4"]})


def test_json_content_type_with_charset_accepted(dp):
    dp.session.get.return_value = _response(200, {"id": "NM_000059.4"},
                                            content_type="application/json; charset=utf-8")
    assert dp._get_transcript("NM_000059.4") == {"id": "NM_000059.4"}


def test_post_405_means_endpoint_absent(dp):
    dp.session.post.return_value = _response(405)
    assert dp._post_to_url("https://example.org/transcripts", {"ids": ["NM_000059.4"]}) is None


def test_prefetch_server_error_propagates(dp):
    dp.session.post.return_value = _response(500)
    with pytest.raises(requests.HTTPError):
        dp.prefetch(["NM_000059.4", "NM_001754.5"])
    assert dp.transcripts == {}


def test_prefetch_timeout_propagates(dp):
    dp.session.post.side_effect = requests.Timeout("read timed out")
    with pytest.raises(requests.Timeout):
        dp.prefetch(["NM_000059.4"])
    assert dp.transcripts == {}
