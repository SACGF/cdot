"""
Tests for RESTDataProvider read-ahead cache warming (prefetch / prefetch_from_hgvs).

Issues SACGF/cdot#108 (batch endpoint) and SACGF/cdot#109 (read-ahead pre-warm).

The REST layer is mocked - no network. We monkeypatch the provider's _post_to_url /
_get_from_url seams (the same seams the existing get_tx_ac_tags_for_gene tests use) so we can
assert on what was requested and what landed in the cache.
"""
import pytest

from cdot.hgvs.dataproviders.json_data_provider import RESTDataProvider


@pytest.fixture
def dp():
    return RESTDataProvider(url="https://example.org")


# ---------------------------------------------------------------------------
# prefetch() - batch endpoint path
# ---------------------------------------------------------------------------

def test_prefetch_batch_populates_cache(dp, monkeypatch):
    captured = {}

    def fake_post(url, json_body):
        captured["url"] = url
        captured["ids"] = json_body["ids"]
        return {"NM_000059.3": {"id": "NM_000059.3"},
                "NM_007294.3": {"id": "NM_007294.3"}}

    monkeypatch.setattr(dp, "_post_to_url", fake_post)

    n = dp.prefetch(["NM_000059.3", "NM_007294.3"])

    assert captured["url"] == "https://example.org/transcripts"
    assert sorted(captured["ids"]) == ["NM_000059.3", "NM_007294.3"]
    assert n == 2
    # Every accession is now a cache hit
    assert dp._get_transcript("NM_000059.3") == {"id": "NM_000059.3"}
    assert dp._get_transcript("NM_007294.3") == {"id": "NM_007294.3"}


def test_prefetch_batch_missing_stored_as_none(dp, monkeypatch):
    # Server returns null for versioned misses (matches the "store None" cache behaviour)
    monkeypatch.setattr(dp, "_post_to_url",
                        lambda url, json_body: {"NM_000059.3": {"id": "NM_000059.3"},
                                                "NM_999999.9": None})
    n = dp.prefetch(["NM_000059.3", "NM_999999.9"])
    assert n == 2
    assert dp.transcripts["NM_999999.9"] is None
    # A cached None must NOT trigger another fetch
    monkeypatch.setattr(dp, "_get_from_url",
                        lambda url: pytest.fail("should not refetch a cached miss"))
    assert dp._get_transcript("NM_999999.9") is None


def test_prefetch_batch_versionless_expands(dp, monkeypatch):
    # Versionless id expands server-side to all versions, keyed by full accession
    monkeypatch.setattr(dp, "_post_to_url",
                        lambda url, json_body: {"NM_000059.3": {"id": "NM_000059.3"},
                                                "NM_000059.4": {"id": "NM_000059.4"}})
    n = dp.prefetch(["NM_000059"])
    assert n == 2  # both versions cached from one request
    assert dp._get_transcript("NM_000059.4") == {"id": "NM_000059.4"}


def test_prefetch_skips_already_cached(dp, monkeypatch):
    dp.transcripts["NM_000059.3"] = {"id": "NM_000059.3"}  # pre-cached

    captured = {}

    def fake_post(url, json_body):
        captured["ids"] = json_body["ids"]
        return {"NM_007294.3": {"id": "NM_007294.3"}}

    monkeypatch.setattr(dp, "_post_to_url", fake_post)

    n = dp.prefetch(["NM_000059.3", "NM_007294.3"])
    assert captured["ids"] == ["NM_007294.3"]  # cached accession not re-requested
    assert n == 1


def test_prefetch_empty_input_no_request(dp, monkeypatch):
    monkeypatch.setattr(dp, "_post_to_url",
                        lambda url, json_body: pytest.fail("no request for empty/all-cached"))
    assert dp.prefetch([]) == 0


def test_prefetch_all_cached_no_request(dp, monkeypatch):
    dp.transcripts["NM_000059.3"] = {"id": "NM_000059.3"}
    monkeypatch.setattr(dp, "_post_to_url",
                        lambda url, json_body: pytest.fail("no request when all cached"))
    assert dp.prefetch(["NM_000059.3"]) == 0


# ---------------------------------------------------------------------------
# prefetch() - fallback to concurrent singles when no batch endpoint
# ---------------------------------------------------------------------------

def test_prefetch_falls_back_when_no_batch_endpoint(dp, monkeypatch):
    # _post_to_url returns None => endpoint absent (404/405) => concurrent fallback
    monkeypatch.setattr(dp, "_post_to_url", lambda url, json_body: None)

    requested = []

    def fake_get(url):
        requested.append(url)
        ac = url.rsplit("/", 1)[1]
        return {"id": ac}

    monkeypatch.setattr(dp, "_get_from_url", fake_get)

    n = dp.prefetch(["NM_000059.3", "NM_007294.3"])
    assert n == 2
    assert sorted(requested) == ["https://example.org/transcript/NM_000059.3",
                                 "https://example.org/transcript/NM_007294.3"]
    assert dp._get_transcript("NM_000059.3") == {"id": "NM_000059.3"}


def test_prefetch_batch_false_uses_concurrent(dp, monkeypatch):
    monkeypatch.setattr(dp, "_post_to_url",
                        lambda url, json_body: pytest.fail("batch=False must not POST"))
    monkeypatch.setattr(dp, "_get_from_url", lambda url: {"id": url.rsplit("/", 1)[1]})
    n = dp.prefetch(["NM_000059.3"], batch=False)
    assert n == 1
    assert dp.transcripts["NM_000059.3"] == {"id": "NM_000059.3"}


def test_prefetch_concurrent_failure_stores_none(dp, monkeypatch):
    monkeypatch.setattr(dp, "_post_to_url", lambda url, json_body: None)

    def fake_get(url):
        raise ConnectionError("boom")

    monkeypatch.setattr(dp, "_get_from_url", fake_get)
    n = dp.prefetch(["NM_000059.3"])  # must not raise
    assert n == 1
    assert dp.transcripts["NM_000059.3"] is None


# ---------------------------------------------------------------------------
# prefetch() - server request-size limits (chunking + 413/414 auto-shrink)
# ---------------------------------------------------------------------------

def test_prefetch_chunks_to_respect_server_limit(dp, monkeypatch):
    # 25 accessions with batch_size=10 -> 3 POSTs (10 + 10 + 5), all cached.
    accs = [f"NM_{i:06d}.1" for i in range(25)]
    sizes = []

    def fake_post(url, json_body):
        ids = json_body["ids"]
        sizes.append(len(ids))
        return {ac: {"id": ac} for ac in ids}

    monkeypatch.setattr(dp, "_post_to_url", fake_post)
    n = dp.prefetch(accs, batch_size=10)
    assert sizes == [10, 10, 5]      # chunked, not one giant request
    assert n == 25
    assert dp._get_transcript("NM_000024.1") == {"id": "NM_000024.1"}


def test_prefetch_halves_chunk_on_payload_too_large(dp, monkeypatch):
    # Server rejects >4 ids with 413; client must auto-shrink and still complete.
    import requests

    def fake_post(url, json_body):
        ids = json_body["ids"]
        if len(ids) > 4:
            resp = requests.Response()
            resp.status_code = 413
            raise requests.HTTPError("Payload Too Large", response=resp)
        return {ac: {"id": ac} for ac in ids}

    monkeypatch.setattr(dp, "_post_to_url", fake_post)
    accs = [f"NM_{i:06d}.1" for i in range(10)]
    n = dp.prefetch(accs, batch_size=8)   # 8 -> 413 -> shrink to 4 -> succeeds
    assert n == 10
    assert all(dp.transcripts[ac] == {"id": ac} for ac in accs)


def test_prefetch_chunked_missing_endpoint_falls_back_for_remainder(dp, monkeypatch):
    # First chunk works, then the endpoint "disappears" (None) -> concurrent fallback for the rest.
    accs = sorted(f"NM_{i:06d}.1" for i in range(15))
    calls = {"post": 0}

    def fake_post(url, json_body):
        calls["post"] += 1
        if calls["post"] == 1:
            return {ac: {"id": ac} for ac in json_body["ids"]}
        return None  # endpoint gone

    monkeypatch.setattr(dp, "_post_to_url", fake_post)
    monkeypatch.setattr(dp, "_get_from_url", lambda url: {"id": url.rsplit("/", 1)[1]})
    n = dp.prefetch(accs, batch_size=10)
    assert n == 15  # 10 via batch + 5 via concurrent singles
    assert all(dp.transcripts[ac] == {"id": ac} for ac in accs)


# ---------------------------------------------------------------------------
# prefetch_from_hgvs()
# ---------------------------------------------------------------------------

def test_prefetch_from_hgvs_extracts_accessions(dp, monkeypatch):
    captured = {}

    def fake_prefetch(accessions, **kwargs):
        captured["accessions"] = set(accessions)
        return len(accessions)

    monkeypatch.setattr(dp, "prefetch", fake_prefetch)

    dp.prefetch_from_hgvs(["NM_000059.3:c.1A>G",
                           "NM_007294.3:c.100del",
                           "NM_000059.3:c.2A>G"])  # duplicate accession
    assert captured["accessions"] == {"NM_000059.3", "NM_007294.3"}


def test_prefetch_from_hgvs_cleans_messy_input(dp, monkeypatch):
    captured = {}
    monkeypatch.setattr(dp, "prefetch",
                        lambda accessions, **kw: captured.setdefault("accessions", set(accessions)))
    # lowercase transcript + whitespace - clean_hgvs() should normalise to a usable accession
    dp.prefetch_from_hgvs([" nm_000059.4:c.316+5G>A "])
    assert captured["accessions"] == {"NM_000059.4"}


def test_prefetch_from_hgvs_strips_gene_parens(dp, monkeypatch):
    captured = {}
    monkeypatch.setattr(dp, "prefetch",
                        lambda accessions, **kw: captured.setdefault("accessions", set(accessions)))
    dp.prefetch_from_hgvs(["NM_000059.4(BRCA2):c.1A>G"])
    assert captured["accessions"] == {"NM_000059.4"}


def test_prefetch_from_hgvs_skips_non_transcripts(dp, monkeypatch):
    captured = {}
    monkeypatch.setattr(dp, "prefetch",
                        lambda accessions, **kw: captured.setdefault("accessions", set(accessions)))
    dp.prefetch_from_hgvs(["BRCA2:c.1A>G",            # gene-only - not a transcript
                           "NC_000017.11:g.43000A>G",  # genomic contig - not a transcript
                           "NM_000059.4:c.1A>G"])       # the only real transcript
    assert captured["accessions"] == {"NM_000059.4"}


def test_prefetch_from_hgvs_versionless(dp, monkeypatch):
    captured = {}
    monkeypatch.setattr(dp, "prefetch",
                        lambda accessions, **kw: captured.setdefault("accessions", set(accessions)))
    dp.prefetch_from_hgvs(["NM_000059:c.1A>G"])  # no version - kept as-is for server expansion
    assert captured["accessions"] == {"NM_000059"}


# ---------------------------------------------------------------------------
# _accession_from_hgvs (the pure extractor)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("hgvs_string,expected", [
    ("NM_000059.3:c.1A>G", "NM_000059.3"),
    ("NM_000059.4(BRCA2):c.1A>G", "NM_000059.4"),
    ("ENST00000380152.7:c.1A>G", "ENST00000380152.7"),
    ("NM_000059:c.1A>G", "NM_000059"),       # versionless kept
    ("BRCA2:c.1A>G", None),                    # gene-only
    ("NC_000017.11:g.43000A>G", None),         # genomic contig
    ("NM_000059.3", None),                      # no allele / no colon
])
def test_accession_from_hgvs(hgvs_string, expected):
    assert RESTDataProvider._accession_from_hgvs(hgvs_string) == expected
