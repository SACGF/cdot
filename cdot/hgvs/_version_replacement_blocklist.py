"""Transcript-version re-placement blocklist (GENERATED - do not edit).

Regenerate with generate_transcript_data/generate_version_replacement_blocklist.py.
Source: cdot-0.2.33.refseq.GRCh38.json.gz (GRCh38).

Each (accession, requested_version, substitute_version) is a substitution whose
intrinsic CDS structure and alignment gaps match (so the structural safety check
would call it coordinate-safe) but whose genomic CDS placement differs, so it is
refused by is_version_substitution_safe even when the requested version is absent
from the loaded data (the case the provider cannot recompute).
"""
KNOWN_REPLACEMENT_SUBSTITUTIONS = frozenset({
    ('NM_000910', 2, 3),
    ('NM_000910', 2, 4),
    ('NM_000910', 3, 2),
    ('NM_000910', 4, 2),
    ('NM_001005567', 2, 3),
    ('NM_001005567', 3, 2),
    ('NM_001142273', 1, 2),
    ('NM_001142273', 2, 1),
    ('NM_001142274', 1, 2),
    ('NM_001142274', 2, 1),
    ('NM_001146181', 1, 2),
    ('NM_001146181', 1, 4),
    ('NM_001146181', 2, 1),
    ('NM_001146181', 4, 1),
    ('NM_001207051', 1, 2),
    ('NM_001207051', 2, 1),
    ('NM_001365429', 1, 2),
    ('NM_001365429', 2, 1),
    ('NM_004066', 1, 2),
    ('NM_004066', 1, 3),
    ('NM_004066', 2, 1),
    ('NM_004066', 3, 1),
    ('NM_015282', 2, 3),
    ('NM_015282', 3, 2),
    ('NM_018263', 4, 6),
    ('NM_018263', 5, 6),
    ('NM_018263', 6, 4),
    ('NM_018263', 6, 5),
})
