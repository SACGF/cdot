"""
Optional typed representations of the cdot JSON data structures.

These classes are an *optional* convenience layer for consumers who want typed,
self-documenting access to cdot JSON (issue #37). The JSON files remain the
canonical data format - nothing in cdot's core runtime depends on these classes,
and they require the optional ``msgspec`` dependency::

    pip install cdot[typed]

``msgspec`` is used (rather than dataclasses/attrs/pydantic) because it provides
fast, zero-boilerplate typed (de)serialisation and - via ``array_like`` structs -
lets us give the positional ``exons`` arrays real field names, which is most of
the documentation value being asked for.

Example::

    from cdot import models

    data = models.load("cdot-0.2.27.refseq.grch38.json.gz")
    tx = data.transcripts["NM_001637.3"]
    print(tx.gene_name, tx.protein)
    for build_name, build in tx.genome_builds.items():
        for exon in build.exons:
            print(exon.alt_start, exon.alt_end, exon.cds_start, exon.gap)
"""
from __future__ import annotations

import gzip
from typing import Dict, List, Optional

try:
    import msgspec
except ImportError as e:  # pragma: no cover - exercised only without the extra installed
    raise ImportError(
        "cdot.models requires the optional 'msgspec' dependency. "
        "Install it with: pip install cdot[typed]"
    ) from e


class Exon(msgspec.Struct, array_like=True, forbid_unknown_fields=False):
    """A single exon alignment.

    Serialised in JSON as a positional array, e.g.::

        [36552548, 36552986, 20, 2001, 2440, "M196 I1 M61 I1 M181"]

    Coordinates: ``alt_start``/``alt_end`` are 0-based genomic coordinates on the
    contig; ``cds_start``/``cds_end`` are 1-based transcript (cDNA) coordinates.
    ``exon_id`` is the ordinal of the exon in stranded (transcript) order.
    """
    alt_start: int
    alt_end: int
    exon_id: int
    cds_start: int
    cds_end: int
    gap: Optional[str] = None
    """Alignment gap as a cdot 'gap' string (e.g. ``'M196 I1 M61'``) or ``None`` if the exon aligns cleanly."""


class GenomeBuild(msgspec.Struct, forbid_unknown_fields=False):
    """A transcript's coordinates on one genome build (e.g. ``GRCh38``)."""
    contig: str
    """Chromosome/contig accession, e.g. ``'NC_000007.14'``."""
    strand: str
    """``'+'`` or ``'-'``."""
    exons: List[Exon]
    url: Optional[str] = None
    """Source annotation file the transcript was extracted from."""
    # cds_* / start / stop are absent for non-coding transcripts and some sources
    cds_start: Optional[int] = None
    cds_end: Optional[int] = None
    start: Optional[int] = None
    stop: Optional[int] = None
    tag: Optional[str] = None
    """Comma-separated tags (e.g. ``'MANE_Select,Ensembl_canonical'``); typically Ensembl only."""


class Transcript(msgspec.Struct, forbid_unknown_fields=False):
    """A single transcript and its per-build coordinates."""
    id: str
    """Transcript accession, e.g. ``'NM_001637.3'``."""
    genome_builds: Dict[str, GenomeBuild]
    """Build name (``'GRCh37'`` / ``'GRCh38'`` / ``'T2T-CHM13v2.0'``) -> coordinates."""
    gene_name: Optional[str] = None
    gene_version: Optional[str] = None
    biotype: Optional[List[str]] = None
    protein: Optional[str] = None
    """Protein accession (coding transcripts only)."""
    start_codon: Optional[int] = None
    """1-based transcript (cDNA) coordinate of the CDS start (coding only)."""
    stop_codon: Optional[int] = None
    """1-based transcript (cDNA) coordinate of the CDS end (coding only)."""
    hgnc: Optional[str] = None
    cdot: Optional[str] = None
    """cdot version that generated/last touched this transcript record."""


class Gene(msgspec.Struct, forbid_unknown_fields=False):
    """Gene-level metadata (present when the source provided gene info)."""
    gene_symbol: Optional[str] = None
    aliases: Optional[str] = None
    biotype: Optional[str] = None
    description: Optional[str] = None
    hgnc: Optional[str] = None
    map_location: Optional[str] = None
    summary: Optional[str] = None
    url: Optional[str] = None


class CdotData(msgspec.Struct, forbid_unknown_fields=False):
    """The top-level contents of a cdot JSON(.gz) file."""
    cdot_version: str
    genome_builds: List[str]
    transcripts: Dict[str, Transcript] = {}
    genes: Dict[str, Gene] = {}


def loads(data: bytes | str) -> CdotData:
    """Decode cdot JSON (``bytes`` or ``str``) into a typed :class:`CdotData`."""
    return msgspec.json.decode(data, type=CdotData)


def load(file_or_filename) -> CdotData:
    """Load a cdot JSON file into a typed :class:`CdotData`.

    Accepts a filename (``.json`` or ``.json.gz``) or an already-open binary
    file object. Mirrors ``JSONDataProvider``'s file handling.
    """
    if isinstance(file_or_filename, str):
        open_func = gzip.open if file_or_filename.endswith(".gz") else open
        with open_func(file_or_filename, "rb") as f:
            return loads(f.read())
    return loads(file_or_filename.read())


def transcript_from_dict(d: dict) -> Transcript:
    """Build a typed :class:`Transcript` from a plain cdot transcript dict.

    Useful for the on-demand path: keep the existing dict store and construct a
    typed object only when a consumer asks for one.
    """
    return msgspec.convert(d, Transcript)


def gene_from_dict(d: dict) -> Gene:
    """Build a typed :class:`Gene` from a plain cdot gene dict."""
    return msgspec.convert(d, Gene)
