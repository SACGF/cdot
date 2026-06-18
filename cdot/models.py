"""
Typed representations of the cdot JSON data structures (issue #37).

``JSONDataProvider`` decodes cdot JSON into these ``msgspec`` structs rather than
plain dicts: it's ~2.5x faster to parse and ~20% smaller in memory on a full
RefSeq file. The structs are dict-read-compatible (``obj["key"]`` / ``obj.get()`` /
positional exon access), so this is invisible to consumers - the public HGVS
interface returns the same output as before. They also double as typed,
self-documenting access for anyone consuming cdot JSON directly.

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
from typing import Any, Dict, List, Optional, Union

import msgspec


class _DictAccessStruct(msgspec.Struct):
    """Base for the object-keyed structs giving them read-only ``dict``-style access.

    cdot's data-provider code (and external subclasses like cdot_rest's Redis
    provider) reads transcript/gene/build data with ``obj["key"]`` / ``obj.get(...)``.
    Supporting that here lets those structs be a drop-in for the plain dicts the
    code used to hold, so the typed storage stays invisible to consumers.
    """
    def __getitem__(self, key):
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError(key) from None

    def get(self, key, default=None):
        return getattr(self, key, default)


class Exon(msgspec.Struct, array_like=True, forbid_unknown_fields=False):
    """A single exon alignment.

    Serialised in JSON as a positional array, e.g.::

        [36552548, 36552986, 20, 2001, 2440, "M196 I1 M61 I1 M181"]

    Coordinates: ``alt_start``/``alt_end`` are 0-based genomic coordinates on the
    contig; ``cds_start``/``cds_end`` are 1-based transcript (cDNA) coordinates.
    ``exon_id`` is the ordinal of the exon in stranded (transcript) order.

    The original positional access (``exon[0]``, tuple unpacking) still works, so
    this is a drop-in for the plain list it replaces.
    """
    alt_start: int
    """0-based genomic start of the exon on the contig."""
    alt_end: int
    """0-based genomic end of the exon on the contig (exclusive)."""
    exon_id: int
    """Exon ordinal in stranded (transcript) order, starting at 0."""
    cds_start: int
    """1-based transcript (cDNA) start coordinate of the exon."""
    cds_end: int
    """1-based transcript (cDNA) end coordinate of the exon."""
    gap: Optional[str] = None
    """Alignment gap as a cdot 'gap' string (e.g. ``'M196 I1 M61'``) or ``None`` if the exon aligns cleanly."""

    def __getitem__(self, i):
        return getattr(self, self.__struct_fields__[i])

    def __len__(self):
        return len(self.__struct_fields__)

    def __iter__(self):
        return (getattr(self, f) for f in self.__struct_fields__)


class GenomeBuild(_DictAccessStruct, forbid_unknown_fields=False):
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
    """0-based genomic CDS start on the contig (coding transcripts only)."""
    cds_end: Optional[int] = None
    """0-based genomic CDS end on the contig (coding transcripts only)."""
    start: Optional[int] = None
    """0-based genomic start of the transcript on the contig."""
    stop: Optional[int] = None
    """0-based genomic end of the transcript on the contig."""
    tag: Optional[str] = None
    """Comma-separated tags (e.g. ``'MANE_Select,Ensembl_canonical'``); typically Ensembl only."""
    note: Optional[str] = None
    other_chroms: Optional[List[str]] = None
    """Other contigs this transcript also aligns to (e.g. PAR/alt loci)."""
    source: Optional[Union[str, List[str]]] = None
    """Annotation source (GTF/GFF column 2, e.g. ``'BestRefSeq'``); data schema >= 0.2.32.
    A single string in early 0.2.32 data, a list (e.g. ``['BestRefSeq']``) from 0.2.33 on."""
    ccds: Optional[str] = None
    """CCDS id, when present; data schema >= 0.2.33."""
    transcript_support_level: Optional[str] = None
    """Ensembl transcript support level (TSL); data schema >= 0.2.33."""


class Transcript(_DictAccessStruct, forbid_unknown_fields=False):
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
    source: Optional[List[str]] = None
    """Annotation source(s) this transcript came from (e.g. ``['NCBI']``)."""
    partial: Optional[int] = None
    """Non-zero if the transcript is annotated as partial/incomplete."""


class Gene(_DictAccessStruct, forbid_unknown_fields=False):
    """Gene-level metadata (present when the source provided gene info)."""
    gene_symbol: Optional[str] = None
    aliases: Optional[str] = None
    # gene biotype is a list[str] in current releases but a comma-separated bare str
    # in older data (<= 0.2.19; it became a list in 0.2.20 - see issue #111). The
    # data schema int wasn't bumped, so old and new releases are indistinguishable by
    # version - we accept both forms and normalise to list[str] in __post_init__.
    biotype: Optional[Union[List[str], str]] = None
    description: Optional[str] = None
    hgnc: Optional[str] = None
    map_location: Optional[str] = None
    summary: Optional[str] = None
    url: Optional[str] = None
    source: Optional[List[str]] = None
    transcripts: Optional[List[str]] = None
    """Transcript accessions belonging to this gene (when provided)."""

    def __post_init__(self):
        # Convert the legacy comma-separated str biotype (<= 0.2.19) to a list so
        # consumers always see list[str], regardless of which release produced the data.
        if isinstance(self.biotype, str):
            self.biotype = [b for b in self.biotype.split(",") if b]


class CdotData(msgspec.Struct, forbid_unknown_fields=False):
    """The top-level contents of a cdot JSON(.gz) file."""
    cdot_version: str
    genome_builds: List[str]
    transcripts: Dict[str, Transcript] = {}
    genes: Dict[str, Gene] = {}
    metadata: Optional[Dict[str, Any]] = None
    """Release provenance (input_files, method, sys.argv, url_counts, ...); freeform, present in merged release files."""


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
