""" Biocommons hgvs data providers and sequence fetchers backed by cdot data.

    Data providers (each implements :class:`hgvs.dataproviders.interface.Interface`):

    - :class:`JSONDataProvider`: cdot JSON(.gz) release files loaded into memory
    - :class:`RESTDataProvider`: a cdot REST server, https://cdotlib.org by default
    - :class:`EnsemblTarkDataProvider`: Ensembl's Tark transcript archive REST API, no cdot data
    - :class:`LocalDataProvider`: base class for providers that hold the data themselves

    Sequence fetchers, passed as ``seqfetcher=`` to a provider:

    - :class:`FastaSeqFetcher`: a local genome FASTA (needs pysam, ``pip install cdot[fasta]``)
    - :class:`ChainedSeqFetcher`: try several fetchers in turn, eg SeqRepo then FASTA
    - :class:`PrefixSeqFetcher`, :class:`VerifyMultipleSeqFetcher`, :class:`EnsemblTarkSeqFetcher`
"""
from .json_data_provider import (
    AbstractJSONDataProvider,
    JSONDataProvider,
    LocalDataProvider,
    RESTDataProvider,
    get_ac_name_map,
    get_name_ac_map,
)
from .seqfetcher import (
    AbstractTranscriptSeqFetcher,
    AlwaysFailSeqFetcher,
    ChainedSeqFetcher,
    MultiSeqFetcher,
    PrefixSeqFetcher,
    VerifyMultipleSeqFetcher,
)
from .ensembl_tark_data_provider import EnsemblTarkDataProvider, EnsemblTarkSeqFetcher

try:
    from .fasta_seqfetcher import (
        ExonsFromGenomeFastaSeqFetcher,
        FastaSeqFetcher,
        GenomeFastaSeqFetcher,
        GenomeMismatchError,
        GenomeMismatchPolicy,
        GenomeMismatchWarning,
        get_ensembl_tark_fasta_seqfetchers,
    )
    _FASTA_NAMES = [
        "FastaSeqFetcher",
        "GenomeFastaSeqFetcher",
        "ExonsFromGenomeFastaSeqFetcher",
        "GenomeMismatchPolicy",
        "GenomeMismatchWarning",
        "GenomeMismatchError",
        "get_ensembl_tark_fasta_seqfetchers",
    ]
except ImportError as e:
    if (e.name or "").partition(".")[0] != "pysam":
        raise

    class FastaSeqFetcher:
        """ Stand-in for :class:`cdot.hgvs.dataproviders.fasta_seqfetcher.FastaSeqFetcher` when pysam
            is not installed. Constructing it raises ImportError saying how to install it. """

        def __init__(self, *args, **kwargs):
            raise ImportError("FastaSeqFetcher needs pysam: pip install cdot[fasta]")

    _FASTA_NAMES = ["FastaSeqFetcher"]

__all__ = [
    "AbstractJSONDataProvider",
    "LocalDataProvider",
    "JSONDataProvider",
    "RESTDataProvider",
    "EnsemblTarkDataProvider",
    "EnsemblTarkSeqFetcher",
    "get_ac_name_map",
    "get_name_ac_map",
    "PrefixSeqFetcher",
    "MultiSeqFetcher",
    "ChainedSeqFetcher",
    "VerifyMultipleSeqFetcher",
    "AlwaysFailSeqFetcher",
    "AbstractTranscriptSeqFetcher",
] + _FASTA_NAMES
