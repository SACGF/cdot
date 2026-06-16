import abc
import re

from pysam.libcfaidx import FastaFile
from hgvs.dataproviders.interface import Interface
from hgvs.exceptions import HGVSDataNotAvailableError
from bioutils.sequences import reverse_complement

from cdot.hgvs.dataproviders.seqfetcher import AbstractTranscriptSeqFetcher, PrefixSeqFetcher


class GenomeFastaSeqFetcher:
    def __init__(self, *args):
        self.source = "Local Fasta file reference"
        self.contig_fastas = {}
        for fasta_filename in args:
            fasta_file = FastaFile(fasta_filename)
            for contig in fasta_file.references:
                self.contig_fastas[contig] = fasta_file

        if not self.contig_fastas:
            raise ValueError("Need to provide at least one of fasta file as argument")

    def fetch_seq(self, ac, start_i=None, end_i=None):
        if fasta_file := self.contig_fastas.get(ac):  # Contig
            return fasta_file.fetch(ac, start_i, end_i).upper()

        raise HGVSDataNotAvailableError(f"Accession '{ac}' not in fasta contigs")


class ExonsFromGenomeFastaSeqFetcher(AbstractTranscriptSeqFetcher):
    """ This produces artificial transcript sequences by pasting together exons from the genome
        It is possible that this does not exactly match the transcript sequences - USE AT OWN RISK! """
    _CIGAR_PATTERN = re.compile(r"(\d+)([=DIX])")

    def __init__(self, *args, cache=True):
        self.cache = cache
        self.transcript_cache = {}
        self.hdp = None  # Set when passed to data provider (via set_data_provider)
        self.source = "Transcript Exons using Genome Fasta file reference"
        self.contig_fastas = {}
        for fasta_filename in args:
            fasta_file = FastaFile(fasta_filename)
            for contig in fasta_file.references:
                self.contig_fastas[contig] = fasta_file

        if not self.contig_fastas:
            raise ValueError("Need to provide at least one of fasta file as argument")
        super().__init__(*args, cache)

    def get_mapping_options(self, ac):
        return self.hdp.get_tx_mapping_options(ac)

    def _get_transcript_seq(self, ac):
        possible_contigs = set()
        for tx_mo in self.get_mapping_options(ac):
            alt_ac = tx_mo["alt_ac"]
            possible_contigs.add(alt_ac)
            if alt_ac in self.contig_fastas:
                return self._fetch_seq_from_fasta(ac, alt_ac, tx_mo["alt_aln_method"])

        msg = f"Failed to fetch {ac} from {self.source}. "
        if possible_contigs:
            possible_contigs = sorted(possible_contigs)
            raise HGVSDataNotAvailableError(f"{msg} No Fasta provided with contigs: {possible_contigs}")
        raise HGVSDataNotAvailableError(f"{msg} Transcript '{ac}' not found.")

    def _fetch_seq_from_fasta(self, ac, alt_ac, alt_aln_method):
        fasta_file = self.contig_fastas[alt_ac]

        exons = self.hdp.get_tx_exons(ac, alt_ac, alt_aln_method)
        exon_sequences = []
        expected_transcript_length = 0
        sorted_exons = list(sorted(exons, key=lambda ex: ex["ord"]))
        first_exon = sorted_exons[0]
        transcript_start_offset = first_exon["tx_start_i"]  # HGVS/UTA starts w/0
        if transcript_start_offset:
            exon_sequences.append("N" * transcript_start_offset)
            expected_transcript_length += transcript_start_offset

        for exon in sorted_exons:
            exon_seq = fasta_file.fetch(alt_ac, exon["alt_start_i"], exon["alt_end_i"])
            exon_seq = exon_seq.upper()

            exon_seq_list = []
            start = 0
            # We are using HGVS cigar
            for (length_str, op) in self._CIGAR_PATTERN.findall(exon["cigar"]):
                length = int(length_str)
                if op == 'D':    # Deletion in reference vs transcript
                    exon_seq_list.append("N" * length)
                    # Don't increment start (as we didn't move along genomic exon)
                elif op == 'I':  # Insertion in reference vs transcript
                    # Leave out of exon_seq
                    start += length  # We do increment through genomic sequence though
                else:  # match/mismatch
                    exon_seq_list.append(exon_seq[start:start+length])
                    start += length

            exon_seq = "".join(exon_seq_list)
            if exon["alt_strand"] == -1:
                exon_seq = reverse_complement(exon_seq)
            exon_sequences.append(exon_seq)
            expected_transcript_length += exon["tx_end_i"] - exon["tx_start_i"]

        transcript_sequence = "".join(exon_sequences)
        if len(transcript_sequence) != expected_transcript_length:
            raise ValueError(f"Error creating {ac} sequence from genome fasta ({alt_ac}): "
                             f"{expected_transcript_length=} != {len(transcript_sequence)=}")
        return transcript_sequence



def get_ensembl_tark_fasta_seqfetchers(*fasta_files, cache=True):
    """ Build ``(nc_seqfetcher, refseq_seqfetcher)`` for ``EnsemblTarkSeqFetcher`` from genome FASTA
        files. Requires pysam (``pip install cdot[fasta]``).

        Tark has no RefSeq alignments (so no gap info), so RefSeq sequences are built from the genome
        and verified against Tark - if they disagree (eg the transcript aligns with gaps) it fails
        rather than return a wrong sequence.

        Example::

            nc_sf, refseq_sf = get_ensembl_tark_fasta_seqfetchers("GRCh38.fna")
            seqfetcher = EnsemblTarkSeqFetcher(nc_seqfetcher=nc_sf, refseq_seqfetcher=refseq_sf)
            hdp = EnsemblTarkDataProvider(seqfetcher=seqfetcher)
    """
    # Imported lazily to keep the (pysam-free) ensembl_tark module importable without cdot[fasta]
    from cdot.hgvs.dataproviders.ensembl_tark_data_provider import _EnsemblTarkTranscriptSeqFetcher
    from cdot.hgvs.dataproviders.seqfetcher import VerifyMultipleSeqFetcher

    nc_seqfetcher = GenomeFastaSeqFetcher(*fasta_files)

    # Tark doesn't have RefSeq alignments, so can't check whether there are gaps - look up the
    # genome reference, and if it doesn't match Tark, throw an error
    class NoValidationExonsFromGenomeFastaSeqFetcher(ExonsFromGenomeFastaSeqFetcher):
        def get_mapping_options(self, ac):
            # Normal 'get_tx_mapping_options' has a check that causes recursion
            return self.hdp.get_tx_mapping_options_without_validation(ac)

    exons_seqfetcher = NoValidationExonsFromGenomeFastaSeqFetcher(*fasta_files, cache=cache)
    refseq_seqfetcher = VerifyMultipleSeqFetcher(_EnsemblTarkTranscriptSeqFetcher(), exons_seqfetcher)
    return nc_seqfetcher, refseq_seqfetcher


class FastaSeqFetcher(PrefixSeqFetcher):
    """ Re-implementing using above - deprecated use """

    def __init__(self, *args, cache=True):
        default_seqfetcher = ExonsFromGenomeFastaSeqFetcher(*args, cache=True)
        super().__init__(default_seqfetcher=default_seqfetcher)
        self.genome_fasta_seq_fetcher = GenomeFastaSeqFetcher(*args)
        self.prefix_seqfetchers.update({
            "NC_": self.genome_fasta_seq_fetcher,
        })

    @property
    def contig_fastas(self):
        """ For backwards compatibility """
        return self.genome_fasta_seq_fetcher.contig_fastas
