import json
import logging
import os
import re
from collections import defaultdict

import requests
from hgvs.dataproviders.seqfetcher import SeqFetcher
from hgvs.exceptions import HGVSDataNotAvailableError

from cdot.hgvs.dataproviders import get_ac_name_map, get_name_ac_map, ExonsFromGenomeFastaSeqFetcher, \
    GenomeFastaSeqFetcher
from cdot.hgvs.dataproviders.seqfetcher import AbstractTranscriptSeqFetcher, PrefixSeqFetcher, \
    VerifyMultipleSeqFetcher, AlwaysFailSeqFetcher
from hgvs.dataproviders.interface import Interface

_REFSEQ_PREFIXES = {"NM_", "NR_"}

class _EnsemblTarkTranscriptSeqFetcher(AbstractTranscriptSeqFetcher):
    """ This retrieves sequences from Tark (but not genome/RefSeq checks etc) """
    def _get_transcript_seq(self, ac):
        if ac.startswith("NC_"):
            raise HGVSDataNotAvailableError()

        return self.hdp.get_transcript_sequence(ac)

    @property
    def source(self):
        return f"EnsemblTarkTranscriptSeqFetcher: hdp={self.hdp}"


class EnsemblTarkSeqFetcher(PrefixSeqFetcher):
    """ Default for EnsemblTarkDataProvider
        You may need to instantiate your own copy to provide fasta_files """
    def __init__(self, *args, fasta_files=None):
        super().__init__()
        tark_seqfetcher = _EnsemblTarkTranscriptSeqFetcher()
        if fasta_files is not None:
            fasta_seqfetcher = GenomeFastaSeqFetcher(*fasta_files)

            # For RefSeq - Tark doesn't have alignments, so can't check whether there are gaps
            # So we'll look up the genome reference, and if they don't match, throw an error
            class NoValidationExonsFromGenomeFastaSeqFetcher(ExonsFromGenomeFastaSeqFetcher):
                def get_mapping_options(self, ac):
                    # Normal 'get_tx_mapping_options' has a check that causes recursion
                    return self.hdp.get_tx_mapping_options_without_validation(ac)

            exons_seqfetcher = NoValidationExonsFromGenomeFastaSeqFetcher(*fasta_files, cache=True)
            refseq_seqfetcher = VerifyMultipleSeqFetcher(tark_seqfetcher, exons_seqfetcher)
        else:
            fasta_seqfetcher = SeqFetcher()  # Default HGVS
            msg = "You need to provide 'fasta_files' to use RefSeq transcripts. RefSeq transcripts can align with " + \
                   "gaps, so need to compare transcript/genome sequences"
            refseq_seqfetcher = AlwaysFailSeqFetcher(msg)

        self.prefix_seqfetchers.update({
            "NC_": fasta_seqfetcher,
            "ENST": tark_seqfetcher,
        })
        self.prefix_seqfetchers.update({rp: refseq_seqfetcher for rp in _REFSEQ_PREFIXES})


class EnsemblTarkDataProvider(Interface):
    """
        Tark - Transcript Archive - https://tark.ensembl.org/
    """
    NCBI_ALN_METHOD = "splign"
    required_version = "1.1"

    def __init__(self, assemblies: list[str] = None, mode=None, cache=None, seqfetcher=None):
        """ assemblies: defaults to ["GRCh37", "GRCh38"]
            seqfetcher defaults to EnsemblTarkSeqFetcher(), which uses hgvs SeqFetcher for fastas
        """
        self.base_url = "https://tark.ensembl.org/api"
        # Local caches
        self.transcript_results = {}

        if assemblies is None:
            assemblies = ["GRCh37", "GRCh38"]

        super().__init__(mode=mode, cache=cache)
        if seqfetcher is None:
            seqfetcher = EnsemblTarkSeqFetcher()

        try:
            seqfetcher.set_data_provider(self)
        except AttributeError:
            pass
        self.seqfetcher = seqfetcher
        self.assembly_maps = {}
        self.name_to_assembly_maps = {}

        for assembly_name in assemblies:
            self.assembly_maps[assembly_name] = get_ac_name_map(assembly_name)
            self.name_to_assembly_maps[assembly_name] = get_name_ac_map(assembly_name)

        self.assembly_by_contig = {}
        for assembly_name, contig_map in self.assembly_maps.items():
            self.assembly_by_contig.update({contig: assembly_name for contig in contig_map.keys()})

    def _get_from_url(self, url):
        response = requests.get(url)
        if response.ok:
            if 'application/json' in response.headers.get('Content-Type'):
                return response.json()
            else:
                raise ValueError("Non-json response received for '%s' - are you behind a firewall?" % url)
        response.raise_for_status()

    def _get_all_paginated_transcript_results(self, url):
        results = []
        while url:
            # Next links are http
            url = re.sub(r'^http://', 'https://', url)
            data = self._get_from_url(url)
            results.extend(data["results"])
            url = data["next"]
        return self._filter_dupes_take_most_recent(results)

    @staticmethod
    def _filter_dupes_take_most_recent(results):
        # There can be multiple results for a transcript accession / Genome build, eg:
        # https://tark.ensembl.org/web/transcript_details/NM_015120.4/NM_015120.4/?assembly_name=GRCh38
        # We want to pick the one with the most recent 'transcript_release_set'
        build_transcripts = defaultdict(lambda: defaultdict(list))
        for r in results:
            genome_build = EnsemblTarkDataProvider._get_genome_build(r)
            transcript_accession = EnsemblTarkDataProvider._get_transcript_accession(r)
            build_transcripts[genome_build][transcript_accession].append(r)

        filtered_results = []
        for genome_build, transcript_results in build_transcripts.items():
            for results in transcript_results.values():
                if len(results) > 1:
                    def _get_most_recent_release_date(result):
                        trs = result["transcript_release_set"]
                        return sorted([t["release_date"] for t in trs], reverse=True)[0]

                    results = sorted(results, key=_get_most_recent_release_date)
                filtered_results.append(results[0])
        return filtered_results

    @staticmethod
    def _get_transcript_id_and_version(transcript_accession: str):
        parts = transcript_accession.split(".")
        if len(parts) == 2:
            identifier = str(parts[0])
            version = int(parts[1])
        else:
            identifier, version = transcript_accession, None
        return identifier, version

    def _get_transcript_results(self, tx_ac):
        """ This can be a list of (1 per build) """

        # We store None for 404 on REST - so return even if false
        if tx_ac in self.transcript_results:
            return self.transcript_results[tx_ac]

        url = os.path.join(self.base_url, "transcript/?")
        stable_id, version = self._get_transcript_id_and_version(tx_ac)
        params = {
            "stable_id": stable_id,
            "stable_id_version": version,
            "expand_all": "true",
        }
        url += "&".join([f"{k}={v}" for k, v in params.items()])
        if results := self._get_all_paginated_transcript_results(url):
            if len(results) >= 1:
                self.transcript_results[tx_ac] = results
                return results
        raise HGVSDataNotAvailableError(f"Data for transcript='{tx_ac}' did not contain 'results': {results}")

    def _get_transcript_for_contig(self, transcript_results, alt_ac):
        assembly = self.assembly_by_contig.get(alt_ac)
        if assembly is None:
            supported_assemblies = ", ".join(self.assembly_maps.keys())
            raise ValueError(f"Contig '{alt_ac}' not supported. Supported assemblies: {supported_assemblies}")

        for transcript in transcript_results:
            genome_build = self._get_genome_build(transcript)
            if genome_build == assembly:
                return transcript
        return None

    def data_version(self):
        return self.required_version

    def schema_version(self):
        return self.required_version

    def get_assembly_map(self, assembly_name):
        """return a list of accessions for the specified assembly name (e.g., GRCh38.p5) """
        assembly_map = self.assembly_maps.get(assembly_name)
        if assembly_map is None:
            supported_assemblies = ", ".join(self.assembly_maps.keys())
            raise ValueError(f"Assembly '{assembly_name}' not supported. Supported assemblies: {supported_assemblies}")

        return assembly_map

    def sequence_source(self):
        return self.seqfetcher.source

    def get_seq(self, ac, start_i=None, end_i=None):
        return self.seqfetcher.fetch_seq(ac, start_i, end_i)

    def get_transcript_sequence(self, ac):
        seq = None
        if results := self._get_transcript_results(ac):
            transcript = results[0]  # any is fine
            seq = transcript["sequence"]["sequence"]
        return seq

    @staticmethod
    def _get_cds_start_end(transcript):
        # TODO: Could get this out of cds_info - three_prime_utr_length etc
        cds_start_i = None
        cds_end_i = None
        three_prime_utr_seq = transcript["three_prime_utr_seq"]
        five_prime_utr_seq = transcript["five_prime_utr_seq"]
        if three_prime_utr_seq and five_prime_utr_seq:
            sequence_length = sum([ex["loc_end"] - ex["loc_start"] + 1 for ex in transcript["exons"]])
            cds_start_i = len(five_prime_utr_seq)
            cds_end_i = sequence_length - len(three_prime_utr_seq)
        return cds_start_i, cds_end_i

    @staticmethod
    def _get_genome_build(transcript):
        assembly = transcript["assembly"]
        # Sometimes this can be expanded
        if isinstance(assembly, str):
            genome_build = assembly
        else:
            genome_build = assembly["assembly_name"]
        return genome_build

    def _get_transcript_contig(self, transcript):
        genome_build = self._get_genome_build(transcript)
        name_to_ac_map = self.name_to_assembly_maps[genome_build]
        name = transcript["loc_region"]
        return name_to_ac_map[name]

    @staticmethod
    def _get_transcript_accession(transcript):
        transcript_id = transcript["stable_id"]
        version = transcript["stable_id_version"]
        return f"{transcript_id}.{version}"

    def _get_chrom_from_contig(self, alt_ac):
        for ac_to_name_map in self.assembly_maps.values():
            if name := ac_to_name_map.get(alt_ac):
                return name
        raise ValueError(f"'{alt_ac}': unknown contig")

    def get_tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        self._check_alt_aln_method(alt_aln_method)
        transcript_results = self._get_transcript_results(tx_ac)
        if not transcript_results:
            return None

        transcript = self._get_transcript_for_contig(transcript_results, alt_ac)
        tx_exons = []  # Genomic order
        alt_strand = transcript["loc_strand"]

        tx_pos = 0
        for exon in transcript["exons"]:
            length = exon["loc_end"] - exon["loc_start"] + 1
            # UTA tx is 0 based
            tx_start_i = tx_pos
            tx_end_i = tx_pos + length
            tx_pos = tx_end_i

            # UTA is 0 based
            exon_data = {
                'tx_ac': tx_ac,
                'alt_ac': alt_ac,
                'alt_strand': alt_strand,
                'alt_aln_method': alt_aln_method,
                'ord': exon["exon_order"] - 1,  # tark is 1 based, UTA 0 based
                'tx_start_i': tx_start_i,
                'tx_end_i': tx_end_i,
                'alt_start_i': exon["loc_start"] - 1,
                'alt_end_i': exon["loc_end"],
                'cigar': str(length) + "=",  # Tark doesn't have alignment gaps
            }
            tx_exons.append(exon_data)

        # UTA wants exons in genomic order
        if alt_strand == -1:
            tx_exons.reverse()

        return tx_exons

    def get_tx_identity_info(self, tx_ac):
        # Get any transcript as it's assembly independent
        transcript_results = self._get_transcript_results(tx_ac)
        if not transcript_results:
            return None

        transcript = transcript_results[0]  # Any will do
        tx_info = self._get_transcript_info(transcript)

        # Only using lengths (same in each build) not coordinates so grab anything
        exons = transcript["exons"]
        tx_info["lengths"] = [ex["loc_end"] + 1 - ex["loc_start"] for ex in exons]
        tx_info["tx_ac"] = tx_ac
        tx_info["alt_ac"] = tx_ac  # Same again
        tx_info["alt_aln_method"] = "transcript"
        return tx_info

    @staticmethod
    def _get_transcript_info(transcript):
        gene_name = None
        if gene := transcript.get("genes"):
            gene_name = gene[0]["name"]
        cds_start_i, cds_end_i = EnsemblTarkDataProvider._get_cds_start_end(transcript)
        return {
            "hgnc": gene_name,
            "cds_start_i": cds_start_i,
            "cds_end_i": cds_end_i,
        }

    def get_tx_info(self, tx_ac, alt_ac, alt_aln_method):
        self._check_alt_aln_method(alt_aln_method)

        if transcript_results := self._get_transcript_results(tx_ac):
            if transcript := self._get_transcript_for_contig(transcript_results, alt_ac):
                tx_info = self._get_transcript_info(transcript)
                tx_info["tx_ac"] = tx_ac
                tx_info["alt_ac"] = alt_ac
                tx_info["alt_aln_method"] = self.NCBI_ALN_METHOD
                return tx_info

        raise HGVSDataNotAvailableError(
            f"No tx_info for (tx_ac={tx_ac},alt_ac={alt_ac},alt_aln_method={alt_aln_method})"
        )

    def get_tx_mapping_options_without_validation(self, tx_ac):
        # We need to be able to call this from NoValidationExonsFromGenomeFastaSeqFetcher
        mapping_options = []
        if transcript_results := self._get_transcript_results(tx_ac):
            for transcript in transcript_results:
                alt_ac = self._get_transcript_contig(transcript)
                mo = {
                    "tx_ac": tx_ac,
                    "alt_ac": alt_ac,
                    "alt_aln_method": self.NCBI_ALN_METHOD,
                }
                mapping_options.append(mo)
        return mapping_options

    def get_tx_mapping_options(self, tx_ac):
        try:
            self._verify_no_alignment_gaps(tx_ac)
            mapping_options = self.get_tx_mapping_options_without_validation(tx_ac)
        except HGVSDataNotAvailableError:
            logging.debug("'%s' transcript/genome sequence mismatch without alignment information - skipping.", tx_ac)
            mapping_options = []  # Can't map
        return mapping_options

    def get_acs_for_protein_seq(self, seq):
        """
            This is not implemented. The only caller has comment: 'TODO: drop get_acs_for_protein_seq'
            And is only ever called as a backup when get_pro_ac_for_tx_ac fails
        """
        return None

    def get_gene_info(self, gene):
        """
            This info is not available in TARK
            return {
                "hgnc": None,
                "maploc": None,
                "descr": None,
                "summary": None,
                "aliases": None, # UTA produces aliases that look like '{DCML,IMD21,MONOMAC,NFE1B}'
                "added": None,  # Don't know where this is stored/comes from (hgnc?)
            }
        """
        raise NotImplementedError()

    def get_pro_ac_for_tx_ac(self, tx_ac):
        pro_ac = None
        if transcript_results := self._get_transcript_results(tx_ac):
            for transcript in transcript_results:
                if translations := transcript.get("translations"):
                    t = translations[0]
                    protein_id = t["stable_id"]
                    version = t["stable_id_version"]
                    pro_ac = f"{protein_id}.{version}"
                    return pro_ac
        return None

    def get_similar_transcripts(self, tx_ac):
        """ UTA specific functionality that uses tx_similarity_v table
            This is not used by the HGVS library """
        raise NotImplementedError()

    def get_alignments_for_region(self, alt_ac, start_i, end_i, alt_aln_method=None):
        """ Prefer to use get_tx_for_region as this may be removed/deprecated
            This is never called externally, only used to implement get_tx_for_region in UTA data provider. """
        if alt_aln_method is None:
            alt_aln_method = self.NCBI_ALN_METHOD
        return self.get_tx_for_region(alt_ac, alt_aln_method, start_i, end_i)

    def _check_alt_aln_method(self, alt_aln_method):
        if alt_aln_method != self.NCBI_ALN_METHOD:
            raise HGVSDataNotAvailableError(f"cdot only supports alt_aln_method={self.NCBI_ALN_METHOD}")

    def _verify_no_alignment_gaps(self, tx_ac):
        # Tark doesn't have alignments, thus RefSeq could be wrong
        # See https://github.com/Ensembl/tark/issues/81
        for prefix in _REFSEQ_PREFIXES:
            if tx_ac.startswith(prefix):
                # see EnsemblTarkSeqFetcher - refseq gets seq from Tark and Exon fastas, die if not the same
                self.get_seq(tx_ac)
                break

    def get_tx_for_gene(self, gene):
        url = os.path.join(self.base_url, "transcript/search/?")
        params = {
            "identifier_field": gene,
            "expand": "exons,genes,sequence",
        }
        url += "&".join([f"{k}={v}" for k, v in params.items()])

        # TODO: We could cache the transcripts from these...
        tx_list = []
        if results := self._get_from_url(url):
            for transcript in results:
                cds_start_i, cds_end_i = self._get_cds_start_end(transcript)
                alt_ac = self._get_transcript_contig(transcript)

                tx_ac = self._get_transcript_accession(transcript)
                tx_list.append({
                    "hgnc": gene,
                    "cds_start_i": cds_start_i,
                    "cds_end_i": cds_end_i,
                    "tx_ac": tx_ac,
                    "alt_ac": alt_ac,
                    "alt_aln_method": "splign"
                })
        return tx_list


    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        assert end_i >= start_i, f"{end_i=} must be greater or equal than {start_i=}"

        url = os.path.join(self.base_url, "transcript/?")
        assembly = self.assembly_by_contig[alt_ac]
        loc_region = self._get_chrom_from_contig(alt_ac)
        params = {
            "assembly_name": assembly,  # Restrict to genome build
            # TODO: Expand all then store transcripts?
            "expand": "transcript_release_set",  # We need this to filter dupes
            "loc_end": end_i,
            "loc_region": loc_region,
            "loc_start": start_i + 1,  # UTA is 0 based, Tark is 1-based
        }
        url += "&".join([f"{k}={v}" for k, v in params.items()])
        tx_list = []
        if results := self._get_all_paginated_transcript_results(url):
            for transcript in results:
                contig = self._get_transcript_contig(transcript)
                tx_start = transcript["loc_start"] - 1
                tx_end = transcript["loc_end"]
                if contig == alt_ac:
                    tx_ac = self._get_transcript_accession(transcript)
                    tx_list.append({
                        "alt_ac": alt_ac,
                        "alt_aln_method": self.NCBI_ALN_METHOD,
                        "alt_strand": transcript["loc_strand"],
                        "start_i": tx_start,
                        "end_i": tx_end,
                        "tx_ac": tx_ac,
                    })
        return tx_list
