import os

import requests
from hgvs.exceptions import HGVSDataNotAvailableError

from cdot.hgvs.dataproviders import get_ac_name_map, get_name_ac_map
from cdot.hgvs.dataproviders.seqfetcher import AbstractTranscriptSeqFetcher
from hgvs.dataproviders.interface import Interface


class EnsemblTarkTranscriptSeqFetcher(AbstractTranscriptSeqFetcher):
    def _get_transcript_seq(self, ac, alt_ac, alt_aln_method):
        pass


class EnsemblTarkDataProvider(Interface):
    """
        Tark - Transcript Archive - https://tark.ensembl.org/
    """
    NCBI_ALN_METHOD = "splign"
    required_version = "1.1"

    def __init__(self, assemblies: list[str] = None, mode=None, cache=None, seqfetcher=None):
        """ assemblies: defaults to ["GRCh37", "GRCh38"]
            seqfetcher defaults to biocommons SeqFetcher()
        """
        self.base_url = "https://tark.ensembl.org/api"
        # Local caches
        self.transcript_results = {}

        if assemblies is None:
            assemblies = ["GRCh37", "GRCh38"]

        super().__init__(mode=mode, cache=cache)
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
        data = None
        response = requests.get(url)
        if response.ok:
            if 'application/json' in response.headers.get('Content-Type'):
                data = response.json()
            else:
                raise ValueError("Non-json response received for '%s' - are you behind a firewall?" % url)
        return data

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
        data = self._get_from_url(url)
        if results := data["results"]:
            if len(results) >= 1:
                self.transcript_results[tx_ac] = results
                return results
        raise HGVSDataNotAvailableError(f"Data for transcript='{tx_ac}' did not contain 'results': {data}")

    def _get_transcript_for_contig(self, transcript_results, alt_ac):
        assembly = self.assembly_by_contig.get(alt_ac)
        if assembly is None:
            supported_assemblies = ", ".join(self.assembly_maps.keys())
            raise ValueError(f"Contig '{alt_ac}' not supported. Supported assemblies: {supported_assemblies}")

        for transcript in transcript_results:
            if transcript["assembly"]["assembly_name"] == assembly:
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

    def _get_transcript_contig(self, transcript):
        assembly = transcript["assembly"]
        # Sometimes this can be expanded
        if isinstance(assembly, str):
            genome_build = assembly
        else:
            genome_build = assembly["assembly_name"]
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
            length = exon["loc_end"] - exon["loc_start"]  # Will be same for both transcript/genomic
            tx_start_i = tx_pos + 1
            tx_end_i = tx_pos + length
            tx_pos += length

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
        gene_name = transcript["genes"][0]["name"]
        cds_start_i, cds_end_i = EnsemblTarkDataProvider._get_cds_start_end(transcript)
        return {
            "hgnc": gene_name,
            "cds_start_i": cds_start_i,
            "cds_end_i": cds_end_i,
        }

    def get_tx_info(self, tx_ac, alt_ac, alt_aln_method):
        self._check_alt_aln_method(alt_aln_method)

        if transcript_results := self._get_transcript_results(tx_ac):
            transcript = self._get_transcript_for_contig(transcript_results, alt_ac)
            tx_info = self._get_transcript_info(transcript)
            tx_info["tx_ac"] = tx_ac
            tx_info["alt_ac"] = alt_ac
            tx_info["alt_aln_method"] = self.NCBI_ALN_METHOD
            return tx_info

        raise HGVSDataNotAvailableError(
            f"No tx_info for (tx_ac={tx_ac},alt_ac={alt_ac},alt_aln_method={alt_aln_method})"
        )

    def get_tx_mapping_options(self, tx_ac):
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
        url = os.path.join(self.base_url, "transcript/?")
        # TODO: Off by 1 errors?
        # I think there is an issue for ensembl tark here about cross genomes
        loc_region = self._get_chrom_from_contig(alt_ac)
        params = {
            "loc_region": loc_region,
            "loc_start": start_i,
            "loc_end": end_i,
            # TODO: Expand all then store transcripts?
            # "expand_all": "false",
        }
        url += "&".join([f"{k}={v}" for k, v in params.items()])
        tx_list = []
        if results := self._get_from_url(url):
            for transcript in results:
                contig = self._get_transcript_contig(transcript)
                # TODO: Off by 1?
                tx_start = transcript["loc_start"]
                tx_end = transcript["loc_end"]
                if contig == alt_ac:
                    tx_ac = self._get_transcript_accession(transcript)
                    tx_list.append({
                        "alt_ac": alt_ac,
                        "alt_aln_method": self.NCBI_ALN_METHOD,
                        "alt_strand": transcript["strand"],
                        "start_i": tx_start,
                        "end_i": tx_end,
                        "tx_ac": tx_ac,
                    })
        return tx_list