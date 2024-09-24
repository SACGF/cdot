import os

import requests
from hgvs.exceptions import HGVSDataNotAvailableError

from cdot.hgvs.dataproviders import get_ac_name_map
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
        self.transcripts = {}
        self.genes = {}

        if assemblies is None:
            assemblies = ["GRCh37", "GRCh38"]

        super().__init__(mode=mode, cache=cache)
        self.assembly_maps = {}
        for assembly_name in assemblies:
            self.assembly_maps[assembly_name] = get_ac_name_map(assembly_name)
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
    def get_transcript_id_and_version(transcript_accession: str):
        parts = transcript_accession.split(".")
        if len(parts) == 2:
            identifier = str(parts[0])
            version = int(parts[1])
        else:
            identifier, version = transcript_accession, None
        return identifier, version

    def _get_transcript(self, tx_ac):
        # We store None for 404 on REST
        if tx_ac in self.transcripts:
            return self.transcripts[tx_ac]

        url = os.path.join(self.base_url, "transcript/?")
        stable_id, version = self.get_transcript_id_and_version(tx_ac)
        params = {
            "stable_id": stable_id,
            "stable_id_version": version,
            "expand_all": "true",
        }
        url += "&".join([f"{k}={v}" for k, v in params.items()])
        transcript = self._get_from_url(url)
        self.transcripts[tx_ac] = transcript
        return transcript

    def _get_gene(self, gene_name):
        # We store None for 404 on REST
        if gene_name in self.genes:
            return self.genes[gene_name]

        gene = self._get_from_url(self.url + "/gene/" + gene_name)
        self.genes[gene_name] = gene
        return gene

    def _get_transcript_coordinates_for_contig(self, transcript, alt_ac):
        assembly = self.assembly_by_contig.get(alt_ac)
        if assembly is None:
            supported_assemblies = ", ".join(self.assembly_maps.keys())
            raise ValueError(f"Contig '{alt_ac}' not supported. Supported assemblies: {supported_assemblies}")

        return transcript["genome_builds"].get(assembly)

    @staticmethod
    def _get_contig_start_end_strand(build_data):
        contig = build_data["contig"]
        strand = 1 if build_data["strand"] == "+" else -1
        start = build_data["exons"][0][0]
        end = build_data["exons"][-1][1]
        return contig, start, end, strand

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

    def get_transcript_seq(self, ac):
        transcript = self.get_transcript(ac)

    def get_tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        self._check_alt_aln_method(alt_aln_method)
        transcript = self._get_transcript(tx_ac)
        if not transcript:
            return None

        assembly_coordinates = self._get_transcript_coordinates_for_contig(transcript, alt_ac)

        tx_exons = []  # Genomic order
        exons = assembly_coordinates["exons"]
        alt_strand = 1 if assembly_coordinates["strand"] == "+" else -1

        for (alt_start_i, alt_end_i, exon_id, cds_start_i, cds_end_i, gap) in exons:
            # cdot tx_start_i/tx_end_i is 1 based while UTA is 0 based
            tx_start_i = cds_start_i - 1
            tx_end_i = cds_end_i
            if gap is not None:
                cigar = self._convert_gap_to_cigar(gap)
            else:
                length = alt_end_i - alt_start_i  # Will be same for both transcript/genomic
                cigar = str(length) + "="

            exon_data = {
                'tx_ac': tx_ac,
                'alt_ac': alt_ac,
                'alt_strand': alt_strand,
                'alt_aln_method': alt_aln_method,
                'ord': exon_id,
                'tx_start_i': tx_start_i,
                'tx_end_i': tx_end_i,
                'alt_start_i': alt_start_i,
                'alt_end_i': alt_end_i,
                'cigar': cigar,
            }
            tx_exons.append(exon_data)

        return tx_exons

    def get_tx_identity_info(self, tx_ac):
        # Get any transcript as it's assembly independent
        transcript = self._get_transcript(tx_ac)
        if not transcript:
            return None

        tx_info = self._get_transcript_info(transcript)

        # Only using lengths (same in each build) not coordinates so grab anything
        exons = []
        for build_coordinates in transcript["genome_builds"].values():
            exons = build_coordinates["exons"]
            break

        stranded_order_exons = sorted(exons, key=lambda e: e[2])  # sort by exon_id
        tx_info["lengths"] = [ex[4] + 1 - ex[3] for ex in stranded_order_exons]
        tx_info["tx_ac"] = tx_ac
        tx_info["alt_ac"] = tx_ac  # Same again
        tx_info["alt_aln_method"] = "transcript"
        return tx_info

    @staticmethod
    def _get_transcript_info(transcript):
        gene_name = transcript["gene_name"]
        cds_start_i = transcript.get("start_codon")
        cds_end_i = transcript.get("stop_codon")
        return {
            "hgnc": gene_name,
            "cds_start_i": cds_start_i,
            "cds_end_i": cds_end_i,
        }

    def get_tx_info(self, tx_ac, alt_ac, alt_aln_method):
        self._check_alt_aln_method(alt_aln_method)

        if transcript := self._get_transcript(tx_ac):
            for build_data in transcript["genome_builds"].values():
                if alt_ac == build_data["contig"]:
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
        if transcript := self._get_transcript(tx_ac):
            for build_coordinates in transcript["genome_builds"].values():
                mo = {
                    "tx_ac": tx_ac,
                    "alt_ac": build_coordinates["contig"],
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
        gene_info = None
        if g := self._get_gene(gene):
            # UTA produces aliases that look like '{DCML,IMD21,MONOMAC,NFE1B}'
            uta_style_aliases = '{' + g["aliases"].replace(" ", "") + '}'
            gene_info = {
                "hgnc": g["gene_symbol"],
                "maploc": g["map_location"],
                "descr": g["description"],
                "summary": g["summary"],
                "aliases": uta_style_aliases,
                "added": None,  # Don't know where this is stored/comes from (hgnc?)
            }
        return gene_info

    def get_pro_ac_for_tx_ac(self, tx_ac):
        pro_ac = None
        if transcript := self._get_transcript(tx_ac):
            pro_ac = transcript.get("protein")
        return pro_ac

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
        raise NotImplementedError()

    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        raise NotImplementedError()
