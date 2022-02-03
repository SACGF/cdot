import abc
import gzip
import json
import os
import requests

from bioutils.assemblies import make_ac_name_map
from hgvs.dataproviders.interface import Interface
from hgvs.dataproviders.seqfetcher import SeqFetcher


class AbstractJSONDataProvider(Interface):
    NCBI_ALN_METHOD = "splign"
    required_version = "1.1"

    def __init__(self, assemblies, mode=None, cache=None):
        super().__init__(mode=mode, cache=cache)
        self.seqfetcher = SeqFetcher()
        self.assembly_maps = {}
        for assembly_name in assemblies:
            self.assembly_maps[assembly_name] = make_ac_name_map(assembly_name)
        self.assembly_by_contig = {}
        for assembly_name, contig_map in self.assembly_maps.items():
            self.assembly_by_contig.update({contig: assembly_name for contig in contig_map.keys()})

    @abc.abstractmethod
    def _get_transcript(self, tx_ac):
        pass

    def _get_transcript_coordinates_for_contig(self, transcript, alt_ac):
        assembly = self.assembly_by_contig.get(alt_ac)
        if assembly is None:
            supported_assemblies = ", ".join(self.assembly_maps.keys())
            raise ValueError(f"Contig '{alt_ac}' not supported. Supported assemblies: {supported_assemblies}")

        return transcript["genome_builds"].get(assembly)

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

    @staticmethod
    def sequence_source():
        seqrepo_dir = os.environ.get("HGVS_SEQREPO_DIR")
        seqrepo_url = os.environ.get("HGVS_SEQREPO_URL")
        if seqrepo_dir:
            return seqrepo_dir
        elif seqrepo_url:
            return seqrepo_url
        else:
            return "seqfetcher"

    def get_seq(self, ac, start_i=None, end_i=None):
        return self.seqfetcher.fetch_seq(ac, start_i, end_i)

    @staticmethod
    def _convert_gap_to_cigar(gap):
        """
                gap = 'M196 I1 M61 I1 M181'
                CIGAR = '194=1D60=1D184='
        """

        # This has to/from sequences inverted, so insertion is a deletion
        OP_CONVERSION = {
            "M": "=",
            "I": "D",
            "D": "I",
        }

        cigar_ops = []
        for gap_op in gap.split():
            gap_code = gap_op[0]
            length = int(gap_op[1:])

            cigar_ops.append(str(length) + OP_CONVERSION[gap_code])

        return "".join(cigar_ops)

    def get_tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        transcript = self._get_transcript(tx_ac)
        if not transcript:
            return None

        assembly_coordinates = self._get_transcript_coordinates_for_contig(transcript, alt_ac)

        tx_exons = []  # Genomic order
        exons = assembly_coordinates["exons"]
        alt_strand = 1 if assembly_coordinates["strand"] == "+" else -1

        for (alt_start_i, alt_end_i, exon_id, cds_start_i, cds_end_i, gap) in exons:
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
        transcript = self._get_transcript(tx_ac)
        if not transcript:
            return None

        tx_info = self._get_transcript_info(transcript)
        tx_info["tx_ac"] = tx_ac
        tx_info["alt_ac"] = alt_ac
        tx_info["alt_aln_method"] = alt_aln_method
        return tx_info

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
        raise NotImplementedError()

    def get_gene_info(self, gene):
        raise NotImplementedError()

    def get_pro_ac_for_tx_ac(self, tx_ac):
        raise NotImplementedError()

    def get_similar_transcripts(self, tx_ac):
        raise NotImplementedError()

    def get_tx_for_gene(self, gene):
        # TODO: We could build this from JSON pretty easily
        raise NotImplementedError()

    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        # TODO: This would be hard to do in Redis but we could do via HTSeq
        raise NotImplementedError()


class JSONDataProvider(AbstractJSONDataProvider):
    def _get_transcript(self, tx_ac):
        return self.transcripts.get(tx_ac)

    def __init__(self, file_or_filename_list, mode=None, cache=None):
        assemblies = set()
        self.transcripts = {}
        for file_or_filename in file_or_filename_list:
            if isinstance(file_or_filename, str):
                if file_or_filename.endswith(".gz"):
                    f = gzip.open(file_or_filename)
                else:
                    f = open(file_or_filename)
            else:
                f = file_or_filename
            data = json.load(f)
            assemblies.update(data["genome_builds"])
            self.transcripts.update(data["transcripts"])
        super().__init__(assemblies=assemblies, mode=mode, cache=cache)


class RESTDataProvider(AbstractJSONDataProvider):

    def _get_transcript(self, tx_ac):
        # We store None for 404 on REST
        if tx_ac in self.transcripts:
            return self.transcripts[tx_ac]

        transcript_url = self.url + "/transcript/" + tx_ac
        response = requests.get(transcript_url)
        if response.ok:
            if 'application/json' in response.headers.get('Content-Type'):
                transcript = response.json()
            else:
                raise ValueError("Non-json response received for '%s' - are you behind a firewall?" % transcript_url)
        else:
            transcript = None
        self.transcripts[tx_ac] = transcript
        return transcript

    def __init__(self, url=None, secure=True, mode=None, cache=None):
        assemblies = ["GRCh37", "GRCh38"]
        super().__init__(assemblies=assemblies, mode=mode, cache=cache)
        if url is None:
            if secure:
                url = "https://cdot.cc"
            else:
                url = "http://cdot.cc"
        self.url = url
        self.transcripts = {}
