import abc
import gzip
import json
import os

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
    def _get_transcript_for_assembly(self, tx_ac, assembly):
        pass

    def _get_transcript(self, tx_ac, alt_ac):
        assembly = self.assembly_by_contig.get(alt_ac)
        if assembly is None:
            supported_assemblies = ", ".join(self.assembly_maps.keys())
            raise ValueError(f"Contig '{alt_ac}' not supported. Supported assemblies: {supported_assemblies}")
        return self._get_transcript_for_assembly(tx_ac, assembly)

    def data_version(self):
        return self.required_version

    def schema_version(self):
        return self.required_version

    def get_acs_for_protein_seq(self, seq):
        raise NotImplementedError("JSON data provider doesn't support get_acs_for_protein_seq")

    def get_assembly_map(self, assembly_name):
        """return a list of accessions for the specified assembly name (e.g., GRCh38.p5) """
        assembly_map = self.assembly_maps.get(assembly_name)
        if assembly_map is None:
            raise ValueError(f"Assembly '{assembly_name}' not supported.")
        return assembly_map

    def get_gene_info(self, gene):
        raise NotImplementedError("JSON data provider doesn't support get_gene_info")

    def get_pro_ac_for_tx_ac(self, tx_ac):
        raise NotImplementedError("JSON data provider doesn't support get_pro_ac_for_tx_ac")

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

    def get_similar_transcripts(self, tx_ac):
        raise NotImplementedError("JSON data provider doesn't support get_similar_transcripts")

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
        transcript = self._get_transcript(tx_ac, alt_ac)
        if not transcript:
            return None

        tx_exons = []  # Genomic order
        exons = transcript["exons"]
        alt_strand = 1 if transcript["strand"] == "+" else -1

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

    def get_tx_for_gene(self, gene):
        # TODO: We could do this via JSON genes data
        raise NotImplementedError("TODO: get_tx_for_gene")

    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        # TODO: We could do this via HTSeq etc
        raise NotImplementedError("TODO")

    def get_tx_identity_info(self, tx_ac):
        # Get any transcript as it's assembly independent
        transcript = None
        for assembly in self.assembly_maps:
            transcript = self._get_transcript_for_assembly(tx_ac, assembly)
            if transcript:
                break

        if not transcript:
            return None

        tx_info = self._get_transcript_info(transcript)
        exons = transcript["exons"]
        stranded_order_exons = sorted(exons, key=lambda e: e[2])  # sort by exon_id
        tx_info["lengths"] = [ex[4] + 1 - ex[3] for ex in stranded_order_exons]
        tx_info["tx_ac"] = tx_ac
        tx_info["alt_ac"] = tx_ac  # Same again
        tx_info["alt_aln_method"] = "transcript"
        return tx_info

    @staticmethod
    def _get_transcript_info(transcript):
        gene_name = transcript["gene_name"]
        cds_start_i = transcript["start_codon"]
        cds_end_i = transcript["stop_codon"]
        return {
            "hgnc": gene_name,
            "cds_start_i": cds_start_i,
            "cds_end_i": cds_end_i,
        }

    def get_tx_info(self, tx_ac, alt_ac, alt_aln_method):
        transcript = self._get_transcript(tx_ac, alt_ac)
        if not transcript:
            return None

        tx_info = self._get_transcript_info(transcript)
        tx_info["tx_ac"] = tx_ac
        tx_info["alt_ac"] = alt_ac
        tx_info["alt_aln_method"] = alt_aln_method
        return tx_info

    def get_tx_mapping_options(self, tx_ac):
        mapping_options = []
        for assembly in self.assembly_maps:
            if transcript := self._get_transcript_for_assembly(tx_ac, assembly):
                mo = {
                    "tx_ac": tx_ac,
                    "alt_ac": transcript["contig"],
                    "alt_aln_method": self.NCBI_ALN_METHOD,
                }
                mapping_options.append(mo)
        return mapping_options


class JSONDataProvider(AbstractJSONDataProvider):
    def _get_transcript_for_assembly(self, tx_ac, assembly):
        return self._assembly_json[assembly].get(tx_ac)

    def __init__(self, assembly_json, mode=None, cache=None):
        super().__init__(assemblies=assembly_json.keys(), mode=mode, cache=cache)

        # No point being lazy as HGVS calls get_tx_mapping_options which requires us to check all builds
        self._assembly_json = {}
        for assembly, file_or_filename in assembly_json.items():
            if isinstance(file_or_filename, str):
                if file_or_filename.endswith(".gz"):
                    f = gzip.open(file_or_filename)
                else:
                    f = open(file_or_filename)
            else:
                f = file_or_filename
            data = json.load(f)
            # TODO: Check version is ok?
            self._assembly_json[assembly] = data["transcripts"]


class RESTDataProvider(AbstractJSONDataProvider):
    def _get_transcript_for_assembly(self, tx_ac, assembly):
        pass

    def __init__(self, url=None, mode=None, cache=None):
        if url is None:
            url = "http://cdot.cc"
        assemblies = ["GRCh37", "GRCh38"]
        super().__init__(assemblies=assemblies, mode=mode, cache=cache)
