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

    def __init__(self, mode=None, cache=None):
        super().__init__(mode=mode, cache=cache)
        self.seqfetcher = SeqFetcher()

    @abc.abstractmethod
    def _get_transcript(self, tx_ac):
        pass

    def data_version(self):
        return self.required_version

    def schema_version(self):
        return self.required_version

    def get_acs_for_protein_seq(self, seq):
        raise NotImplementedError("JSON data provider doesn't support get_acs_for_protein_seq")

    def get_assembly_map(self, assembly_name):
        """return a list of accessions for the specified assembly name (e.g., GRCh38.p5) """
        return make_ac_name_map(assembly_name)

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
        transcript = self.json_data.get(tx_ac)
        if not transcript:
            return None

        tx_exons = []  # Genomic order
        exons = transcript["cdna_match"]  # PyHGVS Needs to change
        alt_strand = 1 if transcript["strand"] == "+" else -1

        for (exon_id, (alt_start_i, alt_end_i, cds_start_i, cds_end_i, gap)) in enumerate(exons):
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
            # print(exon_data.values())
            tx_exons.append(exon_data)

        tx_exons.sort(key=lambda ex: ex["alt_start_i"])  # Sort by genomic order

#         print([d.values() for d in reversed(sorted(tx_exons, key=lambda e: e["ord"]))])
        return tx_exons

    def get_tx_for_gene(self, gene):
        # TODO: We could do this via JSON genes data
        raise NotImplementedError("TODO: get_tx_for_gene")

    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        # TODO: We could do this via HTSeq etc
        raise NotImplementedError("TODO")

    def get_tx_identity_info(self, tx_ac):
        transcript = self.json_data.get(tx_ac)
        if not transcript:
            return None

        tx_info = self.get_tx_info(tx_ac, tx_ac, "transcript")
        exons = transcript["cdna_match"]  # TODO: won't exist in all current PyReference etc
        tx_info["lengths"] = [end+1 - start for (_, _, start, end, _) in exons]  # Stranded order
        return tx_info

    def get_tx_info(self, tx_ac, alt_ac, alt_aln_method):
        transcript = self.json_data.get(tx_ac)
        if not transcript:
            return None

        hgnc = transcript["gene_name"]
        cds_start_i = transcript["start_codon_transcript_pos"]
        cds_end_i = transcript["stop_codon_transcript_pos"]
        return {
            "hgnc": hgnc,
            "cds_start_i": cds_start_i,
            "cds_end_i": cds_end_i,
            "tx_ac": tx_ac,
            "alt_ac": alt_ac,
            "alt_aln_method": alt_aln_method,
        }

    def get_tx_mapping_options(self, tx_ac):
        mapping_options = []
        transcript = self.json_data.get(tx_ac)
        if transcript:
            mo = {
                "tx_ac": tx_ac,
                "alt_ac": transcript["chrom"],
                "alt_aln_method": self.NCBI_ALN_METHOD,
            }
            mapping_options.append(mo)
        return mapping_options


class JSONDataProvider(AbstractJSONDataProvider):
    def _get_transcript(self, tx_ac):
        return self.json_data["transcripts"][tx_ac]

    def __init__(self, file_or_filename=None, json_data=None, mode=None, cache=None):
        super().__init__(mode=mode, cache=cache)
        if json_data:
            self.json_data = json_data
        elif file_or_filename:
            if isinstance(file_or_filename, str):
                if file_or_filename.endswith(".gz"):
                    f = gzip.open(file_or_filename)
                else:
                    f = open(file_or_filename)
            else:
                f = file_or_filename
            self.json_data = json.load(f)
        else:
            raise ValueError("Must pass either 'file_or_filename' or 'json_data'")
