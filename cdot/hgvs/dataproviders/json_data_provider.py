import abc
import gzip
import json
import os
import requests

from collections import defaultdict
from lazy import lazy
from hgvs.dataproviders.interface import Interface
from hgvs.dataproviders.seqfetcher import SeqFetcher
from intervaltree import IntervalTree
from typing import List

from cdot.assembly_helper import get_ac_name_map

class AbstractJSONDataProvider(Interface):
    NCBI_ALN_METHOD = "splign"
    required_version = "1.1"

    def __init__(self, assemblies: List[str] = None, mode=None, cache=None, seqfetcher=None):
        """ assemblies: defaults to ["GRCh37", "GRCh38"]
            seqfetcher defaults to biocommons SeqFetcher()
        """
        if assemblies is None:
            assemblies = ["GRCh37", "GRCh38"]

        super().__init__(mode=mode, cache=cache)
        if seqfetcher:
            try:
                seqfetcher.set_data_provider(self)
            except AttributeError:
                pass
        else:
            seqfetcher = SeqFetcher()
        self.seqfetcher = seqfetcher
        self.assembly_maps = {}
        for assembly_name in assemblies:
            self.assembly_maps[assembly_name] = get_ac_name_map(assembly_name)
        self.assembly_by_contig = {}
        for assembly_name, contig_map in self.assembly_maps.items():
            self.assembly_by_contig.update({contig: assembly_name for contig in contig_map.keys()})

    @abc.abstractmethod
    def _get_transcript(self, tx_ac):
        pass

    @abc.abstractmethod
    def _get_gene(self, gene):
        pass

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


class LocalDataProvider(AbstractJSONDataProvider):
    """ For JSON and Redis providers (implemented in cdot_rest)
        https://github.com/SACGF/cdot_rest - cdot_rest.redis_data_provider.RedisDataProvider """

    @abc.abstractmethod
    def _get_transcript_ids_for_gene(self, gene):
        pass

    @abc.abstractmethod
    def _get_contig_interval_tree(self, alt_ac):
        pass

    def get_tx_for_gene(self, gene):
        """ return transcript info records for supplied gene, in order of decreasing length """

        tx_list = []  # Store in tuples with length, so we can sort before returning
        for transcript_id in self._get_transcript_ids_for_gene(gene):
            transcript_data = self._get_transcript(transcript_id)
            cds_start_i = transcript_data.get("start_codon")
            cds_end_i = transcript_data.get("stop_codon")
            for build_data in transcript_data["genome_builds"].values():
                contig, tx_start, tx_end, _ = self._get_contig_start_end_strand(build_data)
                length = tx_end - tx_start
                tx_data = {
                    "hgnc": gene,
                    "cds_start_i": cds_start_i,
                    "cds_end_i": cds_end_i,
                    "tx_ac": transcript_id,
                    "alt_ac": contig,
                    "alt_aln_method": self.NCBI_ALN_METHOD,
                }
                tx_list.append((length, tx_data))

        return [x[1] for x in sorted(tx_list, key=lambda x: x[0], reverse=True)]

    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        """ return transcripts that overlap given region """
        tx_list = []
        if alt_aln_method == self.NCBI_ALN_METHOD:
            contig_iv_tree = self._get_contig_interval_tree(alt_ac)
            for interval in contig_iv_tree[start_i:end_i+1]:
                transcript_id = interval.data
                transcript_data = self._get_transcript(transcript_id)
                build_data = self._get_transcript_coordinates_for_contig(transcript_data, alt_ac)
                contig, tx_start, tx_end, strand = self._get_contig_start_end_strand(build_data)
                if contig == alt_ac:
                    tx_list.append({
                        "alt_ac": alt_ac,
                        "alt_aln_method": self.NCBI_ALN_METHOD,
                        "alt_strand": strand,
                        "start_i": tx_start,
                        "end_i": tx_end,
                        "tx_ac": transcript_id,
                    })
        return tx_list

    @staticmethod
    def _get_tx_by_gene_and_intervals(transcript_iter_items):
        # The region query works on exons, but storing all of these makes the interval tree huge
        # So we just store the start/end of each transcript ID, then look up the exons at retrieval time

        tx_by_gene = defaultdict(set)
        tx_intervals = defaultdict(IntervalTree)
        for transcript_id, transcript_data in transcript_iter_items:
            if gene_name := transcript_data['gene_name']:
                tx_by_gene[gene_name].add(transcript_id)

            for build_data in transcript_data["genome_builds"].values():
                contig = build_data["contig"]
                tx_start = build_data["exons"][0][0]
                tx_end = build_data["exons"][-1][1]
                tx_intervals[contig][tx_start:tx_end] = transcript_id
        return tx_by_gene, tx_intervals


class JSONDataProvider(LocalDataProvider):
    """ Local JSON file """
    def __init__(self, file_or_filename_list, mode=None, cache=None, seqfetcher=None):
        assemblies = set()
        self.transcripts = {}
        self.genes = {}
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
            if genes := data.get("genes"):
                for g in genes.values():
                    if gene_symbol := g.get("gene_symbol"):
                        self.genes[gene_symbol] = g
            self.cdot_data_version = tuple(int(v) for v in data["cdot_version"].split("."))

        super().__init__(assemblies=assemblies, mode=mode, cache=cache, seqfetcher=seqfetcher)

    def _get_transcript(self, tx_ac):
        return self.transcripts.get(tx_ac)

    def _get_gene(self, gene):
        return self.genes.get(gene)

    def _get_transcript_ids_for_gene(self, gene):
        tx_by_gene, _ = self._tx_by_gene_and_intervals
        return tx_by_gene[gene]

    def _get_contig_interval_tree(self, alt_ac):
        _, tx_intervals = self._tx_by_gene_and_intervals
        return tx_intervals[alt_ac]

    def get_pro_ac_for_tx_ac(self, tx_ac):
        if self.cdot_data_version < (0, 2, 8):
            cdot_version = '.'.join(str(v) for v in self.cdot_data_version)
            msg = f"ProteinID not in your JSON data version '{cdot_version}'. " \
                  "Please use data generated from cdot >= 0.2.8"
            raise NotImplementedError(msg)
        return super().get_pro_ac_for_tx_ac(tx_ac)

    @lazy
    def _tx_by_gene_and_intervals(self):
        return self._get_tx_by_gene_and_intervals(self.transcripts.items())

    def get_gene_info(self, gene):
        if self.cdot_data_version < (0, 2, 10):
            cdot_version = '.'.join(str(v) for v in self.cdot_data_version)
            msg = f"Gene Info not in your JSON data version '{cdot_version}'. " \
                  "Please use data generated from cdot >= 0.2.10"
            raise NotImplementedError(msg)
        return super().get_gene_info(gene)


class RESTDataProvider(AbstractJSONDataProvider):

    def __init__(self, url=None, secure=True, mode=None, cache=None, seqfetcher=None):
        assemblies = ["GRCh37", "GRCh38"]
        super().__init__(assemblies=assemblies, mode=mode, cache=cache, seqfetcher=seqfetcher)
        if url is None:
            if secure:
                url = "https://cdot.cc"
            else:
                url = "http://cdot.cc"
        self.url = url
        self.transcripts = {}
        self.genes = {}

    def _get_from_url(self, url):
        data = None
        response = requests.get(url)
        if response.ok:
            if 'application/json' in response.headers.get('Content-Type'):
                data = response.json()
            else:
                raise ValueError("Non-json response received for '%s' - are you behind a firewall?" % transcript_url)
        return data

    def _get_transcript(self, tx_ac):
        # We store None for 404 on REST
        if tx_ac in self.transcripts:
            return self.transcripts[tx_ac]

        transcript = self._get_from_url(self.url + "/transcript/" + tx_ac)
        self.transcripts[tx_ac] = transcript
        return transcript

    def _get_gene(self, gene_name):
        # We store None for 404 on REST
        if gene_name in self.genes:
            return self.genes[gene_name]

        gene = self._get_from_url(self.url + "/gene/" + gene_name)
        self.genes[gene_name] = gene
        return gene

    def get_tx_for_gene(self, gene_name):
        tx_list = []
        if data := self._get_from_url(f"{self.url}/transcripts/gene/{gene_name}"):
            tx_list = data["results"]
        return tx_list

    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        tx_list = []
        url = f"{self.url}/transcripts/region/{alt_ac}/{alt_aln_method}/{start_i}/{end_i}"
        if data := self._get_from_url(url):
            tx_list = data["results"]
        return tx_list
