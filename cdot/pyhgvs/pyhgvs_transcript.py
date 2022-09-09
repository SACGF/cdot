import abc
import gzip
import json
from importlib import metadata

import requests

from pyhgvs.utils import make_transcript


class AbstractPyHGVSTranscriptFactory(abc.ABC):
    GENERATE_CDS_START_END_KEY = "generate_cds_start_end"

    def __init__(self):
        pass

    @abc.abstractmethod
    def _get_transcript(self, transcript_id):
        pass

    def get_transcript_grch37(self, transcript_id, sacgf_pyhgvs_fork=False):
        return self.get_transcript(transcript_id, "GRCh37", sacgf_pyhgvs_fork=sacgf_pyhgvs_fork)

    def get_transcript_grch38(self, transcript_id, sacgf_pyhgvs_fork=False):
        return self.get_transcript(transcript_id, "GRCh38", sacgf_pyhgvs_fork=sacgf_pyhgvs_fork)

    def get_transcript(self, transcript_id, genome_build, sacgf_pyhgvs_fork=False):
        transcript = None
        if pyhgvs_data := self.get_pyhgvs_data(transcript_id, genome_build, sacgf_pyhgvs_fork=sacgf_pyhgvs_fork):
            transcript = make_transcript(pyhgvs_data)
            if self.GENERATE_CDS_START_END_KEY in pyhgvs_data:
                self._replace_cds_start_end(transcript)
        return transcript

    def get_pyhgvs_data(self, transcript_id, genome_build, sacgf_pyhgvs_fork=False):
        transcript_json = self._get_transcript(transcript_id)
        build_coords = transcript_json["genome_builds"].get(genome_build)
        if build_coords is None:
            return {}

        exons = build_coords['exons']
        start = exons[0][0]
        end = exons[-1][1]

        pyhgvs_data = {
            "id": transcript_json["id"],
            "chrom": build_coords['contig'],
            "start": start,
            "end": end,
            "strand": build_coords["strand"],
            # PyHGVS has cds_start/cds_end equal start/end if non-coding
            "cds_start": build_coords.get('cds_start', start),
            "cds_end": build_coords.get('cds_end', end),
            "gene_name": transcript_json['gene_name'],
        }

        if sacgf_pyhgvs_fork:
            # Remove the 3rd element (exon_number)
            exons = [e[:2] + e[3:] for e in exons]
            pyhgvs_data["cdna_match"] = exons
            pyhgvs_data["start_codon_transcript_pos"] = transcript_json.get("start_codon")
            pyhgvs_data["stop_codon_transcript_pos"] = transcript_json.get("stop_codon")
            if other_chroms := build_coords.get("other_chroms"):
                pyhgvs_data["other_chroms"] = other_chroms
        else:
            # Standard PyHGVS - only keep start/end
            exons = [e[:2] for e in exons]

        pyhgvs_data["exons"] = exons

        # UTA data doesn't have cds_start/cds_end which we need for PyHGVS
        if "start_codon" in transcript_json and "cds_start" not in build_coords:
            pyhgvs_data[self.GENERATE_CDS_START_END_KEY] = True

        return pyhgvs_data

    @staticmethod
    def _replace_cds_start_end(transcript):
        """ Replace cds_start/cds_end with values generated from start/stop codons
            This is because some data sources (eg UTA) do not provide cds_start/cds_end """

        print("FIXING _replace_cds_start_end")
        start_codon = CDNACoord(coord=transcript._start_codon_transcript_pos)
        stop_codon = CDNACoord(coord=transcript._stop_codon_transcript_pos)
        cds_start = transcript.cdna_to_genomic_coord(start_codon)
        cds_end = transcript.cdna_to_genomic_coord(stop_codon)
        # In pyhgvs they are always in genomic order
        if transcript.strand == '-':
            (cds_start, cds_end) = (cds_end, cds_start)
        transcript.chrom_start = cds_start
        transcript.chrom_stop = cds_stop


class PyHGVSTranscriptFactory(AbstractPyHGVSTranscriptFactory):
    def _get_transcript(self, transcript_id):
        return self.transcripts.get(transcript_id)

    def __init__(self, transcripts):
        super().__init__()
        self.transcripts = transcripts


class JSONPyHGVSTranscriptFactory(PyHGVSTranscriptFactory):
    def __init__(self, file_or_filename_list):
        transcripts = {}
        for file_or_filename in file_or_filename_list:
            if isinstance(file_or_filename, str):
                if file_or_filename.endswith(".gz"):
                    f = gzip.open(file_or_filename)
                else:
                    f = open(file_or_filename)
            else:
                f = file_or_filename
            data = json.load(f)
            transcripts.update(data["transcripts"])
        super().__init__(transcripts=transcripts)


class RESTPyHGVSTranscriptFactory(AbstractPyHGVSTranscriptFactory):

    def _get_transcript(self, transcript_id):
        # We store None for 404 on REST
        if transcript_id in self.transcripts:
            return self.transcripts[transcript_id]

        transcript_url = self.url + "/transcript/" + transcript_id
        response = requests.get(transcript_url)
        if response.ok:
            if 'application/json' in response.headers.get('Content-Type'):
                transcript = response.json()
            else:
                raise ValueError("Non-json response received for '%s' - are you behind a firewall?" % transcript_url)
        else:
            transcript = None
        self.transcripts[transcript_id] = transcript
        return transcript

    def __init__(self, url=None, secure=True):
        super().__init__()
        if url is None:
            if secure:
                url = "https://cdot.cc"
            else:
                url = "http://cdot.cc"
        self.url = url
        self.transcripts = {}


def is_sacgf_pyhgvs_fork():
    required_version = (0, 12, 0)  # Bumped version on 24 Nov 2021 - has mito and cDNA_match fixes
    imported_version = [int(v) for v in metadata.version("pyhgvs").split(".")]
    return tuple(imported_version) >= required_version


# Changes from old loading:

# See dot has no cds_start/end if non-coding
#  PyHGVS expects cds_start/cds_end be equal to end/end for non-coding transcripts (so coding length ie end-start = 0)
#     cds_start = transcript_data.get("cds_start", end)
#     cds_end = transcript_data.get("cds_end", end)


# VG loader also expects biotype to be comma sep, now is list
