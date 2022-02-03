import abc
import json

import requests

from pyhgvs.models import Transcript, Position, Exon


class AbstractPyHGVSTranscriptFactory(abc.ABC):
    def __init__(self):
        pass

    @abc.abstractmethod
    def _get_transcript(self, transcript_id):
        pass

    def get_transcript_grch37(self, transcript_id):
        return self.get_transcript(transcript_id, "GRCh37")

    def get_transcript_grch38(self, transcript_id):
        return self.get_transcript(transcript_id, "GRCh38")

    def get_transcript(self, transcript_id, genome_build):
        transcript_json = self._get_transcript(transcript_id)
        build_coords = transcript_json["genome_builds"].get(genome_build)
        if build_coords is None:
            return None

        transcript_name = transcript_json['id']
        if '.' in transcript_name:
            name, version = transcript_name.split('.')
        else:
            name, version = transcript_name, None

        contig = build_coords['contig']
        start = build_coords['exons'][0][0]
        end = build_coords['exons'][-1][1]
        is_forward_strand = build_coords['strand'] == '+'

        transcript = Transcript(
            name=name,
            version=int(version) if version is not None else None,
            gene=transcript_json['gene_name'],
            tx_position=Position(
                contig,
                start,
                end,
                is_forward_strand),
            cds_position=Position(
                contig,
                build_coords['cds_start'],
                build_coords['cds_end'],
                is_forward_strand),
        )

        stranded_exons = sorted([ex[:3] for ex in build_coords['exons']], key=lambda ex: ex[2])
        for (exon_start, exon_end, exon_number) in stranded_exons:
            transcript.exons.append(
                Exon(transcript=transcript,
                     tx_position=Position(
                         contig,
                         exon_start,
                         exon_end,
                         is_forward_strand),
                     exon_number=exon_number))

        return transcript


class JSONPyHGVSTranscriptFactory(AbstractPyHGVSTranscriptFactory):
    def _get_transcript(self, transcript_id):
        return self.transcripts.get(transcript_id)

    def __init__(self, file_or_filename_list):
        super().__init__()
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
            self.transcripts.update(data["transcripts"])


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


# TODO: Write code to load SACGF pyhgvs changes, ie gaps and start/stop codons etc

# Changes from old loading:

# See dot has no cds_start/end if non-coding
#  PyHGVS expects cds_start/cds_end be equal to end/end for non-coding transcripts (so coding length ie end-start = 0)
#     cds_start = transcript_data.get("cds_start", end)
#     cds_end = transcript_data.get("cds_end", end)



# VG loader also expects biotype to be comma sep, now is list