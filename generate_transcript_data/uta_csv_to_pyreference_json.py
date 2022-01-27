#!/bin/env python3

import gzip
import json
import re
from argparse import ArgumentParser
from csv import DictReader


def handle_args():
    parser = ArgumentParser(description='Convert UTA csv (from SQL export) into cdot JSON')
    parser.add_argument('--uta-csv', required=True, help='UTA CSV export')
    parser.add_argument('--output', required=True, help='Output filename')
    parser.add_argument('--url', required=True, help='URL to store for UTA connection')
    return parser.parse_args()


def main():
    args = handle_args()
    genes_by_id = {}
    transcripts_by_id = {}

    for data in DictReader(open(args.uta_csv)):
        transcript_accession = data["ac"]
        gene_name = data["hgnc"]

        # PyReference expects a gene versions etc - we don't know these so will make a fake one.
        # Anything starting with "_" is not copied over
        fake_gene_version = "_" + gene_name
        gene = genes_by_id.get(fake_gene_version)
        if gene is None:
            genes_by_id[fake_gene_version] = {
                "name": gene_name,
                "transcripts": [transcript_accession],
            }
        elif transcript_accession not in gene["transcripts"]:
            gene["transcripts"].append(transcript_accession)

        transcript_version = {
            "contig": data["contig"],
            "strand": "+" if data["strand"] == 1 else "-",
            "gene_version": fake_gene_version,
            "exons": _convert_uta_exons(data["exon_starts"], data["exon_ends"], data["cigars"]),
        }
        cds_start_i = data["cds_start_i"]
        if cds_start_i:
            transcript_version["start_codon"] = int(cds_start_i)
            transcript_version["stop_codon"] = int(data["cds_end_i"])

        transcripts_by_id[transcript_accession] = transcript_version

    print("Writing UTA to pyreference JSON.gz")
    with gzip.open(args.output, 'w') as outfile:
        data = {
            "genes_by_id": genes_by_id,
            "transcripts_by_id": transcripts_by_id,
            "reference_gtf": {"url": args.url},
        }
        json_str = json.dumps(data)
        outfile.write(json_str.encode('ascii'))


def _convert_uta_exons(exon_starts, exon_ends, cigars):
    # UTA is output sorted in exon order (stranded)
    exon_starts = exon_starts.split(",")
    exon_ends = exon_ends.split(",")
    cigars = cigars.split(",")
    exons = []
    ex_ord = 0
    ex_transcript_start = 1  # transcript coords are 1 based
    for ex_start, ex_end, ex_cigar in zip(exon_starts, exon_ends, cigars):
        gap, exon_length = _cigar_to_gap_and_length(ex_cigar)
        ex_transcript_end = ex_transcript_start + exon_length - 1
        exons.append((int(ex_start), int(ex_end), ex_ord, ex_transcript_start, ex_transcript_end, gap))
        ex_transcript_start += exon_length
        ex_ord += 1

    exons.sort(key=lambda e: e[0])  # Genomic order
    return exons


def _cigar_to_gap_and_length(cigar):
    """
            gap = 'M196 I1 M61 I1 M181'
            CIGAR = '194=1D60=1D184='
    """

    # This has to/from sequences inverted, so insertion is a deletion
    OP_CONVERSION = {
        "=": "M",
        "D": "I",
        "I": "D",
        "X": "=",  # TODO: This is probably wrong! check if GTF gap has mismatch?
    }

    cigar_pattern = re.compile("(\d+)([" + "".join(OP_CONVERSION.keys()) + "])")
    gap_ops = []
    exon_length = 0
    for (length_str, cigar_code) in cigar_pattern.findall(cigar):
        exon_length += int(length_str)  # This shouldn't include one of the indels?
        gap_ops.append(OP_CONVERSION[cigar_code] + length_str)

    gap = " ".join(gap_ops)
    if len(gap_ops) == 1:
        gap_op = gap_ops[0]
        if gap_op.startswith("M"):
            gap = None

    return gap, exon_length


if __name__ == '__main__':
    main()

    # record looks like
    """
    {'ac': 'NM_000014.4',
     'hgnc': 'A2M',
     'origin_url': 'http://www.ncbi.nlm.nih.gov/refseq/',
     'contig': 'NC_000012.11',
     'strand': '-1',
     'cds_start_i': '113',
     'cds_end_i': '4538',
     'exon_starts': '9268359,9265955,9264972,9264754,9262909,9262462,9261916,9260119,9259086,9258831,9256834,9254042,9253739,9251976,9251202,9248134,9247568,9246060,9243796,9242951,9242497,9241795,9232689,9232234,9231839,9230296,9229941,9229351,9227155,9225248,9224954,9223083,9222340,9221335,9220778,9220303',
     'exon_ends': '9268558,9266139,9265132,9264807,9262930,9262631,9262001,9260240,9259201,9258941,9256996,9254270,9253803,9252119,9251352,9248296,9247680,9246175,9244025,9243078,9242619,9241847,9232773,9232411,9231927,9230453,9230016,9229532,9227379,9225467,9225082,9223174,9222409,9221438,9220820,9220435',
     'cigars': '199=,184=,160=,53=,21=,169=,85=,121=,115=,110=,162=,228=,64=,143=,150=,63=1X98=,112=,115=,229=,127=,122=,52=,84=,177=,88=,157=,75=,181=,224=,219=,128=,91=,69=,103=,42=,132='}
    """

    """
'NM_001001568.1': {'biotype': ['protein_coding'],
  'cds_end': 44195403,
  'cds_start': 44073976,
  'contig': 'NC_000021.8',
  'exons': [[44073861, 44073993, 0, 1, 132, None],
   [44119077, 44119121, 1, 133, 176, None],
   [44152179, 44152234, 2, 177, 231, None],
   [44153465, 44153536, 3, 232, 302, None],
   [44163885, 44163970, 4, 303, 387, None],
   [44171225, 44171307, 5, 388, 469, None],
   [44174099, 44174174, 6, 470, 544, None],
   [44179108, 44179195, 7, 545, 631, None],
   [44180437, 44180542, 8, 632, 736, None],
   [44180934, 44181017, 9, 737, 819, None],
   [44182192, 44182349, 10, 820, 976, None],
   [44185490, 44185604, 11, 977, 1090, None],
   [44188297, 44188402, 12, 1091, 1195, None],
   [44189136, 44189265, 13, 1196, 1324, None],
   [44190812, 44190908, 14, 1325, 1420, None],
   [44192548, 44192630, 15, 1421, 1502, None],
   [44195389, 44195619, 16, 1503, 1732, None]],
  'start': 44073861,
  'start_codon': 115,
  'stop': 44195619,
  'stop_codon': 1516,
  'strand': '+',
  'id': 'NM_001001568.1',
  'gene_version': '5152',
  'gene_name': 'PDE9A',
  'url': 'http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/105.20190906/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz',
  'hgnc': '8795'},
    """