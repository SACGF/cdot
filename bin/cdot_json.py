#!/usr/bin/env python3

from __future__ import print_function, absolute_import

import gzip
import json
import re
from argparse import ArgumentParser
from collections import defaultdict, Counter
from csv import DictReader
from typing import Dict

import cdot
import ijson
from cdot.gff.gff_parser import GTFParser, GFF3Parser


def handle_args():
    parser = ArgumentParser(description='cdot JSON manipulation tools')

    subparsers = parser.add_subparsers(dest='subcommand', required=True, help='TODO sub-command help')
    parser_gtf = subparsers.add_parser("gtf_to_json", help="Convert GTF to JSON")
    parser_gtf.add_argument("gtf_filename", help="GTF to convert to JSON")

    parser_gff3 = subparsers.add_parser("gff3_to_json", help="Convert GFF3 to JSON")
    parser_gff3.add_argument("gff3_filename", help="GFF3 to convert to JSON")

    for p in [parser_gtf, parser_gff3]:
        p.add_argument("--discard-contigs-with-underscores", action='store_true', default=True)
        p.add_argument('--url', help='URL (source of GFF) to store in "reference_gtf.url"')

    parser_uta = subparsers.add_parser("uta_to_json", help="Convert UTA to JSON")
    parser_uta.add_argument("uta_csv_filename", help="UTA SQL CSV to convert to JSON")
    parser_uta.add_argument('--url', help='UTA URL to store in "reference_gtf.url"')

    parser_historical = subparsers.add_parser("merge_historical", help="Merge multiple JSON files (keeping latest)")
    parser_historical.add_argument('json_filenames', nargs="+", action="extend",
                                   help='cdot JSON.gz - list OLDEST to NEWEST (newest is kept)')
    parser_historical.add_argument("--store-genes", action='store_true', help="Also store gene/version information")
    parser_historical.add_argument('--genome-build', required=True, help="'GRCh37' or 'GRCh38'")

    parser_builds = subparsers.add_parser("combine_builds", help="Merge multiple JSON files from different genome builds")

    parser_builds.add_argument('--grch37', required=True, help='cdot JSON.gz for GRCh37')
    parser_builds.add_argument('--grch38', required=True, help='cdot JSON.gz for GRCh38')

    # I want this to be subcommands rather than global (would need to be listed before subcommand)
    for p in [parser_gtf, parser_gff3, parser_uta, parser_historical, parser_builds]:
        p.add_argument('--output', required=True, help='Output filename')

    args = parser.parse_args()
    return args


class SortedSetEncoder(json.JSONEncoder):
    """ Dump set as list, from: https://stackoverflow.com/a/8230505/295724 """

    def default(self, obj):
        if isinstance(obj, set):
            return list(sorted(obj))
        return json.JSONEncoder.default(self, obj)


def gtf_to_json(args):
    parser = GTFParser(args.gtf_filename, args.discard_contigs_with_underscores)
    write_json(args.output, parser.get_data(), args.url)


def gff3_to_json(args):
    parser = GFF3Parser(args.gff3_filename, args.discard_contigs_with_underscores)
    write_json(args.output, parser.get_data(), args.url)


def uta_to_json(args):
    genes_by_id = {}
    transcripts_by_id = {}

    for data in DictReader(open(args.uta_csv_filename)):
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

    print("Writing UTA to cdot JSON.gz")
    data = {
        "genes_by_id": genes_by_id,
        "transcripts_by_id": transcripts_by_id,
    }
    write_json(args.output, data, args.url)


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

    cigar_pattern = re.compile(r"(\d+)([" + "".join(OP_CONVERSION.keys()) + "])")
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


def write_json(output_filename, data, url=None):
    if url:
        reference_gtf = data.get("reference_gtf", {})
        reference_gtf["url"] = url
        data["reference_gtf"] = reference_gtf

    data["version"] = cdot.get_json_schema_version()

    with gzip.open(output_filename, 'w') as outfile:
        json_str = json.dumps(data, cls=SortedSetEncoder, sort_keys=True)  # Sort so diffs work
        outfile.write(json_str.encode('ascii'))

    print("Wrote:", output_filename)


def merge_historical(args):
    print("merge_historical")
    TRANSCRIPT_FIELDS = ["biotype", "start_codon", "stop_codon"]
    GENOME_BUILD_FIELDS = ["cds_start", "cds_end", "strand", "contig", "exons"]

    gene_versions = {}  # We only keep those that are in the latest transcript version
    transcript_versions = {}

    for filename in args.json_filenames:
        print(f"Loading '{filename}'")
        with gzip.open(filename) as f:
            reference_gtf = next(ijson.items(f, "reference_gtf"))
            url = reference_gtf["url"]

            # PyReference stores transcripts under genes, while PyReference only has transcripts (that contain genes)
            transcript_gene_version = {}

            f.seek(0)  # Reset for next ijson call
            for gene_id, gene in ijson.kvitems(f, "genes_by_id"):
                if version := gene.get("version"):
                    gene_accession = f"{gene_id}.{version}"
                else:
                    gene_accession = gene_id

                gene_version = convert_gene_pyreference_to_gene_version_data(gene)
                gene_version["url"] = url
                gene_versions[gene_accession] = gene_version

                for transcript_accession in gene["transcripts"]:
                    transcript_gene_version[transcript_accession] = gene_accession

            f.seek(0)  # Reset for next ijson call
            for transcript_accession, pyreference_transcript_version in ijson.kvitems(f, "transcripts_by_id"):
                gene_accession = transcript_gene_version[transcript_accession]
                gene_version = gene_versions[gene_accession]

                transcript_version = {
                    "id": transcript_accession,
                    "gene_name": gene_version["gene_symbol"],
                }
                for field in TRANSCRIPT_FIELDS:
                    value = pyreference_transcript_version.get(field)
                    if value is not None:
                        transcript_version[field] = value

                if not gene_accession.startswith("_"):
                    transcript_version["gene_version"] = gene_accession
                if hgnc := gene_version.get("hgnc"):
                    transcript_version["hgnc"] = hgnc

                genome_build_coordinates = {
                    "url": url,
                }
                for field in GENOME_BUILD_FIELDS:
                    value = pyreference_transcript_version.get(field)
                    if value is not None:
                        genome_build_coordinates[field] = value
                transcript_version["genome_builds"] = {args.genome_build: genome_build_coordinates}
                transcript_versions[transcript_accession] = transcript_version

    # Summarise where it's from
    transcript_urls = Counter()
    for tv in transcript_versions.values():
        for build_coordinates in tv["genome_builds"].values():
            transcript_urls[build_coordinates["url"]] += 1

    total = sum(transcript_urls.values())
    print(f"{total} transcript versions from:")
    for url, count in transcript_urls.most_common():
        print(f"{url}: {count} ({count*100 / total:.1f}%)")

    print("Writing cdot data")
    with gzip.open(args.output, 'w') as outfile:
        data = {
            "transcripts": transcript_versions,
            "cdot_version": cdot.__version__,
            "genome_builds": [args.genome_build],
        }
        if args.store_genes:
            data["genes"] = gene_versions

        json_str = json.dumps(data)
        outfile.write(json_str.encode('ascii'))


def convert_gene_pyreference_to_gene_version_data(gene_data: Dict) -> Dict:
    gene_version_data = {
        'description': gene_data.get("description"),
        'gene_symbol': gene_data["name"],
    }

    if biotype_list := gene_data.get("biotype"):
        gene_version_data['biotype'] = ",".join(biotype_list)

    if hgnc_str := gene_data.get("HGNC"):
        # Has HGNC: (5 characters) at start of it
        gene_version_data["hgnc"] = hgnc_str[5:]

    # Only Ensembl Genes have versions
    if version := gene_data.get("version"):
        gene_data["version"] = version

    return gene_version_data


def combine_builds(args):
    print("combine_builds")
    genome_build_data = {
        "GRCh37": json.load(gzip.open(args.grch37)),
        "GRCh38": json.load(gzip.open(args.grch38)),
    }

    all_transcript_ids = set()
    for genome_build, data in genome_build_data.items():
        # TODO: Check cdot versions
        json_builds = data["genome_builds"]
        if json_builds != [genome_build]:
            raise ValueError(f"JSON file provided for {genome_build} needs to have only {genome_build} data (has {json_builds})")
        all_transcript_ids.update(data["transcripts"].keys())

    urls_different_coding = defaultdict(list)
    transcripts = {}
    for transcript_id in all_transcript_ids:
        transcript_data = {}  # Latest transcript data, "genome_builds" will be overwritten at end
        genome_builds = {}
        for genome_build, data in genome_build_data.items():
            if build_transcript := data["transcripts"].get(transcript_id):
                if transcript_data:
                    # Latest always used, but check existing - if codons are different old versions are wrong so remove
                    for field in ["start_codon", "stop_codon"]:
                        old = transcript_data.get(field)
                        new = build_transcript.get(field)
                        if old != new:  # Old relied on different codons so is obsolete
                            for build_coordinates in genome_builds.values():
                                url = build_coordinates["url"]
                                urls_different_coding[url].append(transcript_id)
                            genome_builds = {}

                transcript_data = build_transcript
                genome_builds[genome_build] = build_transcript["genome_builds"][genome_build]
        transcript_data["genome_builds"] = genome_builds
        transcripts[transcript_id] = transcript_data

    print("Writing cdot data")
    with gzip.open(args.output, 'w') as outfile:
        data = {
            "transcripts": transcripts,
            "cdot_version": cdot.__version__,
            "genome_builds": list(genome_build_data.keys()),
        }
        json_str = json.dumps(data)
        outfile.write(json_str.encode('ascii'))

    if urls_different_coding:
        print("Some transcripts were removed as they had different coding coordinates from latest")
        for url, transcript_ids in urls_different_coding.items():
            print(f"{url}: {','.join(transcript_ids)}")


def main():
    args = handle_args()
    subcommands = {
        "gtf_to_json": gtf_to_json,
        "gff3_to_json": gff3_to_json,
        "merge_historical": merge_historical,
        "combine_builds": combine_builds,
        "uta_to_json": uta_to_json,
    }
    subcommands[args.subcommand](args)


if __name__ == '__main__':
    main()
