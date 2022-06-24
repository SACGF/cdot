#!/usr/bin/env python3

from __future__ import print_function, absolute_import

import gzip
import json
import re
import sys
from argparse import ArgumentParser
from collections import defaultdict, Counter
from csv import DictReader

import cdot
import ijson
from cdot.gff.gff_parser import GTFParser, GFF3Parser


def handle_args():
    parser = ArgumentParser(description='cdot JSON manipulation tools')
    parser.add_argument('--version', action='store_true', help='show version number')

    subparsers = parser.add_subparsers(dest='subcommand', help='TODO sub-command help')
    parser_gtf = subparsers.add_parser("gtf_to_json", help="Convert GTF to JSON")
    parser_gtf.add_argument("gtf_filename", help="GTF to convert to JSON")

    parser_gff3 = subparsers.add_parser("gff3_to_json", help="Convert GFF3 to JSON")
    parser_gff3.add_argument("gff3_filename", help="GFF3 to convert to JSON")

    for p in [parser_gtf, parser_gff3]:
        p.add_argument("--discard-contigs-with-underscores", action='store_true', default=True)
        p.add_argument('--url', required=True, help='URL (source of GFF) to store in "reference_gtf.url"')
        p.add_argument('--genome-build', required=True, help="'GRCh37' or 'GRCh38'")

    parser_uta = subparsers.add_parser("uta_to_json", help="Convert UTA to JSON")
    parser_uta.add_argument("uta_csv_filename", help="UTA SQL CSV to convert to JSON")
    parser_uta.add_argument('--url', required=True, help='UTA URL to store in "reference_gtf.url"')

    parser_historical = subparsers.add_parser("merge_historical", help="Merge multiple JSON files (keeping latest)")
    parser_historical.add_argument('json_filenames', nargs="+", action="extend",
                                   help='cdot JSON.gz - list OLDEST to NEWEST (newest is kept)')
    parser_historical.add_argument("--no-genes", action='store_true', help="Save <5% space by not storing gene/version information")
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
    parser = GTFParser(args.gtf_filename, args.genome_build, args.url,
                       discard_contigs_with_underscores=args.discard_contigs_with_underscores)
    genes, transcripts = parser.get_genes_and_transcripts()
    write_cdot_json(args.output, genes, transcripts, [args.genome_build])


def gff3_to_json(args):
    parser = GFF3Parser(args.gff3_filename, args.genome_build, args.url,
                        discard_contigs_with_underscores=args.discard_contigs_with_underscores)
    genes, transcripts = parser.get_genes_and_transcripts()
    write_cdot_json(args.output, genes, transcripts, [args.genome_build])


def uta_to_json(args):
    genes_by_id = {}
    transcripts_by_id = {}

    for data in DictReader(open(args.uta_csv_filename)):
        contig = data["contig"]
        if "," in contig:
            continue  # query returned both chrX/Y - only 8 of these so skip
        transcript_accession = data["ac"]
        if "/" in transcript_accession:
            continue  # transcript accession looks strange eg "NM_001202404.1/111..1620"
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
            "contig": contig,
            "strand": "+" if data["strand"] == "1" else "-",
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
    write_cdot_json(args.output, genes_by_id, transcripts_by_id, [args.genome_build])


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
        UTA performs their own alignment, and stores a CIGAR string of alignments

        We need to convert to the GFF3 gap attribute, see:

        But UTA/GFF3 has to/from sequences inverted, so insertion is a deletion

        Example input:  CIGAR = '194=1D60=1D184='
        Example output:   gap = 'M194 I1 M60 I1 M184'
    """

    OP_CONVERSION = {
        "=": "M",
        "D": "I",
        "I": "D",
        "X": "M",  # X=mismatch, GTF gap doesn't have mismatch, so count as match
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


def write_cdot_json(filename, genes, transcript_versions, genome_builds):
    print("Writing cdot data")
    data = {
        "cdot_version": cdot.__version__,
        "genome_builds": genome_builds,
        "transcripts": transcript_versions,
    }
    if genes:
        data["genes"] = genes

    with gzip.open(filename, 'wt') as outfile:
        json.dump(data, outfile, cls=SortedSetEncoder, sort_keys=True)  # Sort so diffs work


def merge_historical(args):
    """ Loads multiple JSON files, produces a file with the latest transcript version """
    gene_versions = {}  # We only keep those that are in the latest transcript version
    transcript_versions = {}
    gene_accessions_for_symbol = defaultdict(set)

    for filename in args.json_filenames:
        print(f"Loading '{filename}'")
        with gzip.open(filename) as f:
            for gene_accession, gene_version in ijson.kvitems(f, "genes"):
                gene_versions[gene_accession] = gene_version
                if not gene_accession.startswith("_"):
                    gene_accessions_for_symbol[gene_version["gene_symbol"]].add(gene_accession)

            f.seek(0)  # Reset for next ijson call
            for transcript_accession, historical_transcript_version in ijson.kvitems(f, "transcripts"):
                gene_accession = historical_transcript_version["gene_version"]
                gene_version = gene_versions[gene_accession]
                gene_symbol = gene_version["gene_symbol"]

                if gene_accession.startswith("_"):  # Not real - fake UTA so try and grab old one
                    fixed_ga = None
                    if previous_tv := transcript_versions.get(transcript_accession):
                        fixed_ga = previous_tv.get("gene_version")
                    if not fixed_ga:
                        if potential_ga := gene_accessions_for_symbol.get(gene_symbol):
                            if len(potential_ga):
                                fixed_ga = next(iter(potential_ga))
                    if fixed_ga:
                        gene_accession = fixed_ga
                historical_transcript_version["gene_version"] = gene_accession
                transcript_versions[transcript_accession] = historical_transcript_version

    genes = {}  # Only keep those that are used in transcript versions
    # Summarise where it's from
    transcript_urls = Counter()
    for tv in transcript_versions.values():
        if not args.no_genes:
            if gene_accession := tv.get("gene_version"):
                genes[gene_accession] = gene_versions[gene_accession]

        for build_coordinates in tv["genome_builds"].values():
            transcript_urls[build_coordinates["url"]] += 1

    total = sum(transcript_urls.values())
    print(f"{total} transcript versions from:")
    for url, count in transcript_urls.most_common():
        print(f"{url}: {count} ({count*100 / total:.1f}%)")

    write_cdot_json(args.output, genes, transcript_versions, [args.genome_build])


def combine_builds(args):
    print("combine_builds")
    genome_build_file = {
        "GRCh37": gzip.open(args.grch37),
        "GRCh38": gzip.open(args.grch38),
    }

    urls_different_coding = defaultdict(list)
    transcripts = {}
    for genome_build, f in genome_build_file.items():
        # TODO: Check cdot versions
        json_builds = next(ijson.items(f, "genome_builds"))
        if json_builds != [genome_build]:
            raise ValueError(f"JSON file provided for {genome_build} needs to have only {genome_build} data (has {json_builds})")

        f.seek(0)  # Reset for next ijson call
        for transcript_id, build_transcript in ijson.kvitems(f, "transcripts"):
            existing_transcript = transcripts.get(transcript_id)
            genome_builds = {}
            if existing_transcript:
                genome_builds = existing_transcript["genome_builds"]
                # Latest always used, but check existing - if codons are different old versions are wrong so remove
                for field in ["start_codon", "stop_codon"]:
                    old = existing_transcript.get(field)
                    new = build_transcript.get(field)
                    if old != new:  # Old relied on different codons so is obsolete
                        for build_coordinates in genome_builds.values():
                            url = build_coordinates["url"]
                            urls_different_coding[url].append(transcript_id)
                        genome_builds = {}

            genome_builds[genome_build] = build_transcript["genome_builds"][genome_build]
            # Use latest (with merged genome builds)
            build_transcript["genome_builds"] = genome_builds
            transcripts[transcript_id] = build_transcript
        f.close()

    print("Writing cdot data")
    with gzip.open(args.output, 'wt') as outfile:
        data = {
            "transcripts": transcripts,
            "cdot_version": cdot.__version__,
            "genome_builds": list(genome_build_file.keys()),
        }
        json.dump(data, outfile)

    if urls_different_coding:
        print("Some transcripts were removed as they had different coding coordinates from latest")
        for url, transcript_ids in urls_different_coding.items():
            print(f"{url}: {','.join(transcript_ids)}")


def main():
    args = handle_args()
    if args.version:
        print(cdot.__version__)
        sys.exit(0)

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
