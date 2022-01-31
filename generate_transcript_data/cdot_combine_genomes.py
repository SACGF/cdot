#!/bin/env python3

"""
    Stores transcript version data, with a "genome_builds" dict with data for GRCh37/GRCh38

    We could have stored this per contig, which would save a bit of space in the case where contigs are shared
    across genomes (eg MT) and this may have made more sense for BioCommons HGVS - but I think it would be less
    obvious for other users of the data

    I'm also hardcoding the genome builds as they are not likely to change much (HGVS is only for human)

"""

import gzip
import json

from argparse import ArgumentParser
from collections import defaultdict

import cdot


def handle_args():
    parser = ArgumentParser(description='Convert multiple cdot json.gz files into one')
    parser.add_argument('--grch37', required=True, help='cdot JSON.gz for GRCh37')
    parser.add_argument('--grch38', required=True, help='cdot JSON.gz for GRCh38')
    parser.add_argument('--output', required=True, help='Output filename')
    return parser.parse_args()


def main():
    args = handle_args()

    genome_build_data = {
        "GRCh37": json.load(gzip.open(args.grch37)),
        "GRCh38": json.load(gzip.open(args.grch38)),
    }

    all_transcript_ids = set()
    for genome_build, data in genome_build_data.items():
        # TODO: Check cdot versions
        json_builds  = data["genome_builds"]
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


if __name__ == '__main__':
    main()
