"""
    Converts PyReference JSON.gz files (created from RefSeq/Ensembl GTF or GFFs) into format easy to load for HGVS

    @see https://bitbucket.org/sacgf/pyreference/

    Dave Lawrence (davmlaw@gmail.com) on 17/01/2022
"""
import gzip
import json
from argparse import ArgumentParser
from collections import Counter
from typing import Dict

import cdot


def handle_args():
    parser = ArgumentParser(description='Convert multiple PyReference json.gz files into one for cdot')
    parser.add_argument('--pyreference-json', required=True, nargs="+", action="extend",
                        help='PyReference JSON.gz - list OLDEST to NEWEST (newest is kept)')
    parser.add_argument("--store-genes", action='store_true', help="Also store gene/version information")
    parser.add_argument('--output', required=True, help='Output filename')
    return parser.parse_args()


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


def main():
    args = handle_args()

    gene_versions = {}  # We only keep those that are in the latest transcript version
    transcript_versions = {}

    for pyref_filename in args.pyreference_json:
        print(f"Loading '{pyref_filename}'")
        with gzip.open(pyref_filename) as f:
            pyref_data = json.load(f)

            url = pyref_data["reference_gtf"]["url"]

            # PyReference stores transcripts under genes, while PyReference only has transcripts (that contain genes)
            transcript_gene_version = {}

            for gene_id, gene in pyref_data["genes_by_id"].items():
                if version := gene.get("version"):
                    gene_accession = f"{gene_id}.{version}"
                else:
                    gene_accession = gene_id

                gene_version = convert_gene_pyreference_to_gene_version_data(gene)
                gene_version["url"] = url
                gene_versions[gene_accession] = gene_version

                for transcript_accession in gene["transcripts"]:
                    transcript_gene_version[transcript_accession] = gene_accession

            for transcript_accession, transcript_version in pyref_data["transcripts_by_id"].items():
                gene_accession = transcript_gene_version[transcript_accession]
                gene_version = gene_versions[gene_accession]
                transcript_version["id"] = transcript_accession
                if not gene_accession.startswith("_"):
                    transcript_version["gene_version"] = gene_accession
                transcript_version["gene_name"] = gene_version["gene_symbol"]
                transcript_version["url"] = url
                transcript_versions[transcript_accession] = transcript_version
                if hgnc := gene_version.get("hgnc"):
                    transcript_version["hgnc"] = hgnc

    # Summarise where it's from
    transcript_urls = Counter()
    for tv in transcript_versions.values():
        transcript_urls[tv["url"]] += 1

    total = sum(transcript_urls.values())
    print(f"{total} transcript versions from:")
    for url, count in transcript_urls.most_common():
        print(f"{url}: {count} ({count*100 / total:.1f}%)")

    print("Writing cdot data")
    with gzip.open(args.output, 'w') as outfile:
        data = {
            "transcripts": transcript_versions,
            "cdot_version": cdot.__version__,
        }
        if args.store_genes:
            data["genes"] = gene_versions

        json_str = json.dumps(data)
        outfile.write(json_str.encode('ascii'))


if __name__ == '__main__':
    main()
