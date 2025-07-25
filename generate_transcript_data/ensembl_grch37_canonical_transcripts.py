#!/bin/env python3

# We want to be able to flag GRCh37 Ensembl transcripts as canonical
# @see https://github.com/SACGF/cdot/issues/36
# This is not in GTF, in Gencode metadata
# Also not in MySQL database (transcript column is_canonical or canonical attributes not present in GRCh37)
# GRCh37 is not being updated so we can just do this once and store it.
import csv
import gzip
import json
import logging
from argparse import ArgumentParser
# The best way I have found to get it is Ensembl REST
# https://grch37.rest.ensembl.org/lookup/ENSG00000159216?content-type=application/json
#
# {
#   "display_name": "RUNX1",
#   "strand": -1,
#   "canonical_transcript": "ENST00000300305.3",
#   "source": "ensembl_havana",
#   "start": 36160098,
#   "object_type": "Gene",
#   "end": 37376965,
#   "biotype": "protein_coding",
#   "assembly_name": "GRCh37",
#   "db_type": "core",
#   "id": "ENSG00000159216",
#   "seq_region_name": "21",
#   "species": "homo_sapiens",
#   "version": 14,
#   "logic_name": "ensembl_havana_gene_homo_sapiens_37",
#   "description": "runt-related transcription factor 1 [Source:HGNC Symbol;Acc:10471]"
# }


from typing import Iterable, TypeVar, Iterator

T = TypeVar("T")

def batch_iterator(iterable: Iterable[T], batch_size: int = 100) -> Iterator[list[T]]:
    batch: list[T] = []
    for record in iterable:
        batch.append(record)
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch

import requests

def fetch_canonical_transcripts(gene_list: list[str]) -> dict[str, str]:
    url = "https://grch37.rest.ensembl.org/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    payload = {"ids": gene_list}

    response = requests.post(url, headers=headers, json=payload)
    response.raise_for_status()
    data = response.json()

    canonical_transcripts = {}
    for gene_id, v in data.items():
        if t := v.get("canonical_transcript"):
            canonical_transcripts[gene_id] = t
    return canonical_transcripts


def handle_args():
    parser = ArgumentParser(description='cdot Ensembl GRCh37 Gene Accession Canonical Transcript retrieval')
    parser.add_argument('--cdot-json', required=True, help='cdot json file')
    parser.add_argument('--output-csv', required=True, help='output filename')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    args = handle_args()
    data = json.load(gzip.open(args.cdot_json))
    genes = list(data["genes"])
    del data
    logging.info("Retrieving canonical transcripts for %d genes", len(genes))

    num_gene_requests = 0
    num_gene_transcripts = 0
    last_update = 0
    update_numbers = 5000

    with open(args.output_csv, 'wt') as f:
        writer = csv.writer(f)
        writer.writerow(["gene_id", "canonical_transcript"])

        for gene_list in batch_iterator(genes):
            num_gene_requests += len(gene_list)
            canonical_map = fetch_canonical_transcripts(gene_list)
            for gene_id, transcript_id in canonical_map.items():
                num_gene_transcripts += 1
                writer.writerow([gene_id, transcript_id])

            if num_gene_requests - last_update >= update_numbers:
                last_update = num_gene_requests
                logging.info(f"Gene requests: {num_gene_requests}, transcripts: {num_gene_transcripts}")
