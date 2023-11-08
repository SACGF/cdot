#!/usr/bin/env python3

import csv
import gzip
import json
import os
from argparse import ArgumentParser
from csv import DictReader
from datetime import datetime
from typing import Iterable, Iterator, List, TypeVar

import cdot
from Bio import Entrez
from cdot.json_encoders import SortedSetEncoder
from io import BytesIO
from lxml import etree

T = TypeVar("T")


def handle_args():
    parser = ArgumentParser(description='cdot Gene Info retrieval')
    parser.add_argument('--email', required=True, help='Entrez email')
    parser.add_argument('--gene-info', required=True, help='refseq gene info file')
    parser.add_argument('--output', required=True, help='output filename')

    args = parser.parse_args()
    return args


def batch_iterator(iterable: Iterable[T], batch_size: int = 10) -> Iterator[List[T]]:
    batch: List[T] = list()
    for record in iterable:
        batch.append(record)
        if len(batch) >= batch_size:
            yield batch
            batch = list()
    if batch:
        yield batch


def _get_entrez_gene_summary(id_list):
    request = Entrez.epost("gene", id=",".join(id_list))
    result = Entrez.read(request)
    web_env = result["WebEnv"]
    query_key = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=web_env, query_key=query_key)
    document = Entrez.read(data, ignore_errors=True)  # Need recent BioPython
    return document["DocumentSummarySet"]["DocumentSummary"]


def iter_entrez_ids(reader):
    for gi in reader:
        if gi["Symbol_from_nomenclature_authority"] != '-':
            yield gi['GeneID']

def main():
    args = handle_args()
    Entrez.email = args.email  # Stop warning message
    start_date = datetime.now().isoformat()

    # 10k limit of return data from NCBI
    NCBI_BATCH_SIZE = 10000

    gene_info = {}
    with gzip.open(args.gene_info, "rt") as f:
        reader = csv.DictReader(f, dialect='excel-tab')

        for entrez_ids in batch_iterator(iter_entrez_ids(reader), batch_size=NCBI_BATCH_SIZE):
            # We should really store it under the gene Id so dupe symbols don't wipe
            for gene_summary in _get_entrez_gene_summary(entrez_ids):
                gene_id = gene_summary.attributes["uid"]
                gene_info[gene_id] = {
                    "gene_symbol": gene_summary["NomenclatureSymbol"],
                    "map_location": gene_summary["MapLocation"],
                    # Already have description for RefSeq but not Ensembl (will just overwrite)
                    "description": gene_summary["NomenclatureName"],
                    # "added": record["date_name_changed"],
                    "aliases": gene_summary["OtherAliases"],
                    "summary": gene_summary["Summary"],
                }

            print(f"Processed {len(gene_info)} records")

    if gene_info:
        with gzip.open(args.output, 'wt') as outfile:
            gene_info_file_dt = datetime.fromtimestamp(os.stat(args.gene_info).st_ctime)

            data = {
                "cdot_version": cdot.__version__,
                "api_retrieval_date": start_date,
                "gene_info_date": gene_info_file_dt.isoformat(),
                "gene_info": gene_info,
            }
            json.dump(data, outfile, cls=SortedSetEncoder, sort_keys=True)  # Sort so diffs work


if __name__ == '__main__':
    main()
