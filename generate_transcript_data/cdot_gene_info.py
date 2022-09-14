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
    document = Entrez.read(data)
    return document["DocumentSummarySet"]["DocumentSummary"]


def main():
    args = handle_args()
    Entrez.email = args.email  # Stop warning message

    # 10k limit of return data from NCBI
    NCBI_BATCH_SIZE = 2000

    gene_info = {}
    with gzip.open(args.gene_info, "rt") as f:
        reader = csv.DictReader(f, dialect='excel-tab')

        for batch in batch_iterator(reader, batch_size=NCBI_BATCH_SIZE):
            entrez_ids = []
            for gi in batch:
                if gi["Symbol_from_nomenclature_authority"] != '-':
                    entrez_ids.append(gi['GeneID'])

            if entrez_ids:
                for gene_summary in _get_entrez_gene_summary(entrez_ids):
                    gene_symbol = gene_summary["NomenclatureSymbol"]
                    gene_info[gene_symbol] = {
                        "map_location": gene_summary["MapLocation"],
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
                "date": datetime.now().isoformat(),
                "gene_info_date": gene_info_file_dt.isoformat(),
                "gene_info": gene_info,
            }
            json.dump(data, outfile, cls=SortedSetEncoder, sort_keys=True)  # Sort so diffs work


if __name__ == '__main__':
    main()
