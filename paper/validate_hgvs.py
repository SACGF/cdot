#!/bin/env python3

import math
import re
import sys
import pandas as pd
from pysam.libcfaidx import FastaFile
import pyhgvs
from cdot.pyhgvs.pyhgvs_transcript import JSONPyHGVSTranscriptFactory


def get_combined_df():
    servers = ["shariant", "vg-aws", "vg3_upgrade"]

    df_list = []

    for server in servers:
        filename = f"{server}_search_hgvs.csv"
        server_df = pd.read_csv(filename, names=["date", "details"], skiprows=1)
        server_df["server"] = server
        df_list.append(server_df)

    return pd.concat(df_list) 


def add_hgvs_column(df):
    pattern_calculated = re.compile(r"'(.+)' calculated")
    pattern_type = re.compile(r"'(.+)' = type")
    pattern_returned = re.compile(r"'(.+)' returned")

    hgvs_list = []
    for details in df["details"].values:
        for pattern in [pattern_calculated, pattern_type, pattern_returned]:
            if m := pattern.match(details):
                hgvs_list.append(m.group(1))
                break
        else:
            print(f"No match for '{details}'")
            hgvs_list.append("")   
    
    df["hgvs"] = hgvs_list


def can_resolve(genome, factory, hgvs_c):
    try:
        pyhgvs.parse_hgvs_name(hgvs_c, genome, get_transcript=factory.get_transcript_grch37) # 37
        return True
    except Exception as e:
        print(e)
        try:
            pyhgvs.parse_hgvs_name(hgvs_c, genome, get_transcript=factory.get_transcript_grch38) # 38
            return True
        except Exception as e2:
            print(e2)
            pass
    return False


def add_hgvs_validation_columns(df):
    genome = FastaFile("/data/annotation/fasta/GCF_000001405.25_GRCh37.p13_genomic.fna.gz")
    factory = JSONPyHGVSTranscriptFactory(["/home/dlawrence/Downloads/cdot-0.2.12.refseq.grch37_grch38.json.gz",
                                           "/home/dlawrence/Downloads/cdot-0.2.12.ensembl.grch37_grch38.json.gz"])

    valid_hgvs_list = []
    can_resolve_list = []
    for hgvs_c in df["hgvs"].values:
        # print(f"testing... {hgvs_c}")
        resolve_ok = False
        try:
            pyhgvs.HGVSName(hgvs_c)
            valid_hgvs = True
            resolve_ok = can_resolve(genome, factory, hgvs_c)
        except:
            valid_hgvs = False

        valid_hgvs_list.append(valid_hgvs)
        can_resolve_list.append(resolve_ok)
    
    df["valid_hgvs"] = valid_hgvs_list
    df["can_resolve"] = can_resolve_list


def split_df_chunks(data_df,chunk_size):
    """ From https://xhinker.medium.com/python-split-a-dataframe-to-a-chunk-list-fe80bf9d63be  """
    total_length     = len(data_df)
    total_chunk_num  = math.ceil(total_length/chunk_size)
    normal_chunk_num = math.floor(total_length/chunk_size)
    chunks = []
    for i in range(normal_chunk_num):
        chunk = data_df[(i*chunk_size):((i+1)*chunk_size)]
        chunks.append(chunk)
    if total_chunk_num > normal_chunk_num:
        chunk = data_df[(normal_chunk_num*chunk_size):total_length]
        chunks.append(chunk)
    return chunks


def main():
    if len(sys.argv) == 1:
        print("main")
        df = get_combined_df()
        add_hgvs_column(df)
        for i, chunk in enumerate(split_df_chunks(df, 500)):
            filename = f"hgvs_search_{i}.csv"
            print(f"writing {filename}")
            chunk.to_csv(filename)
    else:
        filename = sys.argv[1]
        print(f"Processing {filename}")
        df = pd.read_csv(filename)
        add_hgvs_validation_columns(df)
        df.to_csv(f"validate_{filename}") 


if __name__ == "__main__":
    main() 
