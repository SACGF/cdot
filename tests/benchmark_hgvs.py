#!/bin/env python3

"""
    See instructions at end of file on how to extract test HGVS from clinvar
"""

import time
import pandas as pd
from argparse import ArgumentParser

import hgvs
import hgvs.dataproviders.uta
from hgvs.assemblymapper import AssemblyMapper
from hgvs.exceptions import HGVSDataNotAvailableError, HGVSInvalidVariantError

from cdot.hgvs.dataproviders import JSONDataProvider, RESTDataProvider


def handle_args():
    parser = ArgumentParser(description='Benchmark cdot')
    parser.add_argument("hgvs_file")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--uta', action='store_true')
    group.add_argument('--rest', action='store_true')
    group.add_argument('--rest-insecure', action='store_true')
    parser.add_argument('--json', help='JSON file')
    args = parser.parse_args()
    if not any([args.uta, args.rest, args.rest_insecure, args.json]):
        parser.error("You need to specify at least one of 'uta', 'rest', 'rest-insecure', 'json'")
    return args


def main():
    args = handle_args()

    hgvs_g_c_list = []
    with open(args.hgvs_file) as f:
        for line in f:
            hgvs_g_c_list.append(line.split())

    total = len(hgvs_g_c_list)
    print(f"Using {total} test records")

    if args.uta:
        hdp = hgvs.dataproviders.uta.connect()
    elif args.rest:
        hdp = RESTDataProvider()  # Uses API server at cdot.cc
    elif args.json:
        hdp = JSONDataProvider([args.json])
    elif args.json_insecure:
        hdp = JSONDataProvider([args.json], secure=False)
    else:
        raise ValueError("Unknown data provider method!")

    am = AssemblyMapper(hdp,
                        assembly_name='GRCh38',
                        alt_aln_method='splign', replace_reference=True)

    hp = hgvs.parser.Parser()

    run_times = []
    correct = 0
    incorrect = 0
    no_data = 0
    for hgvs_g, hgvs_c in hgvs_g_c_list:
        start = time.time()
        try:
            var_c = hp.parse_hgvs_variant(hgvs_c)
            if ":c." in hgvs_c:
                converted_hgvs_g = str(am.c_to_g(var_c))
            else:
                converted_hgvs_g = str(am.n_to_g(var_c))
        except HGVSDataNotAvailableError:
            no_data += 1
            continue
        except HGVSInvalidVariantError as ive:
            print(f"{hgvs_c}: {ive}")
            incorrect += 1
            continue

        if converted_hgvs_g == hgvs_g:
            correct += 1
        else:
            incorrect += 1
            print(f"{hgvs_c}: '{hgvs_g}' != '{converted_hgvs_g}' (actual)")
            continue

        # We only keep times for correct data
        end = time.time()
        time_taken = end - start
        run_times.append(time_taken)

    print(f"Total: {total}, correct: {correct}, incorrect: {incorrect}, no data: {no_data}")
    df = pd.DataFrame(run_times)
    print(df.describe())


if __name__ == '__main__':
    main()

"""

* Get a subset of rows from ClinVar VCF
    * zgrep "^#" clinvar.vcf.gz > header.txt
    * zcat -v "^#" clinvar.vcf.gz | shuf -n 1000 > clinvar_1k_records.vcf
    * cat header.txt clinvar_1k_rows.vcf | gzip > clinvar_1k.vcf.gz

* Annotate the VCF to get MANE transcript (via --pick)

    vep -i clinvar_1k.vcf.gz -o clinvar_1k.vep_annotated.vcf.gz --cache --dir /data/annotation/VEP/vep_cache --fasta /data/annotation/fasta/GCF_000001405.39_GRCh38.p13_genomic.fna.gz --assembly GRCh38 --offline --use_given_ref --vcf --compress_output gzip --force_overwrite --pick --no_escape --hgvs --refseq --buffer_size 1000

* Extract out the g.HGVS and c.HGVS




def cyvcf2_header_types(cyvcf2_reader):
    header_types = defaultdict(dict)
    for h in cyvcf2_reader.header_iter():
        info = h.info()
        h_id = info.get("ID")
        if h_id:  # Not much use w/o this
            header_types[h.type][h_id] = info
    return header_types


reader = Reader("./clinvar_1k.vcf.gz")
header_types = cyvcf2_header_types(reader)
description = header_types["INFO"]["CSQ"]["description"]
description = description.replace('"', '')  # Strip double quotes

match = "Format: "
columns_str = description[description.rfind(match) + len(match):]
vep_columns = columns_str.split("|")

hgvs = []
for v in reader:
    csq = v.INFO.get("CSQ")
    td = dict(zip(vep_columns, csq.split("|")))
    g_hgvs = v.INFO.get("CLNHGVS")
    c_hgvs = td.get("HGVSc")
    if g_hgvs and c_hgvs:
        hgvs.append((g_hgvs, c_hgvs))


"""
