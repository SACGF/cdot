#!/bin/bash

set -e

BASE_DIR=$(dirname ${BASH_SOURCE[0]})
CDOT_VERSION=$(${BASE_DIR}/cdot_json.py --version)

if [[ -z ${GENE_INFO_JSON} ]]; then
  echo "You need to set environment variable GENE_INFO_JSON, pointing to the filename produced by cdot_gene_info.py"
  exit 1
fi

merge_args=()

filename=GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
url=https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/GCF_009914755.1-RS_2023_03/${filename}
cdot_file=cdot-${CDOT_VERSION}.$(basename $filename .gz).json.gz

if [[ ! -e ${filename} ]]; then
  wget ${url}
fi
if [[ ! -e ${cdot_file} ]]; then
  ${BASE_DIR}/cdot_json.py gff3_to_json "${filename}" --url "${url}" --genome-build=CHM13v2.0 --output "${cdot_file}" --gene-info-json="${GENE_INFO_JSON}"
fi
merge_args+=(${cdot_file})

# If we add a second one - need a 'merge_historical', see other scripts in this dir