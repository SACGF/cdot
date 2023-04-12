#!/bin/bash

set -e

BASE_DIR=$(dirname ${BASH_SOURCE[0]})
CDOT_VERSION=$(${BASE_DIR}/cdot_json.py --version)

if [[ -z ${GENE_INFO_JSON} ]]; then
  echo "You need to set environment variable GENE_INFO_JSON, pointing to the filename produced by cdot_gene_info.py"
  exit 1
fi

merge_args=()

filename=GCF_009914755.1_T2T-CHM13v2.0_genomic.110.gff.gz
url=https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/110/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
cdot_file=cdot-${CDOT_VERSION}.$(basename $filename .gz).json.gz

if [[ ! -e ${filename} ]]; then
  wget ${url} --output-document=${filename}
fi
if [[ ! -e ${cdot_file} ]]; then
  ${BASE_DIR}/cdot_json.py gff3_to_json "${filename}" --url "${url}" --genome-build=CHM13v2.0 --output "${cdot_file}" --gene-info-json="${GENE_INFO_JSON}"
fi
merge_args+=(${cdot_file})


filename=GCF_009914755.1_T2T-CHM13v2.0_genomic.RS_2023_03.gff.gz
url=https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/GCF_009914755.1-RS_2023_03/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
cdot_file=cdot-${CDOT_VERSION}.$(basename $filename .gz).json.gz

if [[ ! -e ${filename} ]]; then
  wget ${url} --output-document=${filename}
fi
if [[ ! -e ${cdot_file} ]]; then
  ${BASE_DIR}/cdot_json.py gff3_to_json "${filename}" --url "${url}" --genome-build=CHM13v2.0 --output "${cdot_file}" --gene-info-json="${GENE_INFO_JSON}"
fi
merge_args+=(${cdot_file})

merged_file="cdot-${CDOT_VERSION}.refseq.CHM13v2.0.json.gz"
if [[ ! -e ${merged_file} ]]; then
  echo "Creating ${merged_file}"
  ${BASE_DIR}/cdot_json.py merge_historical ${merge_args[@]} --genome-build=CHM13v2.0 --output "${merged_file}"
fi