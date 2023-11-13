#!/bin/bash

set -e

BASE_DIR=$(dirname ${BASH_SOURCE[0]})
CDOT_VERSION=$(${BASE_DIR}/cdot_json.py --version)
GENOME_BUILD=T2T-CHM13v2.0

if [[ -z ${GENE_INFO_JSON} ]]; then
  echo "You need to set environment variable GENE_INFO_JSON, pointing to the filename produced by cdot_gene_info.py"
  exit 1
fi

merge_args=()
for release in 2022_06 2022_07; do
  filename=Homo_sapiens-GCA_009914755.4-${release}-genes.gff3.gz
  url=https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/${release}/${filename}
  cdot_file=cdot-${CDOT_VERSION}.ensembl.$(basename $filename .gz).json.gz

  if [[ ! -e ${filename} ]]; then
    wget ${url}
  fi
  if [[ ! -e ${cdot_file} ]]; then
    ${BASE_DIR}/cdot_json.py gff3_to_json "${filename}" --url "${url}" --genome-build=${GENOME_BUILD} --output "${cdot_file}" --gene-info-json="${GENE_INFO_JSON}"
  fi
  merge_args+=(${cdot_file})
done

merged_file="cdot-${CDOT_VERSION}.ensembl.${GENOME_BUILD}.json.gz"
if [[ ! -e ${merged_file} ]]; then
  ${BASE_DIR}/cdot_json.py merge_historical ${merge_args[@]} --genome-build=${GENOME_BUILD} --output "${merged_file}"
fi
