#!/bin/bash

set -e

BASE_DIR=$(dirname $(dirname ${BASH_SOURCE[0]}))
CDOT_VERSION=$(${BASE_DIR}/cdot_json.py --version)

if [[ -z ${GENE_INFO_JSON} ]]; then
  echo "You need to set environment variable GENE_INFO_JSON, pointing to the filename produced by cdot_gene_info.py"
  exit 1
fi


filename=GCF_000001635.27_GRCm39_genomic.gff.gz
url=https://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/annotation_releases/109/GCF_000001635.27_GRCm39/${filename}
cdot_file=cdot-${CDOT_VERSION}.$(basename $filename .gz).json.gz
if [[ ! -e ${filename} ]]; then
  wget ${url}
fi
if [[ ! -e ${cdot_file} ]]; then
  ${BASE_DIR}/cdot_json.py gff3_to_json "${filename}" --no-contig-conversion --url "${url}" --genome-build=GRCm39 --output "${cdot_file}" --gene-info-json="${GENE_INFO_JSON}"
fi
