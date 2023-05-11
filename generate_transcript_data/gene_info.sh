#!/bin/bash

if [[ -z ${EMAIL} ]]; then
  echo "You need to set the 'EMAIL' shell variable (used for NCBI API calls)"
  exit
fi

BASE_DIR=$(dirname ${BASH_SOURCE[0]})
CDOT_VERSION=$(${BASE_DIR}/cdot_json.py --version)
REFSEQ_DIR=${REFSEQ_DIR:-H_sapiens}
SPECIES=${SPECIES:-Homo_sapiens}

echo "Generating Gene Info for REFSEQ_DIR=${REFSEQ_DIR}, SPECIES=${SPECIES}"

filename=${SPECIES}.gene_info.gz
url=https://ftp.ncbi.nlm.nih.gov/refseq/${REFSEQ_DIR}/${filename}
if [[ ! -e ${filename} ]]; then
  echo "Downloading ${url}"
  wget ${url}
fi

out_json=${SPECIES}.gene-info-${CDOT_VERSION}.json.gz
if [[ ! -e ${out_json} ]]; then
  echo "Processing gene info file..."
  ${BASE_DIR}/cdot_gene_info.py --gene-info ${filename} --output ${out_json} --email ${EMAIL}
fi
