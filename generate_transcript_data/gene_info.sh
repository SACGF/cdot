#!/bin/bash

if [[ -z ${EMAIL} ]]; then
  echo "You need to set the 'EMAIL' shell variable (used for NCBI API calls)"
  exit
fi

BASE_DIR=$(dirname ${BASH_SOURCE[0]})
CDOT_VERSION=$(${BASE_DIR}/cdot_json.py --version)
TODAY=$(date --iso)

filename=Homo_sapiens.gene_info.gz
url=https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/${filename}
if [[ ! -e ${filename} ]]; then
  wget ${url}
fi

out_json=gene_summary-${CDOT_VERSION}-${TODAY}.json.gz
${BASE_DIR}/cdot_gene_info.py --gene-info ${filename} --output ${out_json} --email ${EMAIL}