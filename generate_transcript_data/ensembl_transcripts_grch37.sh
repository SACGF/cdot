#!/bin/bash

set -e

BASE_DIR=$(dirname ${BASH_SOURCE[0]})
CDOT_VERSION=$(${BASE_DIR}/cdot_json.py --version)

if [[ -z ${GENE_SUMMARY} ]]; then
  echo "You need to set environment variable GENE_SUMMARY, pointing to the filename produced by cdot_gene_info.py"
  exit 1
fi

# v81 (points to 75) and earlier at GTFs that don't have transcript versions - just skip them

#82 is first GFF3 for GRCh37
#83 has no data
#84 is 82 again
#86 is 85 again
merge_args=()
for release in 82 85 87; do
  # Switched to using GTFs as they contain protein version
  filename=Homo_sapiens.GRCh37.${release}.gtf.gz
  url=ftp://ftp.ensembl.org/pub/grch37/release-${release}/gtf/homo_sapiens/${filename}
  cdot_file=$(basename $filename .gz).json.gz
  if [[ ! -e ${filename} ]]; then
    wget ${url}
  fi
  if [[ ! -e ${cdot_file} ]]; then
    ${BASE_DIR}/cdot_json.py gtf_to_json "${filename}" --url "${url}" --genome-build=GRCh37 --output "${cdot_file}" --gene-info-json="${GENE_SUMMARY}"
  fi
  merge_args+=(${cdot_file})
done

merged_file="cdot-${CDOT_VERSION}.ensembl.grch37.json.gz"
if [[ ! -e ${merged_file} ]]; then
  ${BASE_DIR}/cdot_json.py merge_historical ${merge_args[@]} --genome-build=GRCh37 --output "${merged_file}"
fi
