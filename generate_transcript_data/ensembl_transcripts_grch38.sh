#!/bin/bash

set -e

BASE_DIR=$(dirname ${BASH_SOURCE[0]})
CDOT_VERSION=$(${BASE_DIR}/cdot_json.py --version)

# Skip earlier GTFs as they don't have versions
#for release in 76 77 78 79 80; do
#  filename=Homo_sapiens.GRCh38.${release}.gtf.gz
#  url=ftp://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens/${filename}
#  if [[ ! -e ${filename} ]]; then
#    wget ${url}
#  fi

#  if [[ ! -e ${filename}.json.gz ]]; then
#    pyreference_gff_to_json.py --url "${url}" --gff3 "${filename}"
#  fi
#done

#81 is first GFF3 for GRCh38
merge_args=()
for release in 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106; do
  # Switched to using GTFs as they contain protein version
  filename=Homo_sapiens.GRCh38.${release}.gtf.gz
  url=ftp://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens/${filename}
  cdot_file=$(basename $filename .gz).json.gz

  if [[ ! -e ${filename} ]]; then
    wget ${url}
  fi
  if [[ ! -e ${cdot_file} ]]; then
    ${BASE_DIR}/cdot_json.py gtf_to_json "${filename}" --url "${url}" --genome-build=GRCh38 --output "${cdot_file}"
  fi
  merge_args+=(${cdot_file})
done

merged_file="cdot-${CDOT_VERSION}.ensembl.grch38.json.gz"
if [[ ! -e ${merged_file} ]]; then
  ${BASE_DIR}/cdot_json.py merge_historical ${merge_args[@]} --genome-build=GRCh38 --output "${merged_file}"
fi
