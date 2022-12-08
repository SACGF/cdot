#!/bin/bash

set -e

BASE_DIR=$(dirname ${BASH_SOURCE[0]})
CDOT_VERSION=$(${BASE_DIR}/cdot_json.py --version)

if [[ -z ${GENE_INFO_JSON} ]]; then
  echo "You need to set environment variable GENE_INFO_JSON, pointing to the filename produced by cdot_gene_info.py"
  exit 1
fi

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
for release in 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108; do
  # Switched to using GTFs as they contain protein version
  filename=Homo_sapiens.GRCh38.${release}.gff3.gz
  url=ftp://ftp.ensembl.org/pub/release-${release}/gff3/homo_sapiens/${filename}
  cdot_file=cdot-${CDOT_VERSION}.$(basename $filename .gz).json.gz

  if [[ ! -e ${filename} ]]; then
    wget ${url}
  fi
  if [[ ! -e ${cdot_file} ]]; then
    ${BASE_DIR}/cdot_json.py gff3_to_json "${filename}" --url "${url}" --genome-build=GRCh38 --output "${cdot_file}" --gene-info-json="${GENE_INFO_JSON}"
  fi
  merge_args+=(${cdot_file})
done

merged_file="cdot-${CDOT_VERSION}.ensembl.grch38.json.gz"
if [[ ! -e ${merged_file} ]]; then
  ${BASE_DIR}/cdot_json.py merge_historical ${merge_args[@]} --genome-build=GRCh38 --output "${merged_file}"
fi
