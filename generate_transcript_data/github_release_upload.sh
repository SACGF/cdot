#!/bin/bash

set -e

if [[ -z ${CDOT_DATA_DIR} ]]; then
  echo "You need to set environment variable CDOT_DATA_DIR, pointing to where you ran 'all_transcripts.sh'"
  exit 1
fi

FULL_PATH_TO_SCRIPT="$(realpath "${BASH_SOURCE[-1]}")"
BASE_DIR=$(dirname ${FULL_PATH_TO_SCRIPT})

# Python scripts will import via generate_transcript_data
export PYTHONPATH=${BASE_DIR}/..

CDOT_DATA_VERSION=$(${BASE_DIR}/cdot_json.py --version)

CDOT_RELEASE_NAME=data_v${CDOT_DATA_VERSION}
echo "For the rest of the script to work, it assumes you have tagged + pushed a data release of ${CDOT_DATA_VERSION}"
echo "then run: gh release create ${CDOT_RELEASE_NAME} --title=${CDOT_RELEASE_NAME} --notes 'release notes...'"

gh release upload ${CDOT_RELEASE_NAME} \
  ${CDOT_DATA_DIR}/ensembl/GRCh37/cdot-${CDOT_DATA_VERSION}.ensembl.grch37.json.gz \
  ${CDOT_DATA_DIR}/ensembl/GRCh37/cdot-${CDOT_DATA_VERSION}.ensembl.Homo_sapiens.GRCh37.87.gff3.json.gz \
  ${CDOT_DATA_DIR}/ensembl/GRCh38/cdot-${CDOT_DATA_VERSION}.ensembl.grch38.json.gz \
  ${CDOT_DATA_DIR}/ensembl/GRCh38/cdot-${CDOT_DATA_VERSION}.ensembl.Homo_sapiens.GRCh38.110.gff3.json.gz \
  ${CDOT_DATA_DIR}/ensembl/GRCh38/cdot-${CDOT_DATA_VERSION}.ensembl.Homo_sapiens.GRCh38.111.gff3.json.gz \
  ${CDOT_DATA_DIR}/ensembl/T2T-CHM13v2.0/cdot-${CDOT_DATA_VERSION}.ensembl.T2T-CHM13v2.0.json.gz \
  ${CDOT_DATA_DIR}/refseq/GRCh37/cdot-${CDOT_DATA_VERSION}.refseq.grch37.json.gz \
  ${CDOT_DATA_DIR}/refseq/GRCh37/cdot-${CDOT_DATA_VERSION}.GCF_000001405.25_GRCh37.p13_genomic.105.20201022.gff.json.gz \
  ${CDOT_DATA_DIR}/refseq/GRCh37/cdot-${CDOT_DATA_VERSION}.GCF_000001405.25_GRCh37.p13_genomic.105.20220307.gff.json.gz \
  ${CDOT_DATA_DIR}/refseq/GRCh38/cdot-${CDOT_DATA_VERSION}.refseq.grch38.json.gz \
  ${CDOT_DATA_DIR}/refseq/GRCh38/cdot-${CDOT_DATA_VERSION}.GCF_000001405.40_GRCh38.p14_genomic.110.gff.json.gz \
  ${CDOT_DATA_DIR}/refseq/GRCh38/cdot-${CDOT_DATA_VERSION}.GCF_000001405.40_GRCh38.p14_genomic.RS_2023_10.gff.json.gz \
  ${CDOT_DATA_DIR}/refseq/T2T-CHM13v2.0/cdot-${CDOT_DATA_VERSION}.refseq.T2T-CHM13v2.0.json.gz

