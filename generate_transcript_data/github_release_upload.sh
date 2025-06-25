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
  ${CDOT_DATA_DIR}/Ensembl/cdot-${CDOT_DATA_VERSION}-all-builds-Ensembl-grch37_grch38_t2t-chm13v2.0.json.gz \
  ${CDOT_DATA_DIR}/Ensembl/cdot-${CDOT_DATA_VERSION}-Ensembl-merged-GRCh37.json.gz \
  ${CDOT_DATA_DIR}/Ensembl/GRCh37/cdot-${CDOT_DATA_VERSION}-Homo_sapiens_GRCh37_Ensembl_87.gtf.json.gz \
  ${CDOT_DATA_DIR}/Ensembl/cdot-${CDOT_DATA_VERSION}-Ensembl-merged-GRCh38.json.gz \
  ${CDOT_DATA_DIR}/Ensembl/GRCh38/cdot-${CDOT_DATA_VERSION}-Homo_sapiens_GRCh38_Ensembl_110.gtf.json.gz \
  ${CDOT_DATA_DIR}/Ensembl/GRCh38/cdot-${CDOT_DATA_VERSION}-Homo_sapiens_GRCh38_Ensembl_111.gtf.json.gz \
  ${CDOT_DATA_DIR}/Ensembl/GRCh38/cdot-${CDOT_DATA_VERSION}-Homo_sapiens_GRCh38_Ensembl_112.gtf.json.gz \
  ${CDOT_DATA_DIR}/Ensembl/GRCh38/cdot-${CDOT_DATA_VERSION}-Homo_sapiens_GRCh38_Ensembl_113.gtf.json.gz \
  ${CDOT_DATA_DIR}/Ensembl/GRCh38/cdot-${CDOT_DATA_VERSION}-Homo_sapiens_GRCh38_Ensembl_114.gtf.json.gz \
  ${CDOT_DATA_DIR}/Ensembl/cdot-${CDOT_DATA_VERSION}-Ensembl-T2T-CHM13v2.0.json.gz \
  ${CDOT_DATA_DIR}/RefSeq/cdot-${CDOT_DATA_VERSION}-RefSeq-merged-GRCh37.json.gz \
  ${CDOT_DATA_DIR}/RefSeq/GRCh37/cdot-${CDOT_DATA_VERSION}-Homo_sapiens_GRCh37_RefSeq_105.20201022.gff.json.gz \
  ${CDOT_DATA_DIR}/RefSeq/GRCh37/cdot-${CDOT_DATA_VERSION}-Homo_sapiens_GRCh37_RefSeq_105.20220307.gff.json.gz \
  ${CDOT_DATA_DIR}/RefSeq/cdot-${CDOT_DATA_VERSION}-RefSeq-merged-GRCh38.json.gz \
  ${CDOT_DATA_DIR}/RefSeq/GRCh38/cdot-${CDOT_DATA_VERSION}-Homo_sapiens_GRCh38_RefSeq_110.gff.json.gz \
  ${CDOT_DATA_DIR}/RefSeq/GRCh38/cdot-${CDOT_DATA_VERSION}-Homo_sapiens_GRCh38_RefSeq_RS_2023_10.gff.json.gz \
  ${CDOT_DATA_DIR}/RefSeq/cdot-${CDOT_DATA_VERSION}-RefSeq-T2T-CHM13v2.0.json.gz

