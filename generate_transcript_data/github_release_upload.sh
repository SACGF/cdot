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

files=(
  "${CDOT_DATA_DIR}/ensembl/cdot-${CDOT_DATA_VERSION}.ensembl.GRCh37.json.gz"
  "${CDOT_DATA_DIR}/ensembl/GRCh37/cdot-${CDOT_DATA_VERSION}.Homo_sapiens_GRCh37_Ensembl_87.gtf.json.gz"
  "${CDOT_DATA_DIR}/ensembl/cdot-${CDOT_DATA_VERSION}.ensembl.GRCh38.json.gz"
  "${CDOT_DATA_DIR}/ensembl/GRCh38/cdot-${CDOT_DATA_VERSION}.Homo_sapiens_GRCh38_Ensembl_112.gtf.json.gz"
  "${CDOT_DATA_DIR}/ensembl/GRCh38/cdot-${CDOT_DATA_VERSION}.Homo_sapiens_GRCh38_Ensembl_113.gtf.json.gz"
  "${CDOT_DATA_DIR}/ensembl/GRCh38/cdot-${CDOT_DATA_VERSION}.Homo_sapiens_GRCh38_Ensembl_114.gtf.json.gz"
  "${CDOT_DATA_DIR}/ensembl/GRCh38/cdot-${CDOT_DATA_VERSION}.Homo_sapiens_GRCh38_Ensembl_115.gtf.json.gz"
  "${CDOT_DATA_DIR}/ensembl/cdot-${CDOT_DATA_VERSION}.ensembl.T2T-CHM13v2.0.json.gz"
  "${CDOT_DATA_DIR}/refseq/cdot-${CDOT_DATA_VERSION}.all-builds-refseq-grch37_grch38_t2t-chm13v2.0.json.gz"
  "${CDOT_DATA_DIR}/refseq/cdot-${CDOT_DATA_VERSION}.refseq.GRCh37.json.gz"
  "${CDOT_DATA_DIR}/refseq/GRCh37/cdot-${CDOT_DATA_VERSION}.Homo_sapiens_GRCh37_RefSeq_105.20201022.gff.json.gz"
  "${CDOT_DATA_DIR}/refseq/GRCh37/cdot-${CDOT_DATA_VERSION}.Homo_sapiens_GRCh37_RefSeq_105.20220307.gff.json.gz"
  "${CDOT_DATA_DIR}/refseq/cdot-${CDOT_DATA_VERSION}.refseq.GRCh38.json.gz"
  "${CDOT_DATA_DIR}/refseq/GRCh38/cdot-${CDOT_DATA_VERSION}.Homo_sapiens_GRCh38_RefSeq_110.gff.json.gz"
  "${CDOT_DATA_DIR}/refseq/GRCh38/cdot-${CDOT_DATA_VERSION}.Homo_sapiens_GRCh38_RefSeq_RS_2025_08.gff.json.gz"
  "${CDOT_DATA_DIR}/refseq/cdot-${CDOT_DATA_VERSION}.refseq.T2T-CHM13v2.0.json.gz"
)

git fetch --tags
if ! git tag -l "${CDOT_RELEASE_NAME}" | grep -q .; then
  echo "Git repo has no tag of '${CDOT_RELEASE_NAME}'"
  exit 1
fi


if gh release view ${CDOT_RELEASE_NAME}; then
  echo "Release ${CDOT_RELEASE_NAME} exists"
else
  echo "Creating release ${CDOT_RELEASE_NAME}"
  echo "Generating notes..."
  RELEASE_NOTES_FILENAME="/tmp/${CDOT_RELEASE_NAME}.txt"
  echo > ${RELEASE_NOTES_FILENAME} # Clear
  for f in "${files[@]}"; do
    ${BASE_DIR}/cdot_json.py release_notes $f --show-urls >> ${RELEASE_NOTES_FILENAME}
    echo "" >> ${RELEASE_NOTES_FILENAME}  # New line
  done
  echo "Creating on GitHub"
  gh release create ${CDOT_RELEASE_NAME} --title=${CDOT_RELEASE_NAME} --notes-file="${RELEASE_NOTES_FILENAME}"
fi

gh release upload "${CDOT_RELEASE_NAME}" "${files[@]}"
