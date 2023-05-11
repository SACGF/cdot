#!/bin/bash

set -e

FULL_PATH_TO_SCRIPT="$(realpath "${BASH_SOURCE[-1]}")"
BASE_DIR=$(dirname ${FULL_PATH_TO_SCRIPT})

# Python scripts will import via generate_transcript_data
export PYTHONPATH=${BASE_DIR}/..

CDOT_VERSION=$(${BASE_DIR}/cdot_json.py --version)


# This needs to be passed to called bash scripts, so they are invoked with "." to use these variables
export GENE_INFO_JSON=$(pwd)/Homo_sapiens.gene-info-${CDOT_VERSION}.json.gz

if [[ ! -e ${GENE_INFO_JSON} ]]; then
  ${BASE_DIR}/gene_info.sh
fi

echo "Gene summary variable = ${GENE_INFO_JSON}"

# RefSeq
mkdir -p refseq
cd refseq

mkdir -p GRCh37
cd GRCh37
${BASE_DIR}/refseq_transcripts_grch37.sh
cd ..

mkdir -p GRCh38
cd GRCh38
${BASE_DIR}/refseq_transcripts_grch38.sh
cd ..

mkdir -p CHM13v2.0
cd CHM13v2.0
${BASE_DIR}/refseq_transcripts_chm13v2.sh
cd ..

# Combine genome builds (we're in refseq dir)
REFSEQ_COMBO=cdot-${CDOT_VERSION}.refseq.grch37_grch38.json.gz
if [[ ! -e ${REFSEQ_COMBO} ]]; then
  ${BASE_DIR}/cdot_json.py combine_builds \
      --grch37 GRCh37/cdot-${CDOT_VERSION}.refseq.grch37.json.gz \
      --grch38 GRCh38/cdot-${CDOT_VERSION}.refseq.grch38.json.gz \
      --output ${REFSEQ_COMBO}
fi

cd ..

# Ensembl
mkdir -p ensembl
cd ensembl

mkdir -p GRCh37
cd GRCh37
${BASE_DIR}/ensembl_transcripts_grch37.sh
cd ..

mkdir -p GRCh38
cd GRCh38
${BASE_DIR}/ensembl_transcripts_grch38.sh
cd ..

mkdir -p CHM13v2.0
cd CHM13v2.0
${BASE_DIR}/ensembl_transcripts_chm13v2.sh
cd ..


# Combine genome builds (we're in ensembl dir)
ENSEMBL_COMBO=cdot-${CDOT_VERSION}.ensembl.grch37_grch38.json.gz
if [[ ! -e ${ENSEMBL_COMBO} ]]; then
  ${BASE_DIR}/cdot_json.py combine_builds \
      --grch37 GRCh37/cdot-${CDOT_VERSION}.ensembl.grch37.json.gz \
      --grch38 GRCh38/cdot-${CDOT_VERSION}.ensembl.grch38.json.gz \
      --output ${ENSEMBL_COMBO}
fi

cd ..
