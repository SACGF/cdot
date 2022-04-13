#!/bin/bash

set -e

FULL_PATH_TO_SCRIPT="$(realpath "${BASH_SOURCE[-1]}")"
BASE_DIR=$(dirname ${FULL_PATH_TO_SCRIPT})
CDOT_VERSION=$(${BASE_DIR}/cdot_json.py --version)

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

# Combine genome builds (we're in refseq dir)
${BASE_DIR}/cdot_json.py combine_builds \
    --grch37 GRCh37/cdot-${CDOT_VERSION}.refseq.grch37.json.gz \
    --grch38 GRCh38/cdot-${CDOT_VERSION}.refseq.grch38.json.gz \
    --output cdot-${CDOT_VERSION}.refseq.grch37_grch38.json.gz

cd ..

# RefSeq
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

# Combine genome builds (we're in ensembl dir)
${BASE_DIR}/cdot_json.py combine_builds \
    --grch37 GRCh37/cdot-${CDOT_VERSION}.ensembl.grch37.json.gz \
    --grch38 GRCh38/cdot-${CDOT_VERSION}.ensembl.grch38.json.gz \
    --output cdot-${CDOT_VERSION}.ensembl.grch37_grch38.json.gz

cd ..
