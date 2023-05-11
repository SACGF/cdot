#!/bin/bash

set -e

FULL_PATH_TO_SCRIPT="$(realpath "${BASH_SOURCE[-1]}")"
BASE_DIR=$(dirname ${FULL_PATH_TO_SCRIPT})

# Python scripts will import via generate_transcript_data
export PYTHONPATH=${BASE_DIR}/..

CDOT_VERSION=$(${BASE_DIR}/cdot_json.py --version)

# This needs to be passed to called bash scripts, so they are invoked with "." to use these variables
export REFSEQ_DIR=M_musculus
export SPECIES=Mus_musculus
export GENE_INFO_JSON=$(pwd)/${SPECIES}.gene-info-${CDOT_VERSION}.json.gz

if [[ ! -e ${GENE_INFO_JSON} ]]; then
  ${BASE_DIR}/gene_info.sh
fi

echo "Gene summary variable = ${GENE_INFO_JSON}"

# RefSeq
mkdir -p refseq
cd refseq

mkdir -p GRCm38
cd GRCm38
${BASE_DIR}/Mus_musculus/refseq_transcripts_grcm38.sh
cd ..

mkdir -p GRCm39
cd GRCm39
${BASE_DIR}/Mus_musculus/refseq_transcripts_grcm39.sh
cd ..
