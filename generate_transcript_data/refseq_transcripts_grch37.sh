#!/bin/bash

set -e

BASE_DIR=$(dirname ${BASH_SOURCE[0]})
CDOT_VERSION=$(${BASE_DIR}/cdot_json.py --version)
GENOME_BUILD=grch37
UTA_VERSION=20210129

# Having troubles with corrupted files downloading via FTP from NCBI via IPv6, http works ok
# NOTE: RefSeq transcripts in GRCh37 before p13 did not have alignment gap information

merge_args=()

filename=ref_GRCh37.p5_top_level.gff3.gz
url=http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/BUILD.37.3/GFF/${filename}
cdot_file=$(basename $filename .gz).json.gz
if [[ ! -e ${filename} ]]; then
  wget ${url}
fi
if [[ ! -e ${cdot_file} ]]; then
  ${BASE_DIR}/cdot_json.py gff3_to_json "${filename}" --url "${url}" --genome-build=GRCh37 --output "${cdot_file}"
fi
merge_args+=(${cdot_file})


filename=ref_GRCh37.p9_top_level.gff3.gz
url=http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.103/GFF/${filename}
cdot_file=$(basename $filename .gz).json.gz
if [[ ! -e ${filename} ]]; then
  wget ${url}
fi
if [[ ! -e ${cdot_file} ]]; then
  ${BASE_DIR}/cdot_json.py gff3_to_json "${filename}" --url "${url}" --genome-build=GRCh37 --output "${cdot_file}"
fi
merge_args+=(${cdot_file})


filename=ref_GRCh37.p10_top_level.gff3.gz
url=http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.104/GFF/${filename}
cdot_file=$(basename $filename .gz).json.gz
if [[ ! -e ${filename} ]]; then
  wget ${url}
fi
if [[ ! -e ${cdot_file} ]]; then
  ${BASE_DIR}/cdot_json.py gff3_to_json "${filename}" --url "${url}" --genome-build=GRCh37 --output "${cdot_file}"
fi
merge_args+=(${cdot_file})

# UTA transcripts have gaps, so they should overwrite the earlier refseq transcripts (without gaps)
# But will be overwritten by newer (post p13) official transcripts
cdot_file="cdot.uta_${UTA_VERSION}.${GENOME_BUILD}.json.gz"
${BASE_DIR}/uta_transcripts.sh ${UTA_VERSION} ${GENOME_BUILD}
merge_args+=(${cdot_file})


filename=ref_GRCh37.p13_top_level.gff3.gz
url=http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/${filename}
cdot_file=$(basename $filename .gz).json.gz
if [[ ! -e ${filename} ]]; then
  wget ${url}
fi
if [[ ! -e ${cdot_file} ]]; then
  ${BASE_DIR}/cdot_json.py gff3_to_json "${filename}" --url "${url}" --genome-build=GRCh37 --output "${cdot_file}"
fi
merge_args+=(${cdot_file})


# These all have the same name, so rename them based on release ID
for release in 105.20190906 105.20201022 105.20220307; do
  filename=GCF_000001405.25_GRCh37.p13_genomic.${release}.gff.gz
  url=http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/${release}/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz
  cdot_file=$(basename $filename .gz).json.gz
  if [[ ! -e ${filename} ]]; then
    wget ${url} --output-document=${filename}
  fi
  if [[ ! -e ${cdot_file} ]]; then
    ${BASE_DIR}/cdot_json.py gff3_to_json "${filename}" --url "${url}" --genome-build=GRCh37 --output "${cdot_file}"
  fi
  merge_args+=(${cdot_file})
done

merged_file="cdot-${CDOT_VERSION}.refseq.grch37.json.gz"
if [[ ! -e ${merged_file} ]]; then
  echo "Creating ${merged_file}"
  ${BASE_DIR}/cdot_json.py merge_historical ${merge_args[@]} --genome-build=GRCh37 --output "${merged_file}"
fi
