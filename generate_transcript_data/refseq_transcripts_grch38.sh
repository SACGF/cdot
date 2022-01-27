#!/bin/bash

set -e

BASE_DIR=$(dirname ${BASH_SOURCE[0]})
GENOME_BUILD=grch38
UTA_VERSION=20210129

# Having troubles with corrupted files downloading via FTP from NCBI via IPv6, http works ok

pyreference_args=()

# All GRCh38 transcripts have alignments gaps, so use UTA first (and override with official releases)
uta_cdot_file="cdot.uta_${UTA_VERSION}.${GENOME_BUILD}.json.gz"
${BASE_DIR}/uta_transcripts.sh ${UTA_VERSION} ${GENOME_BUILD}
pyreference_args+=(--pyreference-json ${uta_cdot_file})


filename=ref_GRCh38_top_level.gff3.gz
url=http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.106/GFF/${filename}
pyreference_file=$(basename $filename .gz).json.gz

if [[ ! -e ${filename} ]]; then
  wget ${url}
fi
if [[ ! -e ${pyreference_file} ]]; then
  pyreference_gff_to_json.py --url "${url}" --gff3 "${filename}"
fi
pyreference_args+=(--pyreference-json ${pyreference_file})


filename=ref_GRCh38.p2_top_level.gff3.gz
url=http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.107/GFF/${filename}
pyreference_file=$(basename $filename .gz).json.gz

if [[ ! -e ${filename} ]]; then
  wget ${url}
fi
if [[ ! -e ${pyreference_file} ]]; then
  pyreference_gff_to_json.py --url "${url}" --gff3 "${filename}"
fi
pyreference_args+=(--pyreference-json ${pyreference_file})


filename=ref_GRCh38.p7_top_level.gff3.gz
url=http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.108/GFF/${filename}
pyreference_file=$(basename $filename .gz).json.gz

if [[ ! -e ${filename} ]]; then
  wget ${url}
fi
if [[ ! -e ${pyreference_file} ]]; then
  pyreference_gff_to_json.py --url "${url}" --gff3 "${filename}"
fi
pyreference_args+=(--pyreference-json ${pyreference_file})


filename=ref_GRCh38.p12_top_level.gff3.gz
url=http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.109/GFF/${filename}
pyreference_file=$(basename $filename .gz).json.gz

if [[ ! -e ${filename} ]]; then
  wget ${url}
fi
if [[ ! -e ${pyreference_file} ]]; then
  pyreference_gff_to_json.py --url "${url}" --gff3 "${filename}"
fi
pyreference_args+=(--pyreference-json ${pyreference_file})


filename=GCF_000001405.38_GRCh38.p12_genomic.gff.gz
url=http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109/GCF_000001405.38_GRCh38.p12/${filename}
pyreference_file=$(basename $filename .gz).json.gz

if [[ ! -e ${filename} ]]; then
  wget ${url}
fi
if [[ ! -e ${pyreference_file} ]]; then
  pyreference_gff_to_json.py --url "${url}" --gff3 "${filename}"
fi
pyreference_args+=(--pyreference-json ${pyreference_file})


# These all have the same name, so rename them based on release ID
# 28 Jan 2022 - HTSeq GFF Parser currently dies on 109.20211119 - removed, waiting for it to get fixed
for release in 109.20190607 109.20190905 109.20191205 109.20200228 109.20200522 109.20200815 109.20201120 109.20210226 109.20210514; do
  filename=GCF_000001405.39_GRCh38.p13_genomic.${release}.gff.gz
  url=http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/${release}/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz
  pyreference_file=$(basename $filename .gz).json.gz
  if [[ ! -e ${filename} ]]; then
    wget ${url} --output-document=${filename}
  fi
  if [[ ! -e ${pyreference_file} ]]; then
    pyreference_gff_to_json.py --url "${url}" --gff3 "${filename}"
  fi
  pyreference_args+=(--pyreference-json ${pyreference_file})
done

merged_file="cdot-$(date --iso).refseq.grch38.json.gz"
if [[ ! -e ${merged_file} ]]; then
  echo "Creating ${merged_file}"
  python3 ${BASE_DIR}/pyreference_to_cdot_json.py ${pyreference_args[@]} --output ${merged_file}
fi
