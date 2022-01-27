#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Usage ${BASH_SOURCE[0]} uta_version genome_build"
  exit 1;
fi

UTA_VERSION=${1}
GENOME_BUILD=${2}

BASE_DIR=$(dirname ${BASH_SOURCE[0]})
export PGPASSWORD=anonymous

uta_csv_filename=uta_${UTA_VERSION}_${GENOME_BUILD}.csv
if [[ ! -e ${uta_csv_filename} ]]; then
  SQL=${BASE_DIR}/uta_${UTA_VERSION}_${GENOME_BUILD}.sql

  # can't have newlines in \copy command
  cat ${SQL} | tr -s '\n' ' ' | psql -h uta.invitae.com -U anonymous -d uta
fi

cdot_file="cdot.uta_${UTA_VERSION}.${GENOME_BUILD}.json.gz"
if [[ ! -e ${cdot_file} ]]; then
  POSTGRES_URL=postgresql://uta.invitae.com/uta_${UTA_VERSION}
  python3 ${BASE_DIR}/uta_csv_to_pyreference_json.py --uta-csv=${uta_csv_filename} --output ${cdot_file} --url ${POSTGRES_URL}
fi
