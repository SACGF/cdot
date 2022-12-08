#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Usage ${BASH_SOURCE[0]} uta_version genome_build"
  exit 1;
fi

BASE_DIR=$(dirname ${BASH_SOURCE[0]})
UTA_BASE_URL=uta.biocommons.org  # uta.invitae.com moved here
CDOT_VERSION=$(${BASE_DIR}/cdot_json.py --version)
UTA_VERSION=${1}
GENOME_BUILD=${2}

export PGPASSWORD=anonymous

uta_csv_filename=uta_${UTA_VERSION}_${GENOME_BUILD,,}.csv
if [[ ! -e ${uta_csv_filename} ]]; then
  SQL=${BASE_DIR}/uta_${UTA_VERSION}_${GENOME_BUILD,,}.sql  # Lowercase filename

  # can't have newlines in \copy command
  cat ${SQL} | tr -s '\n' ' ' | psql -h ${UTA_BASE_URL} -U anonymous -d uta
fi

cdot_file="cdot-${CDOT_VERSION}.uta_${UTA_VERSION}.${GENOME_BUILD}.json.gz"
if [[ ! -e ${cdot_file} ]]; then
  POSTGRES_URL=postgresql://${UTA_BASE_URL}/uta_${UTA_VERSION}
  ${BASE_DIR}/cdot_json.py uta_to_json "${uta_csv_filename}" --url "${POSTGRES_URL}" --output "${cdot_file}" --genome-build=${GENOME_BUILD}
fi
