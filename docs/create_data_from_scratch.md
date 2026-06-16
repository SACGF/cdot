# Create data from scratch

Most users use the REST API or download existing JSON files, so the data-generation dependencies are
not installed by default with the package.

## Install dependencies

```bash
sudo apt -y install postgresql-client-12  # Need psql to extract UTA transcripts
python3 -m pip install --upgrade htseq     # Need version after 2.0.2 to handle 109.20211119 GFF3
```

## Download and generate all transcripts

This requires ~3 GB of disk space and 4 GB of RAM.

```bash
export EMAIL=your@email.com            # change to your email - used for NCBI API calls
export CDOT_DATA_DIR=/data/gene_annotation  # change to your location

cd ${CDOT_DATA_DIR}
git clone https://github.com/SACGF/cdot   # generation scripts are not installed by the PyPI package
python3 -m pip install .                  # install the package

${CDOT_DATA_DIR}/cdot/generate_transcript_data/all_transcripts.sh
```

The transcript sources (Ensembl/RefSeq/UTA/T2T releases) are defined in
[`generate_transcript_data/cdot_transcripts.yaml`](../generate_transcript_data/cdot_transcripts.yaml).
