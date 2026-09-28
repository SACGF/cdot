# Create data from scratch

Most users use the REST API or download existing JSON files, so the data-generation dependencies are
not installed by default with the package. The generation scripts live in
[`generate_transcript_data/`](../generate_transcript_data/) and are not part of the PyPI package, so
you need a checkout of the repo.

## Install dependencies

```bash
git clone https://github.com/SACGF/cdot
cd cdot
python3 -m pip install -r generate_transcript_data/requirements.txt "snakemake>=8,<9"
sudo apt -y install postgresql-client   # psql, only needed for the UTA sources
```

## Convert a single GTF/GFF3

`cdot_json.py` converts one annotation file at a time. The gene info JSON (Entrez gene summaries)
is built once from the NCBI `gene_info` file:

```bash
export PYTHONPATH=$(pwd)  # so the scripts can import each other
curl -L -o Homo_sapiens.gene_info.gz https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/Homo_sapiens.gene_info.gz
generate_transcript_data/cdot_gene_info.py --gene-info Homo_sapiens.gene_info.gz --output gene_info.json.gz --email your@email.com

# RefSeq GFF3
generate_transcript_data/cdot_json.py gff3_to_json GCF_000001405.40_GRCh38.p14_genomic.gff.gz \
    --url "https://ftp.ncbi.nlm.nih.gov/.../GCF_000001405.40_GRCh38.p14_genomic.gff.gz" \
    --genome-build GRCh38 --gene-info-json gene_info.json.gz --output cdot.refseq.grch38.json.gz

# Ensembl GTF (Ensembl GFF3 is not supported, its CDS rows have no protein version)
generate_transcript_data/cdot_json.py gtf_to_json Homo_sapiens.GRCh38.116.gtf.gz \
    --url "https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz" \
    --genome-build GRCh38 --gene-info-json gene_info.json.gz --output cdot.ensembl.grch38.json.gz
```

Run `cdot_json.py gtf_to_json --help` for the other options, including `--annotation-consortium`
(RefSeq/Ensembl conventions are otherwise detected from the file) and `--no-contig-conversion` for
genome builds without a biocommons bioutils assembly (the output then won't work with biocommons HGVS).

## Build all the release files

The [Snakemake pipeline](../generate_transcript_data/Snakefile) downloads every source in
[`cdot_transcripts.yaml`](../generate_transcript_data/cdot_transcripts.yaml) (Ensembl GTF, RefSeq
GFF3, UTA and T2T-CHM13v2.0), converts each, then merges them into the per-build and all-builds
files that are published as data releases. It needs roughly 15 GB of disk and 6 GB of RAM, and a
full run takes hours (the official releases are built on GitHub Actions, see
[building a data release](data_release_workflow.md)).

```bash
mkdir data && cd data
snakemake --snakefile ../generate_transcript_data/Snakefile --cores 4
```

Outputs land under `data/refseq/` and `data/ensembl/`. To build a single file, name it as the target,
eg `refseq/GRCh38/cdot-<data version>.Homo_sapiens_GRCh38_RefSeq_RS_2025_08.gff.json.gz`
(`cdot_json.py --version` prints the data version).
