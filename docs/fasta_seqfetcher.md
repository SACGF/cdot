# FastaSeqFetcher — local FASTA sequence fetching

Needs the `fasta` extra: `pip install cdot[fasta]` (brings in `pysam`).

## Background

biocommons HGVS uses [SeqRepo](https://github.com/biocommons/biocommons.seqrepo) to retrieve genome
and transcript sequences.

cdot collects as many transcripts as possible, and some of these are not in SeqRepo — meaning you
can't resolve them using the [biocommons HGVS](https://github.com/biocommons/hgvs/) library.

`FastaSeqFetcher` is a SeqRepo replacement that builds transcript sequences by pasting together exons
from the genome. This means every cdot transcript becomes resolvable.

> **Warning: transcript sequences may differ from the genome sequence — do this at your own risk!**
> (See [Transcript vs genome differences](#warnings--transcript-vs-genome-differences) below.)

## Install

The FASTA files need to be bgzipped then indexed with samtools so they support random access. You
need `samtools` and `bgzip` installed.

To do this as you download them:

```bash
# GRCh37
wget --quiet -O - https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz | gzip -d | bgzip > GCF_000001405.25_GRCh37.p13_genomic.fna.gz
samtools faidx GCF_000001405.25_GRCh37.p13_genomic.fna.gz
```

```bash
# GRCh38
wget --quiet -O - https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz | gzip -d | bgzip > GCF_000001405.39_GRCh38.p13_genomic.fna.gz
samtools faidx GCF_000001405.39_GRCh38.p13_genomic.fna.gz
```

## Usage

Create an instance referencing a bgzipped, samtools-indexed FASTA, then pass it via the `seqfetcher`
keyword argument on the cdot data provider classes:

```python
from cdot.hgvs.dataproviders import JSONDataProvider
from cdot.hgvs.dataproviders.fasta_seqfetcher import FastaSeqFetcher

seqfetcher = FastaSeqFetcher("/data/fasta/GCF_000001405.39_GRCh38.p13_genomic.fna.gz")
hdp = JSONDataProvider(["./cdot-0.2.32.refseq.grch38.json.gz"], seqfetcher=seqfetcher)
```

`ChainedSeqFetcher` is a helper that lets you use the default `SeqFetcher` first and only fall back to
`FastaSeqFetcher` when no transcript sequence is available:

```python
from cdot.hgvs.dataproviders.fasta_seqfetcher import FastaSeqFetcher, ChainedSeqFetcher
from hgvs.dataproviders.seqfetcher import SeqFetcher

seqfetcher = ChainedSeqFetcher(SeqFetcher(), FastaSeqFetcher(fasta_filename))
```

## Warnings / Transcript vs genome differences

Ensembl transcripts always match the genome reference — this only affects RefSeq, where the
transcript and genome sequences can differ.

The errors that could occur from a genome/transcript mismatch are:

* Reference base changes.
* Indels normalised differently (the repeats may be different), leading to a coordinate change.

When choosing which genome reference to use, cdot takes the first one that has a contig we have a
mapping for. If you want to ensure a certain build is used, instantiate `FastaSeqFetcher()` with just
that one genome, then only convert HGVS using that genome.

It's not clear how often this will fail. The `test/benchmark_hgvs.py` script uses 500 random ClinVar
matched g.HGVS and c.HGVS records, and all converted without any issues.
