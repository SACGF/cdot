# cdot

[![PyPi version](https://img.shields.io/pypi/v/cdot.svg)](https://pypi.org/project/cdot/) [![Python versions](https://img.shields.io/pypi/pyversions/cdot.svg)](https://pypi.org/project/cdot/) [![tests](https://github.com/SACGF/cdot/actions/workflows/tests.yml/badge.svg)](https://github.com/SACGF/cdot/actions/workflows/tests.yml) [![DOI](https://zenodo.org/badge/448753921.svg)](https://zenodo.org/doi/10.5281/zenodo.13324621)


**cdot** provides the transcript data needed to map and validate
[HGVS](http://varnomen.hgvs.org/) variants - the gene/transcript coordinates, exon structure and
genome alignments - for the two most popular Python HGVS libraries:
[biocommons HGVS](https://github.com/biocommons/hgvs) and [PyHGVS](https://github.com/counsyl/hgvs).
It also provides tools to repair the malformed or outdated HGVS strings common in real-world data.

To do HGVS work (e.g. convert `NM_001637.3:c.1582G>A` to genomic coordinates) those libraries need a
transcript data source. The usual source, [UTA](https://github.com/biocommons/uta), is a PostgreSQL
database that's slow and heavy to run. cdot instead **converts the official RefSeq/Ensembl annotation
files (GTF/GFF3) into compact JSON** and ships fast loaders for the HGVS libraries. You can use it via:

* **Local `JSON.gz` files** - [download a release](https://github.com/SACGF/cdot/releases) and load it into RAM (fastest).
* **REST API** - query a [cdot_rest](https://github.com/SACGF/cdot_rest) server (no local data needed).

Because it reads the released annotation files directly, cdot covers **1.58 million transcript/genome
alignments**, including historical transcript versions - vs ~141k in UTA (v.20210129) - which matters
when resolving legacy HGVS. See [cdot vs UTA](https://github.com/SACGF/cdot/blob/main/docs/cdot_vs_uta.md) for the trade-offs.

Transcript data alone doesn't get real-world HGVS resolved. Strings pasted from reports, spreadsheets and
old lab archives are often slightly malformed, name a gene instead of a transcript, or cite a transcript
version that has since been retired. cdot also includes tools for this in `cdot.hgvs` (see
[Fixing real-world HGVS](#fixing-real-world-hgvs)):

* **Cleaning** - `clean_hgvs()` repairs common formatting errors (whitespace, case, punctuation typos,
  swapped gene/transcript, trailing `p.` annotations) and reports every change it made.
* **Gene symbols** - `BRCA2:c.68_69del` resolves to the MANE Select transcript.
* **Version substitution** - opt-in: if the cited transcript version isn't available, use an adjacent
  version, but only when a check confirms the variant's genomic coordinate won't change.

Recent changes are in the [changelog](https://github.com/SACGF/cdot/blob/main/CHANGELOG.md), and changes to the published transcript data in the [data changelog](https://github.com/SACGF/cdot/blob/main/CHANGELOG-data.md).

## Install

```
pip install cdot
```

Optional extras:

```
pip install 'cdot[fasta]'   # local genome FASTA sequence fetching (pysam) - needed for the PyHGVS example below. Quotes are required in zsh (macOS default).
```

(`hgvs` is a core dependency, so the biocommons HGVS examples work out of the box.)

## Examples

[Biocommons HGVS](https://github.com/biocommons/hgvs) example:

```
import hgvs
from hgvs.assemblymapper import AssemblyMapper
from cdot.hgvs.dataproviders import JSONDataProvider, RESTDataProvider

hdp = RESTDataProvider()  # Uses API server at cdotlib.org
# hdp = JSONDataProvider(["./cdot-0.2.32.refseq.grch37.json.gz"])  # Uses local JSON file

am = AssemblyMapper(hdp,
                    assembly_name='GRCh37',
                    alt_aln_method='splign', replace_reference=True)

hp = hgvs.parser.Parser()
var_c = hp.parse_hgvs_variant('NM_001637.3:c.1582G>A')
am.c_to_g(var_c)
```

[more Biocommons examples](https://github.com/SACGF/cdot/blob/main/docs/examples_biocommons.md):

> **Tip:** cdot provides many transcripts that aren't in SeqRepo, so the default sequence fetcher will
> raise `HGVSDataNotAvailableError` for them. You almost always want to supply a
> [`FastaSeqFetcher`](https://github.com/SACGF/cdot/blob/main/docs/fasta_seqfetcher.md) (chained after SeqRepo) so every cdot transcript
> resolves against a local genome FASTA.

For fast bulk processing over the REST API, see [read-ahead batch retrieval](https://github.com/SACGF/cdot/blob/main/docs/advanced_usage.md#read-ahead-batch-retrieval).

[PyHGVS](https://github.com/counsyl/hgvs) example (needs `pip install 'cdot[fasta]'` for pysam):

```
import pyhgvs
from pysam.libcfaidx import FastaFile
from cdot.pyhgvs.pyhgvs_transcript import JSONPyHGVSTranscriptFactory, RESTPyHGVSTranscriptFactory

genome = FastaFile("/data/annotation/fasta/GCF_000001405.25_GRCh37.p13_genomic.fna.gz")
factory = RESTPyHGVSTranscriptFactory()
# factory = JSONPyHGVSTranscriptFactory(["./cdot-0.2.32.refseq.grch37.json.gz"])  # Uses local JSON file
pyhgvs.parse_hgvs_name('NM_001637.3:c.1582G>A', genome, get_transcript=factory.get_transcript_grch37)
```

[more PyHGVS examples](https://github.com/SACGF/cdot/blob/main/docs/examples_pyhgvs.md):

## Fixing real-world HGVS

biocommons HGVS rejects strings that are almost, but not quite, valid. `fix_hgvs()` cleans them up first
and returns a list of `HGVSFix` records (severity, code, message) so you can decide whether to trust,
log or reject the result:

```python
from cdot.hgvs import fix_hgvs

result, fixes = fix_hgvs("BRCA2(NM_000059.4):c.68_69DEL p.(Glu23fs)")
# result = "NM_000059.4(BRCA2):c.68_69del"
for fix in fixes:
    print(fix)
# Removed trailing protein (p.) annotation
# Lowercased mutation type 'DEL'
# Swapped gene/transcript
```

Cleaning is a pure string operation, needing no data or parser, so you can run it ahead of any HGVS
library. Give `fix_hgvs()` a data provider and genome build and it can also resolve gene symbols and
(if you opt in) substitute a missing transcript version:

```python
from cdot.hgvs import fix_hgvs, VersionStrategy

# Gene symbol -> MANE Select transcript
result, fixes = fix_hgvs("BRCA2:c.68_69del", hdp, "GRCh38")
# result = "NM_000059.4:c.68_69del"

# Retired transcript version -> adjacent version, only if coordinate-safe
result, fixes = fix_hgvs("NM_000059.2:c.36del", hdp, "GRCh38",
                         version_fallback=VersionStrategy.UP_THEN_DOWN)
```

Version substitution is off by default. When it is on, a substitution that can't be verified
coordinate-safe (same coding structure, alignment gaps and relevant UTR length as the cited version) is
refused by default, returning an ERROR fix and leaving the string unchanged. See [Advanced usage](https://github.com/SACGF/cdot/blob/main/docs/advanced_usage.md)
for all the options, and [transcript-version safety](https://github.com/SACGF/cdot/blob/main/docs/transcript_version_safety.md)
for how the safety check works and the evidence behind it.

## Documentation

See [docs/](https://github.com/SACGF/cdot/tree/main/docs) for reference and how-to guides:

* [JSON data format](https://github.com/SACGF/cdot/blob/main/docs/json_data_format.md) - every field in a cdot JSON(.gz) file
* [Coordinates & exon alignments](https://github.com/SACGF/cdot/blob/main/docs/coordinates_and_exons.md) - how exon coordinates and gap strings work
* [Advanced usage](https://github.com/SACGF/cdot/blob/main/docs/advanced_usage.md) - fixing messy HGVS input, and bulk read-ahead retrieval
* [Transcript-version safety](https://github.com/SACGF/cdot/blob/main/docs/transcript_version_safety.md) - when substituting a different transcript version is coordinate-safe

See the [docs index](https://github.com/SACGF/cdot/blob/main/docs/README.md) for the full list (examples, FastaSeqFetcher, creating data, cdot vs UTA, …).

## Q. What's the performance like?

Resolving real ClinVar c.HGVS to genomic coordinates (GRCh38, biocommons HGVS, local sequence fetching):

* UTA public DB: ~1-1.5 seconds / transcript
* cdot REST service: ~30 HGVS/second sequential, **~500 HGVS/second** with [`prefetch()`](https://github.com/SACGF/cdot/blob/main/docs/advanced_usage.md) batch cache-warming
* cdot JSON.gz (local): 500-1k/second

`prefetch()` warms every transcript in one batch round-trip up front, so bulk resolution over the REST
service runs almost entirely from cache - closing most of the gap to local JSON.gz (~16x faster end-to-end
on 500 variants). Reproduce with `paper/scripts/benchmark_resolution.py`.

## Q. Where can I download the JSON.gz files?

[Download from GitHub releases](https://github.com/SACGF/cdot/releases) - RefSeq (37/38) - 72M, Ensembl (37/38) 61M

Details on what the files contain [here](https://github.com/SACGF/cdot/blob/main/docs/release_files.md)

## Q. How does this compare to Universal Transcript Archive?

Both projects have similar goals of providing transcripts for loading HGVS, but they approach it from different ways

* UTA aligns sequences, then stores coordinates in an SQL database. 
* cdot convert existing Ensembl/RefSeq GTFs into JSON

See [cdot vs UTA](https://github.com/SACGF/cdot/blob/main/docs/cdot_vs_uta.md) for more details

## Q. How do you store transcripts in JSON?

See the **[JSON data format reference](https://github.com/SACGF/cdot/blob/main/docs/json_data_format.md)** for a full description of every field, with a machine-readable [JSON Schema](https://github.com/SACGF/cdot/blob/main/docs/cdot-json-schema.json) alongside it. [Coordinates & exon alignments](https://github.com/SACGF/cdot/blob/main/docs/coordinates_and_exons.md) explains how exon coordinates and the alignment gap strings work. See also [design notes](https://github.com/SACGF/cdot/blob/main/docs/design_notes.md) on why the format looks the way it does.

You can also read the data with typed Python objects (no extra install required):

```python
from cdot import models

data = models.load("cdot-0.2.32.refseq.GRCh38.json.gz")
tx = data.transcripts["NM_001637.3"]
print(tx.gene_name, tx.protein)
```

We think a standard for JSON gene/transcript information would be a great thing, and am keen to collaborate to make it happen!

## Q. What does cdot stand for?

cdot, pronounced "see dot", is a play on the HGVS coding-sequence prefix ```:c.```

This was developed for the [Australian Genomics](https://www.australiangenomics.org.au/) [Shariant](https://shariant.org.au/) project, due to the need to load historical HGVS from lab archives.   
