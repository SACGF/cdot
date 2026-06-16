# cdot

[![PyPi version](https://img.shields.io/pypi/v/cdot.svg)](https://pypi.org/project/cdot/) [![Python versions](https://img.shields.io/pypi/pyversions/cdot.svg)](https://pypi.org/project/cdot/) [![tests](https://github.com/SACGF/cdot/actions/workflows/tests.yml/badge.svg)](https://github.com/SACGF/cdot/actions/workflows/tests.yml) [![DOI](https://zenodo.org/badge/448753921.svg)](https://zenodo.org/doi/10.5281/zenodo.13324621)


**cdot** provides the transcript data needed to map and validate
[HGVS](http://varnomen.hgvs.org/) variants - the gene/transcript coordinates, exon structure and
genome alignments - for the two most popular Python HGVS libraries:
[biocommons HGVS](https://github.com/biocommons/hgvs) and [PyHGVS](https://github.com/counsyl/hgvs).

To do HGVS work (e.g. convert `NM_001637.3:c.1582G>A` to genomic coordinates) those libraries need a
transcript data source. The usual source, [UTA](https://github.com/biocommons/uta), is a PostgreSQL
database that's slow and heavy to run. cdot instead **converts the official RefSeq/Ensembl annotation
files (GTF/GFF3) into compact JSON** and ships fast loaders for the HGVS libraries. You can use it via:

* **Local `JSON.gz` files** - [download a release](https://github.com/SACGF/cdot/releases) and load it into RAM (fastest).
* **REST API** - query a [cdot_rest](https://github.com/SACGF/cdot_rest) server (no local data needed).

Because it reads the released annotation files directly, cdot covers **1.58 million transcript/genome
alignments**, including historical transcript versions - vs ~141k in UTA (v.20210129) - which matters
when resolving legacy HGVS. See [cdot vs UTA](docs/cdot_vs_uta.md) for the trade-offs.

Recent changes are in the [changelog](https://github.com/SACGF/cdot/blob/main/CHANGELOG.md).

## Install

```
pip install cdot
```

Optional extras:

```
pip install cdot[fasta]   # local genome FASTA sequence fetching (pysam) - needed for the PyHGVS example below
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

[more Biocommons examples](docs/examples_biocommons.md):

For fixing messy HGVS input and fast bulk processing, see [Advanced usage](docs/advanced_usage.md).

[PyHGVS](https://github.com/counsyl/hgvs) example (needs `pip install cdot[fasta]` for pysam):

```
import pyhgvs
from pysam.libcfaidx import FastaFile
from cdot.pyhgvs.pyhgvs_transcript import JSONPyHGVSTranscriptFactory, RESTPyHGVSTranscriptFactory

genome = FastaFile("/data/annotation/fasta/GCF_000001405.25_GRCh37.p13_genomic.fna.gz")
factory = RESTPyHGVSTranscriptFactory()
# factory = JSONPyHGVSTranscriptFactory(["./cdot-0.2.32.refseq.grch37.json.gz"])  # Uses local JSON file
pyhgvs.parse_hgvs_name('NM_001637.3:c.1582G>A', genome, get_transcript=factory.get_transcript_grch37)
```

[more PyHGVS examples](docs/examples_pyhgvs.md):

## Documentation

See [docs/](docs/) for reference and how-to guides:

* [JSON data format](docs/json_data_format.md) - every field in a cdot JSON(.gz) file
* [Coordinates & exon alignments](docs/coordinates_and_exons.md) - how exon coordinates and gap strings work
* [Advanced usage](docs/advanced_usage.md) - fixing messy HGVS input, and bulk read-ahead retrieval

See the [docs index](docs/README.md) for the full list (examples, FastaSeqFetcher, creating data, cdot vs UTA, …).

## Q. What's the performance like?

Resolving real ClinVar c.HGVS to genomic coordinates (GRCh38, biocommons HGVS, local sequence fetching):

* UTA public DB: ~1-1.5 seconds / transcript
* cdot REST service: ~30 HGVS/second sequential, **~500 HGVS/second** with [`prefetch()`](docs/advanced_usage.md) batch cache-warming
* cdot JSON.gz (local): 500-1k/second

`prefetch()` warms every transcript in one batch round-trip up front, so bulk resolution over the REST
service runs almost entirely from cache - closing most of the gap to local JSON.gz (~16x faster end-to-end
on 500 variants). Reproduce with `analysis/benchmark_resolution.py`.

## Q. Where can I download the JSON.gz files?

[Download from GitHub releases](https://github.com/SACGF/cdot/releases) - RefSeq (37/38) - 72M, Ensembl (37/38) 61M

Details on what the files contain [here](docs/release_files.md)

## Q. How does this compare to Universal Transcript Archive?

Both projects have similar goals of providing transcripts for loading HGVS, but they approach it from different ways

* UTA aligns sequences, then stores coordinates in an SQL database. 
* cdot convert existing Ensembl/RefSeq GTFs into JSON

See [cdot vs UTA](docs/cdot_vs_uta.md) for more details

## Q. How do you store transcripts in JSON?

See the **[JSON data format reference](docs/json_data_format.md)** for a full description of every field, with a machine-readable [JSON Schema](docs/cdot-json-schema.json) alongside it. [Coordinates & exon alignments](docs/coordinates_and_exons.md) explains how exon coordinates and the alignment gap strings work. See also [design notes](docs/design_notes.md) on why the format looks the way it does.

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
