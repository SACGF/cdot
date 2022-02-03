# cdot

[![PyPi version](https://img.shields.io/pypi/v/cdot.svg)](https://pypi.org/project/cdot/) [![Python versions](https://img.shields.io/pypi/pyversions/cdot.svg)](https://pypi.org/project/cdot/)

cdot provides transcripts for the 2 most popular Python [HGVS](http://varnomen.hgvs.org/) libraries.

It works by:

* Converting RefSeq/Ensembl GTFs to JSON
* Providing loaders from JSON.gz files, or REST API via [cdot_rest](https://github.com/SACGF/cdot_rest))

We currently support ~800k transcripts, with API responses under 0.1 second

## Examples

[Biocommons HGVS](https://github.com/biocommons/hgvs) example:

```
from cdot.hgvs.dataproviders import JSONDataProvider, RESTDataProvider

hdp = RESTDataProvider()  # Uses API server at cdot.cc
# hdp = JSONDataProvider(["./cdot-0.2.1.refseq.grch38.json.gz"])  # Uses local JSON file

am = AssemblyMapper(hdp,
                    assembly_name='GRCh37',
                    alt_aln_method='splign', replace_reference=True)

hp = hgvs.parser.Parser()
var_c = hp.parse_hgvs_variant('NM_001637.3:c.1582G>A')
am.c_to_g(var_c)
```

[PyHGVS](https://github.com/counsyl/hgvs) example:

```
from cdot.pyhgvs.pyhgvs_transcript import JSONPyHGVSTranscriptFactory, RESTPyHGVSTranscriptFactory

factory = RESTPyHGVSTranscriptFactory()
# factory = JSONPyHGVSTranscriptFactory(["./cdot-0.2.1.refseq.grch38.json.gz"])  # Uses local JSON file
pyhgvs.parse_hgvs_name(hgvs_c, genome, get_transcript=factory.get_transcript_grch37)

```

## Download data

TODO

## Philosophical differences from Universal Transcript Archive

cdot aims to be as simple as possible: convert existing Ensembl/RefSeq GTFs into JSON format

Universal transcript archive is an excellent and ambitious project that:

* Performs its own mapping of transcript sequences to reference genomes
* Stores the transcript version data (exons etc) in a SQL database

This has some advantages, namely that you can resolve a GRCh37 coordinate for a transcript which was never officially released for that build.

However the complexity causes a few downsides:

* Alignments may not exactly match those in official Ensembl/RefSeq releases
* Local install requires a PostgreSQL installation
* Internet hosted UTA is a PostgreSQL server, so requires client Postgres libraries, is inaccessible behind firewalls. They have been planning on building a [REST server since 2014](https://github.com/biocommons/uta/issues/164)
* High complexity manual process for releases means they [do not support Ensembl](https://github.com/biocommons/uta/issues/233) and take a while to make RefSeq releases.

## What does cdot stand for?

cdot, pronounced "see dot" stands for Complete Dict of Transcripts
