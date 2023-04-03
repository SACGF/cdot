# cdot

[![PyPi version](https://img.shields.io/pypi/v/cdot.svg)](https://pypi.org/project/cdot/) [![Python versions](https://img.shields.io/pypi/pyversions/cdot.svg)](https://pypi.org/project/cdot/)

cdot provides transcripts for the 2 most popular Python [HGVS](http://varnomen.hgvs.org/) libraries.

It works by:

* Converting RefSeq/Ensembl GTFs to JSON 
* Providing loaders for the HGVS libraries, via JSON.gz files, or REST API via [cdot_rest](https://github.com/SACGF/cdot_rest))

We currently support ~905k transcripts (vs ~141k in UTA v.20210129)

## New 

See [changelog](https://github.com/SACGF/cdot/blob/main/CHANGELOG.md)

2023-04-03:
* #41 - Support for T2T CHM13v2.0 [example code](https://github.com/SACGF/cdot/wiki/Biocommons-T2T-CHM13v2.0-example-code)

2023-03-31:
* #31 - Fasta file Biocommons HGVS SeqFetcher implementation
* #38 - bugfix for Biocommons HGVS get_tx_for_region

## Install

```
pip install cdot
```

## Examples

[Biocommons HGVS](https://github.com/biocommons/hgvs) example:

```
import hgvs
from hgvs.assemblymapper import AssemblyMapper
from cdot.hgvs.dataproviders import JSONDataProvider, RESTDataProvider

hdp = RESTDataProvider()  # Uses API server at cdot.cc
# hdp = JSONDataProvider(["./cdot-0.2.14.refseq.grch37.json.gz"])  # Uses local JSON file

am = AssemblyMapper(hdp,
                    assembly_name='GRCh37',
                    alt_aln_method='splign', replace_reference=True)

hp = hgvs.parser.Parser()
var_c = hp.parse_hgvs_variant('NM_001637.3:c.1582G>A')
am.c_to_g(var_c)
```

[more Biocommons examples](https://github.com/SACGF/cdot/wiki/Biocommons-HGVS-example-code):

[PyHGVS](https://github.com/counsyl/hgvs) example:

```
import pyhgvs
from pysam.libcfaidx import FastaFile
from cdot.pyhgvs.pyhgvs_transcript import JSONPyHGVSTranscriptFactory, RESTPyHGVSTranscriptFactory

genome = FastaFile("/data/annotation/fasta/GCF_000001405.25_GRCh37.p13_genomic.fna.gz")
factory = RESTPyHGVSTranscriptFactory()
# factory = JSONPyHGVSTranscriptFactory(["./cdot-0.2.14.refseq.grch37.json.gz"])  # Uses local JSON file
pyhgvs.parse_hgvs_name('NM_001637.3:c.1582G>A', genome, get_transcript=factory.get_transcript_grch37)
```

[more PyHGVS examples](https://github.com/SACGF/cdot/wiki/PyHGVS-example-code):

## Q. What's the performance like?

* UTA public DB: 1-1.5 seconds / transcript
* cdot REST service: 10/second
* cdot JSON.gz: 500-1k/second

## Q. Where can I download the JSON.gz files?

[Download from GitHub releases](https://github.com/SACGF/cdot/releases) - RefSeq (37/38) - 72M, Ensembl (37/38) 61M

Details on what the files contain [here](https://github.com/SACGF/cdot/wiki/GitHub-release-file-details)

## Q. How does this compare to Universal Transcript Archive?

Both projects have similar goals of providing transcripts for loading HGVS, but they approach it from different ways

* UTA aligns sequences, then stores coordinates in an SQL database. 
* cdot convert existing Ensembl/RefSeq GTFs into JSON

See [wiki for more details](https://github.com/SACGF/cdot/wiki/cdot-vs-UTA)

## Q. How do you store transcripts in JSON?

See [wiki page](https://github.com/SACGF/cdot/wiki/Transcript-JSON-format) for the format.

We think a standard for JSON gene/transcript information would be a great thing, and am keen to collaborate to make it happen!

## Q. What does cdot stand for?

cdot, pronounced "see dot" stands for Complete Dict of Transcripts

This was developed for the [Australian Genomics](https://www.australiangenomics.org.au/) [Shariant](https://shariant.org.au/) project, due to the need to load historical HGVS from lab archives.   
