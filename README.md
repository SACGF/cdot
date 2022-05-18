# cdot

[![PyPi version](https://img.shields.io/pypi/v/cdot.svg)](https://pypi.org/project/cdot/) [![Python versions](https://img.shields.io/pypi/pyversions/cdot.svg)](https://pypi.org/project/cdot/)

cdot provides transcripts for the 2 most popular Python [HGVS](http://varnomen.hgvs.org/) libraries.

It works by:

* Converting RefSeq/Ensembl GTFs to JSON 
* Providing loaders for the HGVS libraries, via JSON.gz files, or REST API via [cdot_rest](https://github.com/SACGF/cdot_rest))

We currently support ~893k transcripts (vs ~141k in UTA v.20210129)

## Install

```
pip install cdot
```

## Examples

[Biocommons HGVS](https://github.com/biocommons/hgvs) example:

```
from cdot.hgvs.dataproviders import JSONDataProvider, RESTDataProvider

hdp = RESTDataProvider()  # Uses API server at cdot.cc
# hdp = JSONDataProvider(["./cdot-0.2.6.refseq.grch38.json.gz"])  # Uses local JSON file

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
# factory = JSONPyHGVSTranscriptFactory(["./cdot-0.2.6.refseq.grch38.json.gz"])  # Uses local JSON file
pyhgvs.parse_hgvs_name(hgvs_c, genome, get_transcript=factory.get_transcript_grch37)

```

## Q. What's the performance like?

* UTA public DB: 1-1.5 seconds / transcript
* cdot REST service: 10/second
* cdot JSON.gz: 500-1k/second

## Q. Where can I download the JSON.gz files?

[RefSeq 37+38](https://cdot.cc/download/cdot-0.2.6.refseq.grch37_grch38.json.gz) - 70Mb
[Ensembl 37+38](https://cdot.cc/download/cdot-0.2.6.ensembl.grch37_grch38.json.gz) - 53Mb

See also [Download JSON.gz files](https://github.com/SACGF/cdot/wiki/Download-JSON.gz-files) if you only want individual builds.

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
