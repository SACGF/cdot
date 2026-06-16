# cdot vs UTA

cdot and the [Universal Transcript Archive (UTA)](https://github.com/biocommons/uta) have similar
goals of providing transcripts for loading HGVS, but they approach it in different ways:

* UTA aligns sequences, then stores coordinates in an SQL database.
* cdot converts existing Ensembl/RefSeq GTFs into JSON.

## Alignment gaps

RefSeq transcript sequences can differ from the genome sequence, which means they can
[align with gaps](https://hgvs.readthedocs.io/en/stable/examples/using-hgvs.html#projecting-in-the-presence-of-a-genome-transcript-gap).
Prior to v105 (GRCh37.p13) RefSeq did not provide alignment gap information, so UTA was forced to do
their own alignment to get CIGAR strings, in order to correctly handle these gaps.

From v105 onwards, RefSeq provide
[these gaps](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md#the-gap-attribute) -
making it possible to use the GFFs. See [Coordinates & exon alignments](coordinates_and_exons.md) for
how cdot stores them.

## Advantages of aligning sequences (UTA)

* UTA can map GRCh37 sequences to GRCh38 and vice-versa.
* UTA can account for alignment gaps in earlier RefSeq releases (cdot uses these UTA transcripts - thanks!).

## Advantages of using existing GTFs (cdot)

* Drastically simpler workflow - meaning we can load more transcripts.
* Alignments exactly match those in the official releases.

## JSON vs SQL

There's a bit of redundancy in JSON, but:

* You can copy flat files around without dealing with Docker/PostgreSQL/database schemas etc.
* It's trivial to write a [REST server](https://github.com/SACGF/cdot_rest) and the client already consumes JSON.
* It's lightning fast to load into RAM in Python.
