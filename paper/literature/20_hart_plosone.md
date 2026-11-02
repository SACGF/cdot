---
title: "SeqRepo: A system for managing local collections of biological sequences"
authors: Reece K Hart, Andreas Prlić
journal: PLoS ONE
year: 2020
volume: 15
pages: e0239883
doi: 10.1371/journal.pone.0239883
pmid: 33270643
---

## Summary

Introduces SeqRepo, a system for building and managing local high-performance non-redundant biological sequence collections. Allows sequence retrieval by standard database identifiers (NCBI, Ensembl, GRCh, LRG) and digest identifiers (sha512t24u, GA4GH). Provides Python API and REST interface with GA4GH refget protocol support. Demonstrates up to 1300× performance improvement over remote sequence retrieval, with 50× throughput gain for variant validation pipelines.

## Relevance to cdot

SeqRepo is the sequence storage companion to biocommons/hgvs. While cdot provides transcript *coordinate* data (exon positions, CDS coordinates), SeqRepo stores the actual *nucleotide sequences* needed for variant validation and normalisation. cdot's `FastaSeqFetcher` is an alternative to SeqRepo for users who have local FASTA files — it provides similar offline capability without the SeqRepo infrastructure. The two projects are complementary: cdot extends transcript coverage; SeqRepo extends sequence coverage. The 1300× speedup claim parallels cdot's own speed advantage over remote UTA access.
