---
title: "hgvs: A Python package for manipulating sequence variants using HGVS nomenclature: 2018 Update"
authors: Meng Wang, Keith M Callenberg, Raymond Dalgleish, Alexandre Fedtsov, Naomi K Fox, Peter J Freeman, Kevin B Jacobs, Piotr Kaleta, Andrew J McMurry, Andreas Prlić, Veena Rajaraman, Reece K Hart
journal: Human Mutation
year: 2018
volume: 39
pages: 1803-1813
doi: 10.1002/humu.23615
pmid: 30129167
---

## Summary

Major update to the biocommons `hgvs` Python package (original: Hart 2015). Documents substantial additions since 2014: parsing/formatting/validating/normalising variants on genome, transcript, and protein sequences; projecting variants between aligned sequences including those with gapped alignments; flexible installation with remote or local data; extensive automated testing involving eight organisations. Validated against clinically relevant variants from ClinVar and HGMD.

## Relevance to cdot

This is the primary downstream consumer of cdot's data provider interface. The "gapped alignments" capability is directly relevant — cdot stores and serves CIGAR gap strings so that biocommons/hgvs can correctly handle transcripts that do not align perfectly to the reference genome. The paper is cited in cdot's own docs/notes.txt as the key reference.
