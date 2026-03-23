---
title: "A joint NCBI and EMBL-EBI transcript set for clinical genomics and research"
authors: Joannella Morales, Shashikant Pujar, Jane E Loveland, Alex Astashyn, et al. (Flicek, Birney, Pruitt, Frankish, Cunningham, Murphy)
journal: Nature
year: 2022
volume: 604
pages: 310-315
doi: 10.1038/s41586-022-04558-8
pmid: 35388217
---

## Summary

Establishes the Matched Annotation from NCBI and EMBL-EBI (MANE) transcript sets: MANE Select (one representative transcript per human protein-coding gene) and MANE Plus Clinical (additional transcripts for reporting known clinical variants). Achieves MANE Select coverage for 97% of human protein-coding genes, including all ACMG secondary findings genes. Unified annotation removes ambiguity arising from RefSeq and Ensembl historically using different transcripts for the same gene.

## Relevance to cdot

cdot already stores MANE Select and MANE Plus Clinical tags in transcript build-specific data (since v0.2.12), enabling cdot users to identify canonical transcripts. The cdot issue #36 (canonical transcripts) and docs/todo.txt both cite MANE as a desired lookup feature. The paper is strong motivation for why canonical transcript lookup is a high-value feature to implement. cdot currently supports MANE for GRCh38 (where tags appear in RefSeq/Ensembl GFFs) but not GRCh37 (MANE is GRCh38-only by design).
