---
title: "Choice of transcripts and software has a large effect on variant annotation"
authors: Davis J McCarthy, Peter Humburg, Alexander Kanapin, Manuel A Rivas, Kyle Gaulton, Jean-Baptiste Cazier, Peter Donnelly
journal: Genome Medicine
year: 2014
volume: 6
pages: 26
doi: 10.1186/gm543
pmid: 24944579
---

## Summary

Landmark study demonstrating that annotation software and transcript set selection dramatically affect variant interpretation outcomes. Compared 80 million variants using RefSeq vs Ensembl transcript sets (via ANNOVAR) and ANNOVAR vs VEP (using Ensembl). Key findings: only 44% agreement for putative loss-of-function variants between RefSeq and Ensembl; only 65% agreement between ANNOVAR and VEP for loss-of-function variants. Splicing variant annotations showed the greatest discordance.

## Relevance to cdot

Provides powerful motivation for cdot's multi-source approach (combining RefSeq and Ensembl) and for the importance of using specified transcript versions. The paper shows that inconsistent transcript selection is not a minor detail but has major consequences for research conclusions. cdot's comprehensive, versioned transcript coverage directly addresses this by (1) including both RefSeq and Ensembl, (2) preserving historical versions, and (3) storing which source each transcript came from. Also supports the ClinVar benchmark idea in cdot issue #5.
