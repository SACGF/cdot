---
title: "APPRIS principal isoforms and MANE Select transcripts define reference splice variants"
authors: Fernando Pozo, José Manuel Rodriguez, Laura Martínez Gómez, Jesús Vázquez, Michael L Tress
journal: Bioinformatics
year: 2022
volume: 38(Suppl_2)
pages: ii89-ii94
doi: 10.1093/bioinformatics/btac473
pmid: 36124785
---

## Summary

Evaluates APPRIS principal isoforms and MANE Select transcripts as reference splice variants using proteomics data, genetic variation, and transcript information. Finds that most coding genes have a single main protein isoform and that APPRIS and MANE Select "coincide with the main proteomics isoform for over 98.2% of genes." Using the longest transcript as the reference is shown to be a poor strategy — longest-transcript-specific exons lack evolutionary conservation. Roughly 95% of genes have available MANE or APPRIS representatives.

## Relevance to cdot

Provides strong biological evidence supporting the use of MANE Select and APPRIS for canonical transcript selection — the same tags cdot already stores. Confirms that these approaches outperform simple heuristics (longest transcript, first alphabetically). Supports the design direction for cdot's canonical transcript selector (issue #36).
