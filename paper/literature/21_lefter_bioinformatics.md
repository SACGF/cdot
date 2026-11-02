---
title: "Mutalyzer 2: next generation HGVS nomenclature checker"
authors: Mihai Lefter, Jonathan K Vis, Martijn Vermaat, Johan T den Dunnen, Peter E M Taschner, Jeroen F J Laros
journal: Bioinformatics
year: 2021
volume: 37
pages: 2811-2817
doi: 10.1093/bioinformatics/btab051
pmid: 33538839
---

## Summary

Mutalyzer 2 is a complete rewrite of the original Mutalyzer HGVS nomenclature checker (see Wildeman 2008). Over five years the system processed over 133 million variant descriptions. Analysis found ~50% of submitted descriptions were syntactically/semantically correct, ~41% had errors, and ~7% could be automatically corrected. The new version is open-source (GNU AGPL) and available at mutalyzer.nl.

## Relevance to cdot

Demonstrates the scale of demand for HGVS tooling and the high error rate in real-world variant descriptions (~50% of 133M submissions incorrect). Mutalyzer focuses on validation and correction of HGVS syntax; cdot focuses on providing the underlying transcript data for coordinate conversion. These are complementary tools. The error rate finding also motivates cdot's #27 issue (HGVS string cleaning framework).
