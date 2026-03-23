---
title: "TransVar: a multilevel variant annotator for precision genomics"
authors: Wanding Zhou, Tenghui Chen, Zechen Chong, Mary A Rohrdanz, James M Melott, Chris Wakefield, Jia Zeng, John N Weinstein, Funda Meric-Bernstam, Gordon B Mills, Ken Chen
journal: Nature Methods
year: 2015
volume: 12
pages: 1002-1003
doi: 10.1038/nmeth.3622
pmid: 26513549
---

## Summary

TransVar is a command-line tool for converting between genomic, transcript, and protein-level variant representations (e.g., genomic VCF ↔ HGVS transcript notation ↔ protein change). Supports multiple annotation databases including RefSeq, Ensembl, GENCODE, and CCDS. Handles both simple substitutions and complex structural variants.

## Relevance to cdot

TransVar is a functional peer to the cdot + biocommons/hgvs pipeline for coordinate conversion, but with a different implementation approach (standalone command-line tool rather than a Python library with a pluggable data layer). Both solve the genomic ↔ HGVS conversion problem. TransVar bundles its own annotation rather than depending on an external data provider, which limits its transcript coverage to what is pre-packaged and makes historical version support harder. cdot's architecture of separating data from code is more flexible.
