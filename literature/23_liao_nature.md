---
title: "A draft human pangenome reference"
authors: Wen-Wei Liao, Mobin Asri, Jana Ebler, Daniel Doerr, Marina Haukness, et al. (Human Pangenome Reference Consortium)
journal: Nature
year: 2023
volume: 617
pages: 312-324
doi: 10.1038/s41586-023-05896-x
pmid: 37165242
---

## Summary

The Human Pangenome Reference Consortium presents 47 phased, diploid assemblies from genetically diverse individuals, covering >99% of expected sequence per genome at >99% accuracy. Adds 119 million base pairs of euchromatic polymorphic sequences and 1,115 gene duplications relative to GRCh38. Pangenome-based short-read analysis reduces small variant errors by 34% and increases structural variant detection by 104% compared to GRCh38 workflows.

## Relevance to cdot

The pangenome era creates new challenges for HGVS tooling — variants called on pangenome graphs need to be expressible in HGVS nomenclature, requiring transcript data for diverse reference sequences. cdot's current architecture supports multiple genome builds; extending to pangenome haplotypes is a future direction. The Ensembl 2023 paper notes pangenome haplotype annotations are being added to Ensembl, which would eventually appear in Ensembl GTFs that cdot consumes.
