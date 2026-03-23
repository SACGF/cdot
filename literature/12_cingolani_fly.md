---
title: "A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3"
authors: Pablo Cingolani, Adrian Platts, Le Lily Wang, Melissa Coon, Tung Nguyen, Luan Wang, Susan J Land, Xiangyi Lu, Douglas M Ruden
journal: Fly (Austin)
year: 2012
volume: 6
pages: 80-92
doi: 10.4161/fly.19695
pmid: 22728672
---

## Summary

Introduces SnpEff, a tool for rapidly annotating and predicting the effects of genetic variants in whole-genome sequencing data. Annotates based on genomic location and predicts coding effects (synonymous/non-synonymous changes, frameshift, stop codon alterations, splice site effects). Demonstrated on Drosophila melanogaster; identified ~356,660 SNPs with ~15,842 synonymous and ~4,467 non-synonymous variants.

## Relevance to cdot

SnpEff is one of the two dominant open-source variant annotation tools (alongside VEP), and like VEP, uses its own bundled transcript data rather than cdot. The Park 2022 study showing 85% agreement between SnpEff and ANNOVAR highlights that different tools produce different HGVS annotations from the same variants. cdot's role as a shared transcript data layer — if adopted widely — could reduce this inconsistency by ensuring tools use the same coordinate data.
