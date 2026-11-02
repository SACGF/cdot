---
title: "The complete sequence of a human genome"
authors: Sergey Nurk, Sergey Koren, Arang Rhie, Mikko Rautiainen, Andrey V Bzikadze, Alla Mikheenko, Mitchell R Vollger, Nicolas Altemose, et al. (T2T Consortium)
journal: Science
year: 2022
volume: 376
pages: 44-53
doi: 10.1126/science.abj6987
pmid: 35357919
---

## Summary

Presents the first complete (telomere-to-telomere) sequence of a human genome: T2T-CHM13, a 3.055 billion-base pair assembly with gapless coverage of all chromosomes except Y. Fills the previously missing ~8% of the human genome, including all centromeric satellite arrays, segmental duplications, and the short arms of acrocentric chromosomes. Adds nearly 200 million base pairs of sequence with 1,956 gene predictions, 99 predicted protein-coding.

## Relevance to cdot

cdot is among the first HGVS tools to support T2T-CHM13v2.0 as a genome build (alongside GRCh37 and GRCh38). This is potentially a novel claim — no other HGVS coordinate conversion library has documented T2T support. The paper establishes why T2T matters: it corrects errors in prior references and unlocks previously inaccessible genomic regions to variant analysis. cdot's T2T support enables HGVS conversion for transcripts annotated on this assembly, using both Ensembl and RefSeq annotation sources.
