---
title: "CSN and CAVA: variant annotation tools for rapid, robust next-generation sequencing analysis in the clinical setting"
authors: Márton Münz, Elise Ruark, Anthony Renwick, Emma Ramsay, Matthew Clarke, Shazia Mahamdallie, Victoria Cloke, Sheila Seal, Ann Strydom, Gerton Lunter, Nazneen Rahman
journal: Genome Medicine
year: 2015
volume: 7
pages: 76
doi: 10.1186/s13073-015-0195-6
pmid: 26315209
---

## Summary

Introduces Clinical Sequencing Nomenclature (CSN) — a deterministic HGVS-aligned variant naming scheme — and CAVA, a fast command-line annotation tool. CAVA processed over 10 million ExAC variants in 13.44 hours and achieved 100% concordance with Sanger data for BRCA1/BRCA2 mutations, versus most other tools correctly annotating only 32% of BRCA2 variants. Addresses clinical annotation challenges including indel handling and strand orientation issues.

## Relevance to cdot

Highlights a persistent problem: most annotation tools have correctness issues for indels and strand-specific HGVS naming. CAVA's superior performance is partly due to stricter adherence to HGVS rules. cdot provides the transcript data layer; the correctness of HGVS conversion still depends on the consuming library. The 32% correct rate for BRCA2 with other tools is a dramatic illustration of the problem that accurate transcript data + correct HGVS libraries (biocommons/hgvs) are meant to solve.
