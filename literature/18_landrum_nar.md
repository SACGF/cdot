---
title: "ClinVar: improving access to variant interpretations and supporting evidence"
authors: Melissa J Landrum, Jennifer M Lee, Mark Benson, Garth R Brown, Chen Chao, Shanmuga Chitipiralla, Baoshan Gu, Jennifer Hart, Douglas Hoffman, Wonhee Jang, Karen Karapetyan, Kenneth Katz, Chunlei Liu, Zenith Maddipatla, Adriana Malheiro, Kurt McDaniel, Michael Ovetsky, George Riley, George Zhou, J Bradley Holmes, Brandi L Kattman, Donna R Maglott
journal: Nucleic Acids Research
year: 2018
volume: 46(D1)
pages: D1062-D1067
doi: 10.1093/nar/gkx1153
pmid: 29165669
---

## Summary

Update to ClinVar, the freely available public archive of human genetic variant interpretations and their significance to disease. Enhancements include support for phenotypic submissions from clinical providers and patient registries, improved search indexing, and refined filtering. ClinVar aggregates clinical significance interpretations submitted by testing laboratories, research institutions, and expert panels.

## Relevance to cdot

ClinVar is the primary real-world dataset for benchmarking HGVS conversion tools. cdot issue #5 proposes using ClinVar HGVS entries to benchmark coverage and accuracy (resolve all ClinVar HGVS against the known genomic coordinates in the VCF). Most variants in ClinVar are described using HGVS transcript notation — they can only be converted to genomic coordinates using a tool like cdot + biocommons/hgvs. The ~50,000+ variants in ClinVar represent the practical demand that cdot serves.
