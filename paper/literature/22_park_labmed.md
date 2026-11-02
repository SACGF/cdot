---
title: "Variations in Nomenclature of Clinical Variants between Annotation Tools"
authors: Kyoung-Jin Park, Jong-Ho Park
journal: Laboratory Medicine
year: 2022
volume: 53
pages: 242-245
doi: 10.1093/labmed/lmab074
pmid: 34612497
---

## Summary

Analysed over 218,000 clinical variants using ANNOVAR and SnpEff, finding ~85% HGVS nomenclature agreement between tools. SnpEff outperformed ANNOVAR for coding variants (99.3% vs 84.9% accuracy) and protein variants (94.3% vs 79.8%). Running both tools together achieved 99.5% accuracy. Concludes that substantial differences between tools warrant using multiple annotation tools in combination.

## Relevance to cdot

Demonstrates in a clinical context that HGVS naming inconsistency between tools is a real and ongoing problem. The 15% disagreement rate found here (using the same underlying data but different software) motivates the value of a standardised, tool-agnostic transcript data source like cdot. If annotation tools shared the same transcript data layer (via cdot), at least the coordinate mapping step would be consistent — disagreements would then only arise from consequence prediction differences, not transcript selection.
