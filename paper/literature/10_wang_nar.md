---
title: "ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data"
authors: Kai Wang, Mingyao Li, Hakon Hakonarson
journal: Nucleic Acids Research
year: 2010
volume: 38
pages: e164
doi: 10.1093/nar/gkq603
pmid: 20601685
---

## Summary

Introduces ANNOVAR, a tool for annotating SNVs and indels from high-throughput sequencing data. Supports gene-based, region-based, and filter-based annotations. Can annotate functional consequences, cytogenetic bands, conservation scores, and cross-reference to population variation databases. Processes 4.7 million variants in 4-15 minutes on a desktop. Successfully applied to identify the causal gene for Miller syndrome.

## Relevance to cdot

ANNOVAR is one of the most widely used variant annotation tools and a key comparator in the McCarthy 2014 study that showed transcript set choice dramatically affects annotations. ANNOVAR does not use HGVS as its primary format and bundles its own transcript databases. cdot's approach — providing a shared, versioned transcript data layer with HGVS-aware coordinate conversion — is philosophically different: rather than a monolithic annotation tool, cdot is infrastructure that enables consistent HGVS conversion across any tool.
