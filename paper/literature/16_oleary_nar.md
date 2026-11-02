---
title: "Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation"
authors: Nuala A O'Leary, Mathew W Wright, J Rodney Brister, Stacy Ciufo, Diana Haddad, Rich McVeigh, Bhanu Rajput, et al. (Murphy, Pruitt)
journal: Nucleic Acids Research
year: 2016
volume: 44(D1)
pages: D733-D745
doi: 10.1093/nar/gkv1189
pmid: 26553804
---

## Summary

Comprehensive description of the NCBI RefSeq database, covering its curated, non-redundant reference sequences for genomes, transcripts, and proteins across >55,000 organisms. Details computational and manual curation processes, institutional partnerships, and expansion into RNA-Seq integration and functional annotation. The database is the authoritative source for human transcript accessions (NM_, NR_, NP_ prefixes) that form the backbone of clinical HGVS nomenclature.

## Relevance to cdot

RefSeq GFF3 files (multiple annotation releases through RS_2025_08) are the primary source for RefSeq transcript data in cdot. cdot parses RefSeq GFF3 using its `GFF3Parser` class, extracting transcript/exon coordinates, CDS boundaries, CCDS IDs, HGNC IDs from Dbxref, and alignment gap CIGAR strings from cDNA_match features. The versioning scheme (NM_XXXXX.N) that makes historical HGVS resolution necessary is a RefSeq-defined system.
