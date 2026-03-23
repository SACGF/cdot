---
title: "GENCODE: reference annotation for the human and mouse genomes in 2023"
authors: Adam Frankish, Sílvia Carbonell-Sala, Mark Diekhans, Irwin Jungreis, Jane E Loveland, Jonathan M Mudge, Cristina Sisu, James C Wright, et al. (Flicek)
journal: Nucleic Acids Research
year: 2023
volume: 51(D1)
pages: D942-D949
doi: 10.1093/nar/gkac1071
pmid: 36420896
---

## Summary

Annual update to GENCODE, the high-quality reference annotation for human and mouse genomes, supported entirely by experimental evidence. 2023 updates include: characterisation of non-canonical open reading frames, LRGASP initiative evaluating long-read transcriptomic data for transcript modelling, enhanced convergence with RefSeq and UniProt on protein-coding gene annotation, expansion across the human pan-genome, and new regulatory feature annotation tools. Available through Ensembl, UCSC, and gencodegenes.org.

## Relevance to cdot

GENCODE is the annotation underpinning Ensembl's gene/transcript data — Ensembl is supposed to use GENCODE transcripts for human. cdot already downloads GENCODE HGNC gz files for HGNC ID/symbol mapping. Open issue #90 asks whether GENCODE GTFs should also be used directly as a transcript source. The convergence of GENCODE, RefSeq, and Ensembl annotation is directly relevant to cdot's merge strategy.
