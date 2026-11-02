---
title: "The Ensembl Variant Effect Predictor"
authors: William McLaren, Laurent Gil, Sarah E Hunt, Harpreet Singh Riat, Graham R S Ritchie, Anja Thormann, Paul Flicek, Fiona Cunningham
journal: Genome Biology
year: 2016
volume: 17
pages: 122
doi: 10.1186/s13059-016-0974-4
pmid: 27268795
---

## Summary

Describes the Ensembl Variant Effect Predictor (VEP), a comprehensive toolset for annotating and prioritising genomic variants in coding and non-coding regions. Provides access to extensive genomic annotation through multiple interfaces (web, API, command line). Open source and fully reproducible. Supports both Ensembl/GENCODE and RefSeq gene sets for annotation. Widely used in research and clinical genomics pipelines.

## Relevance to cdot

VEP is the dominant variant annotation tool that cdot data might interact with. VEP uses its own internal annotation rather than cdot data for variant consequence prediction, but cdot and VEP use the same source data (RefSeq GFFs, Ensembl GTFs). VEP-compatible cdot releases are a planned feature (issue #63) — providing cdot JSON matched to specific VEP versions would let users cross-reference tools. The McCarthy 2014 paper showed that transcript set choice (RefSeq vs Ensembl, and ANNOVAR vs VEP) dramatically affects annotation outcomes, directly motivating cdot's multi-source approach.
