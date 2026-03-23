---
title: "A Python package for parsing, validating, mapping and formatting sequence variants using HGVS nomenclature"
authors: Reece K Hart, Rudolph Rico, Emily Hare, John Garcia, Jody Westbrook, Vincent A Fusaro
journal: Bioinformatics
year: 2015
volume: 31
pages: 268-70
doi: 10.1093/bioinformatics/btu630
pmid: 25273102
---

## Summary

The original publication of the biocommons `hgvs` Python library. Introduces an open-source package for parsing, manipulating, formatting, and validating biological sequence variants according to the HGVS specification. Designed for clinical diagnostics with high-throughput sequencing in mind. Released under Apache 2.0 licence. Prior to this work, HGVS handling tools were largely unavailable in a freely reusable form despite widespread adoption of HGVS nomenclature in clinical practice.

## Relevance to cdot

cdot is a data provider for this library. The `JSONDataProvider` and `RESTDataProvider` in cdot implement the `hgvs.dataproviders.interface.Interface` defined here. cdot's primary raison d'être is to replace UTA as the transcript data backend for this package, offering better coverage, no PostgreSQL dependency, and Ensembl support.
