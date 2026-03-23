---
title: "VariantValidator: Accurate validation, mapping, and formatting of sequence variation descriptions"
authors: Peter J Freeman, Reece K Hart, Liam J Gretton, Anthony J Brookes, Raymond Dalgleish
journal: Human Mutation
year: 2018
volume: 39
pages: 61-68
doi: 10.1002/humu.23348
pmid: 28967166
---

## Summary

VariantValidator is a web-based tool and Python library for validating, mapping, and formatting HGVS variant descriptions. It leverages the biocommons `hgvs` package and adds a user-friendly interface, error correction, and format conversion. Claims "a degree of accuracy that surpasses most competing solutions" for converting between variant description formats, with particular emphasis on genomic ↔ transcript coordinate conversion.

## Relevance to cdot

VariantValidator is a direct functional peer of the cdot-powered biocommons/hgvs pipeline — both solve the HGVS coordinate conversion problem. VariantValidator's approach is a monolithic service with its own data management; cdot's is a lightweight data layer that plugs into existing libraries. An important distinction: VariantValidator requires its own infrastructure to self-host, whereas cdot's JSON files can be loaded locally with no server. Both cite the 2016 HGVS standard and biocommons/hgvs.
