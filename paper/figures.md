# Figures

![](paper/figures/figure1.svg)

**Figure 1. cdot architecture and data flow.** **(A)** Data generation (build time):
historical RefSeq GFF3 and Ensembl GTF releases, plus UTA alignments absent from any
annotation file, are parsed, normalised, and merged chronologically (newest wins) by a
Snakemake pipeline into gzip-compressed JSON, hosted on GitHub and also loaded into the
cdot_rest REST API (cdotlib.org). **(B)** Variant resolution (client): an incoming HGVS
string, possibly malformed, is repaired by `clean_hgvs()`, then parsed and mapped
between c. and g. coordinates by biocommons/hgvs, which draws exons and alignments from
one of cdot's data providers (in-memory `JSONDataProvider`, `RESTDataProvider`, or
`EnsemblTarkDataProvider`) and transcript sequence from the sequence layer (SeqRepo
and/or `FastaSeqFetcher`, tried in order via `ChainedSeqFetcher`). Opt-in helpers
rewrite the input string: `resolve_gene_hgvs()` maps a gene symbol to its MANE
transcript, and `get_best_transcript_version()` substitutes an adjacent transcript
version after a coordinate-safety check.
