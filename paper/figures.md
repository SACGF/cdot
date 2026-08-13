# Figures

![](paper/figures/figure1.svg)

**Figure 1. cdot architecture and data flow.** **(A)** Data generation (build
time). Historical RefSeq GFF3 and Ensembl GTF releases, plus UTA alignments
absent from any annotation file, are parsed, normalised, and merged in
chronological order (newest wins) by a Snakemake pipeline into gzip-compressed
JSON, hosted on GitHub and also loaded into the cdot_rest REST API
(cdotlib.org). **(B)** Variant resolution (client). An incoming HGVS string,
possibly malformed, is repaired by `clean_hgvs()`, then parsed and mapped
between c. and g. coordinates by biocommons/hgvs, which draws exons and
alignments from one of cdot's data providers (in-memory `JSONDataProvider`,
`RESTDataProvider`, or `EnsemblTarkDataProvider`) and transcript sequence from
the sequence layer (SeqRepo and/or `FastaSeqFetcher`, tried in order via
`ChainedSeqFetcher`). Opt-in helpers rewrite the input string:
`resolve_gene_hgvs()` maps a gene symbol to its MANE transcript, and
`get_best_transcript_version()` substitutes an adjacent transcript version
after a coordinate-safety check. Shapes distinguish node types (rounded
rectangle: process or code; folded corner: file; cylinder: database;
parallelogram: string input or output); blue nodes are cdot components, grey
nodes are external.
