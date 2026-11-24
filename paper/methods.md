# Methods

## Data sources and generation

cdot data collection works via a Snakemake pipeline downloading all of the official RefSeq
GFF3 [@OLeary2016] and Ensembl GTF [@Martin2023] releases for GRCh37, GRCh38, and T2T
(Figure 1A; Supplementary Tables S1 and S2). We also download UTA, which, as it computes
its own alignments, has a few accessions that were never in any of the official releases.
We also provide a cdot release which includes all of the historical releases RefSeq or
Ensembl have made for a particular genome build. The historical versions are loaded with
UTA first and then chronologically, so the official releases will override any custom
alignments performed by UTA, and older versions will be overwritten by newer ones which
contain new information such as gaps.

## JSON format

The JSON format has one entry per transcript version with build-specific alignment data in
a `genome_builds` dictionary keyed by build. The layout of the JSON data is the same
whether `genome_builds` contains one or multiple genome builds. The exon data contains a
`gap` field that records the indels found during alignment, and can be applied for
handling gaps during conversion. We copy tags from the original GFF or GTF, such as MANE
Select, MANE Plus Clinical, RefSeq Select, and Ensembl canonical. We include a top-level
`schema_version`, which can be used to check what information it contains and whether it
can be parsed by older clients (Supplementary Table S3).

## String cleaning

For cleaning HGVS strings, we provide a function `fix_hgvs()`. You pass in an HGVS string
and get back the repaired string plus a list of the fixes applied. Some examples include
whitespace, quoting, double dots (`NM_000059..4`), uppercase, lowercase (`c.123delg`), and
the swapped gene and transcript (`BRCA2(NM_000059.4)` → `NM_000059.4(BRCA2)`)
(Supplementary Table S7). If you pass a data provider to `fix_hgvs()`, it can also fix
some transcript errors, such as a dropped accession prefix (`000059.4:c.68del` →
`NM_000059.4:c.68del`). It only does this when exactly one candidate is possible. If you
also pass a genome build, it resolves a gene-only string, such as `BRCA2:c.68_69del`, by
looking up the MANE transcript [@Morales2022; @Wright2023]. Every change is recorded as an
`HGVSFix`, which has severity, code, and before/after, so a caller can audit what was
changed or present them to the user. What it can't fix is flagged as an error-level fix.
We also have `clean_hgvs()`, which is called by `fix_hgvs()`, but it is a pure string
conversion function and doesn't load any data.

## Version substitution

When the requested transcript version isn't available, `fix_hgvs()` can substitute in
another version. This is off by default because it can be wrong. You can pass in a
strategy, `version_fallback`, which you can set to try higher versions then lower,
closest, or the latest. It's reported back as an `HGVSFix`, so the user can see what
happened. We have a check that each transcript version has an identical CDS structure: the
CDS length plus the length of coding exons. This is the same on every genome build, so the
two versions can be compared even when one of them is not aligned or not available for the
genome build you request. It's safe to swap only if the structures match. The alignment
gaps must also match, and for a UTR variant, the UTR lengths must be the same
(Supplementary Methods). If the version substitution doesn't match, by default we just
leave the string unchanged and it returns an error `HGVSFix`. You can override this to
change the string and return a warning that the coordinates may have moved.

## Access and client libraries

cdot implements a biocommons HGVS data provider, so it can be swapped in to load the cdot
transcript versions (Figure 1B). We provide three providers: `JSONDataProvider`, which
loads the gzipped JSON into memory; `RESTDataProvider`, which makes network calls to
cdotlib.org (it also has a prefetch, which batches things into more efficient calls); and
`EnsemblTarkDataProvider`, which reads from Ensembl's official transcript archive, TARK.
biocommons provides SeqRepo [@Hart2020], but we also make a seq fetcher
(`FastaSeqFetcher`) that can load a local FASTA file, which pulls out the exon sequences
and pastes them together to make a transcript. This works perfectly for Ensembl
transcripts, which exactly match the genome, but might have some troubles with RefSeq,
where the transcript reference and the genome reference can have gaps and mismatches.
We've also made a chained seq fetcher (`ChainedSeqFetcher`), which can try several data
sources and then fall back on others. You can use this to take the SeqRepo sequences when
available and then fall back on the other sources.

## Evaluation

Results from this paper that involve public data can be reproduced from a pipeline in this
repository (`paper/scripts`, described in the Supplementary Methods). For HGVS resolution,
we use ClinVar version {{ clinvar_vcf._prov.tools.clinvar }} and the VCV XML as historical
data, which contain the original c.HGVS with the transcript versions provided by labs at
submission time. We score converting c.HGVS to genomic coordinates against ClinVar's own
VCF coordinate (CHROM/POS/REF/ALT). For a performance benchmark, we run
{{ benchmark.n_pairs | int }} ClinVar pairs through the cdot local JSON, REST and the UTA
local. For the UTA remote, we run only {{ benchmark.uta_remote_n | int }} pairs because it
is slow. The HGVS engine and sequence source are held consistent between runs. We use two
data sets to test cleaning functionality: our own production search bar collection, which
is {{ cleaning_corpus.corpus_n | commas }} queries and private data, and the
{{ cleaning.inject_sample_size | commas }} clean public ClinVar strings with synthetically
injected errors, {{ lovd_comparison.n_cases | commas }} cases in total. The synthetic data
is also run on third-party cleaners LOVD [@LovdHgvsChecker], VariantValidator
[@Freeman2018], and Mutalyzer [@Lefter2021]. We check transcript version substitution
safety by comparing consecutive versions of {{ version_stability.sample_n | commas }}
transcripts each from both RefSeq and Ensembl GRCh38. The invalid HGVS from production
searches, used as examples of real-world data cleaning, are private because they contain
clinical information. The HGVS sent from clinical labs into Shariant [@Tudini2022] are
also used as a benchmark. We cannot provide this as it is private clinical data. We have
converted this to aggregate counts which are included in the repository.
