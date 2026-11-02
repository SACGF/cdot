## [unreleased]

### Added

- #148 - HGVS cleaning: a gene symbol written before the transcript accession, separated by whitespace, is moved into the standard `ACCESSION(SYMBOL)` form (`APC NM_001407446.1:c.-40G>A` → `NM_001407446.1(APC):c.-40G>A`). Previously whitespace stripping glued them together (`APCNM_001407446.1:c.-40G>A`). New op `HGVSCleanOp.FIX_SYMBOL_BEFORE_ACCESSION`, reported as a WARNING with the new fix code `WRAPPED_GENE_SYMBOL`
- #148 - HGVS cleaning: a chromosome name before its own chromosomal accession is dropped (`chr13 NC_000013.11:g.32316528del` → `NC_000013.11:g.32316528del`), for GRCh37, GRCh38 and T2T-CHM13v2.0 accessions. Only exact `chr1`-`chr22`/`chrX`/`chrY`/`chrM`/`chrMT` names are recognised (gene symbols such as `CHRNA7` are not affected), and the name is only dropped when it matches the accession's chromosome, so `chr7 NC_000013.11:…` is left alone. New op `HGVSCleanOp.STRIP_CONTIG_NAME`, reported as a WARNING with the new fix code `STRIPPED_CONTIG_NAME`
- #145 - `models.Transcript.fake`: `True` for the RefSeq mitochondrial fake transcripts (data schema >= 0.2.37), `None` for real ones
- #166 - `clean_hgvs` warns with the new `HGVSCleanOpWarning` (also exported from `cdot.hgvs`) when the `ops` set selects an op without an op it relies on, eg `LOWERCASE_MUTATION_TYPE` without `RECONSTRUCT_STRUCTURE`, which inserts the colon it looks after. The ops still run as before, and cleaning output and fixes are unchanged, except that `FIX_GENE_TRANSCRIPT` run without `STRIP_WHITESPACE` no longer drops a newline after the gene parentheses (`NM_000059.4(BRCA2)\n:c.68del`)
- #181 - `EnsemblTarkDataProvider` and `EnsemblTarkSeqFetcher` are exported from `cdot.hgvs.dataproviders` alongside `JSONDataProvider` and `RESTDataProvider`, and the README and biocommons examples now mention the Tark provider
- #181 - Without pysam, `from cdot.hgvs.dataproviders import FastaSeqFetcher` still works and constructing it raises `ImportError("FastaSeqFetcher needs pysam: pip install cdot[fasta]")` instead of the name being missing (`NameError`/`ImportError` with no hint)
- #181 - `__all__` on `cdot.hgvs`, `cdot.hgvs.dataproviders` and its `json_data_provider`, `seqfetcher` and `fasta_seqfetcher` modules, so `from cdot.hgvs.dataproviders import *` and `dir()` give the public names only (not `abc`, `requests`, `IntervalTree` and so on). Every name imported from these modules by the README, docs, tests and cdot_rest is still importable. `JSONDataProvider` and `RESTDataProvider` have docstrings describing every constructor argument, including the biocommons `mode` and `cache` cache-mode arguments, and type-hinted constructors
- #171 - `get_transcript(tx_ac)` on the biocommons HGVS data providers (`JSONDataProvider`, `RESTDataProvider`): the public way to read cdot's transcript record for a versioned accession (`None` if not held). `fix_c_to_p` now reads the record through it instead of the provider's private `_get_transcript`, so a custom provider only needs to implement `get_transcript` to give `fix_c_to_p` its translation data

### Changed

- #170 - `LocalDataProvider.get_tx_for_region` returns the record's versioned accession (`transcript_data["id"]`, eg `NM_000059.4`) as `tx_ac`, as `get_tx_for_gene` and `get_tx_ac_tags_for_gene` already do, instead of the id the provider stores the transcript under. `JSONDataProvider` is unaffected (its ids are already versioned); a provider keyed by versionless id (eg cdot_rest's Redis provider) now returns `NM_000059.4` rather than `NM_000059`
- #170 - A contig outside the provider's assemblies (eg `chr7`) raises `HGVSDataNotAvailableError` from `get_tx_exons` instead of `ValueError`, matching `get_tx_info`, so biocommons hgvs reports it as missing data like any other
- #171 - `resolve_gene_hgvs`, `rank_transcripts_for_gene` and `fix_hgvs`: the consortium argument is now `annotation_consortium`, the name the data pipeline and REST API already use (the old name `prefer_consortium` said preference but it has always been a hard filter, and still is). Same position in the signature and same values (`Consortium.REFSEQ` default, `ENSEMBL`, or `None` for either), so positional callers are unaffected; the "no transcript in this consortium" ERROR message now says to pass `annotation_consortium=...`
- #171 - `fallback_to_longest` now compares the candidates' lengths explicitly (spliced transcript length, the sum of the exon lengths from the provider's `get_tx_identity_info`) instead of trusting the data provider to return transcripts longest-first. With cdot's providers the result is the same transcript; a provider that can't give a length (no `get_tx_identity_info`) gets the first transcript it returned and the `RESOLVED_GENE_TO_TRANSCRIPT` WARNING now says "first transcript returned by the data provider" rather than "longest transcript" in that case
- #171 - `fix_c_to_p` returns the input `var_c` unchanged (not `None`) as the first element when there is an ERROR-level fix, matching `clean_hgvs`/`fix_hgvs` which return their input string on an error; check `fixes` for the ERROR. The return type is annotated `tuple[SequenceVariant, list[HGVSFix]]`
- #171 - `fix_hgvs` with a `data_provider` but no `genome_build` records a WARNING with the new fix code `SKIPPED_GENE_RESOLUTION` when the input is gene-only (eg `BRCA2:c.36del`), as gene resolution needs a build, instead of leaving the string alone with no fix. An input that already names a transcript gets no such fix, and nothing is recorded without a data provider
- #171 - `tag_priority` on `resolve_gene_hgvs`, `rank_transcripts_for_gene` and `fix_hgvs` defaults to `None`, meaning a copy of `DEFAULT_TAG_PRIORITY`, rather than the shared list itself
- #171 - The `USED_ADJACENT_VERSION_COORD_UNVERIFIED` message reads `...; coordinate may have moved (...)` (a semicolon in place of a dash)
- #180 - Removed the `hgvs` extra. `hgvs` is a core dependency, so the extra did nothing; `pip install cdot[hgvs]` still installs cdot (pip only warns about the unknown extra)

### Changed

- #168 - `JSONDataProvider.get_tx_versions()` looks the versions up in a per-accession index built lazily on first use, instead of scanning every loaded transcript on each call. The version fallback in `resolve_transcript_version()` (and `fix_hgvs()` with it) no longer pays a full scan for every versionless or unavailable transcript version. Results are unchanged
- #179 - Package metadata lists the supported Python versions (3.10 to 3.14) in its classifiers

### Deprecated

- #171 - `prefer_consortium` on `resolve_gene_hgvs`, `rank_transcripts_for_gene` and `fix_hgvs` is a deprecated alias of `annotation_consortium` and emits a `DeprecationWarning`; passing both is a `TypeError`. It will be removed in a future release
- #186 - `AbstractJSONDataProvider.is_version_substitution_safe()` is deprecated and emits a `DeprecationWarning`. The check is now the pure `cdot.hgvs.version_safety.is_version_substitution_safe(transcript, substitute_transcript, genome_build=None, cited_position=None)` on the two versions' transcript records, with `genomic_bracket_safe(transcript_below, transcript_above, genome_build)` as the fallback for a version held in no build. `is_version_substitution_safe_by_accession(data_provider, transcript_accession, substitute_transcript_accession, genome_build=None, cited_position=None)` is the same check on two versioned accessions, fetching the records from a data provider; the deprecated provider method delegates to it and will be removed in a future release

### Fixed

- #158 - Transcript version fallback safety check (`is_version_substitution_safe` and the `cdot.hgvs.version_safety` helpers): a malformed record (exon transcript coordinates `None`) now gives "cannot verify" instead of raising `TypeError`. Non-coding transcripts (`NR_`) with both versions loaded are now checked on their exon structure (per-exon lengths, alignment gaps, genomic placement) instead of being refused as "requested version absent", and a coding/non-coding pair is refused. An uncertain range starting in the 5'UTR (eg `c.(-30_5)del`) is now treated as a 5'UTR position, and a strand change is reported as such
- #156 - `cdot.data_release.get_latest_combo_file_urls` now finds the all-builds files (`cdot-<version>.all-builds-refseq-grch37_grch38_t2t-chm13v2.0.json.gz`), which it previously never returned. Ask for them with the genome build `"all-builds"` (`cdot.data_release.ALL_BUILDS`), the full `"grch37_grch38_t2t-chm13v2.0"`, or an old style combination such as `"grch37_grch38"`. Per-build requests (eg `"GRCh38"`) and older release names still work
- #154 - `RESTDataProvider` only caches a 404 as a missing transcript/gene. Other HTTP errors (eg 429 or 503) now raise `requests.HTTPError` instead of returning `None`, and are not cached, so a transient server error no longer makes a transcript permanently missing for the provider's lifetime. `prefetch()`'s concurrent fallback leaves failed accessions uncached so a later lookup retries them
- #154 - `RESTDataProvider` uses a `requests.Session` (connection reuse) and has a `timeout` argument (default 30 seconds) so a hung connection no longer blocks forever. The `session` argument lets you pass your own. `EnsemblTarkDataProvider`, `RESTPyHGVSTranscriptFactory` and the data release lookup also set a 30 second timeout
- #154 - A response with no `Content-Type` header raised `TypeError` in `RESTDataProvider`, `EnsemblTarkDataProvider` and `RESTPyHGVSTranscriptFactory`. It now gives the usual "Non-json response" `ValueError`
- #155 - `RESTDataProvider` rejected T2T-CHM13v2.0 contigs with `ValueError("Contig ... not supported")`. It now supports GRCh37, GRCh38 and T2T-CHM13v2.0 by default, set with the new `assemblies` argument
- #169 - `JSONDataProvider` and `JSONPyHGVSTranscriptFactory` close the JSON files they open (they leaked the handle, giving `ResourceWarning: unclosed file` under `python -W error` or `-X dev`). Both now load through `cdot.models`, which has a new `open_json(file_or_filename)` context manager for the gzip-or-plain open. `FastaSeqFetcher` opens each FASTA file once instead of twice

### Removed

- #103 - The mouse (Mus musculus GRCm38/GRCm39) RefSeq bash scripts under `generate_transcript_data/Mus_musculus/`. They have been dead since the move to Snakemake (#92) removed the wrapper and gene_info download they depended on, and mouse data was never part of a data release. cdot does not currently support non-human species: `bioutils` has no mouse assemblies for contig name mapping, so the data provider can't load a mouse file. #103 tracks what adding it would take

### Fixed

- #150 - HGVS cleaning: structure reconstruction no longer turns an input with a missing colon into a wrong but valid-looking string with only a WARNING (eg `NM_000059.4 123A>G` became `NM_000059.4123(A>):g.` and `NM_000059.4(BRCA2) 100G>A` became `NM_000059.4(BR):c.A2)100G>A`). Inputs that can't be repaired unambiguously are now left alone and reported as a `NO_COLON` ERROR, and a kind letter glued to the version is repaired correctly (`NM_005629.3c1403A>C` → `NM_005629.3:c.1403A>C`)
- #151 - HGVS cleaning: fixed ops that ran before the input they needed existed. Wrapping quotes are stripped first (`"BRCA2 NM_000059.4:c.68del"` → `NM_000059.4(BRCA2):c.68del`), a dropped `N` or underscore in the accession prefix is restored before the separator and empty-version fixes (`NM000059.4;c.68del` → `NM_000059.4:c.68del`, `M_001754.:c.749G>A` → `NM_001754:c.749G>A`), an uppercase mutation type is lowercased after the colon is inserted (`NM_000059.4 c.68_69DELINSAT` → `NM_000059.4:c.68_69delinsAT`), and a protein annotation after a tab is stripped (`NM_000059.4:c.68del\tp.(Glu23fs)`). `HGVSCleanOp` members are reordered to match the new pipeline order
- #152 - HGVS cleaning handles every accession prefix and kind letter. Lowercase genomic and protein accessions are uppercased (`ng_012337.3:g.5000A>G` → `NG_012337.3:g.5000A>G`, likewise `np_`/`nc_`) instead of being reported as `MISSING_REFERENCE_SEQUENCE`, `m.` and `r.` kinds are repaired like the others (`NC_012920.1 m.3243A>G` → `NC_012920.1:m.3243A>G`), `r.` bases are no longer uppercased (they are lowercase in HGVS), and a protein allele after a space is no longer dropped as a trailing protein annotation (`NP_000050.2 p.Arg1Gly` → `NP_000050.2:p.Arg1Gly`, was `NP_000050.2`)
- #165 - HGVS cleaning: a base count after a ranged del/dup is checked against the range before it is dropped. When both positions are plain integers (`c.`/`g.`/`n.`/`m.`, no intronic offset or UTR position) and the count agrees, it is dropped as before (`NM_000059.4:c.123_127del5` → `NM_000059.4:c.123_127del`); when it disagrees the string is left alone and reported as an ERROR with the new fix code `DEL_DUP_COUNT_MISMATCH` (`NM_000059.4:c.123_127del3`, a 5 bp range, was silently turned into `NM_000059.4:c.123_127del`). A range whose length needs the transcript (`NM_000059.4:c.100+5_200-3del5`) still has its count dropped with a WARNING
- #165 - HGVS cleaning: an unbalanced bracket is only stripped where it cannot change the meaning: around the reference (`(NM_004991.4:c.2577+4A>T`, `NM_000059.4(BRCA2:c.68del` → `NM_000059.4(BRCA2):c.68del`) or as a wrapping bracket that closes nothing (`NM_000059.4:c.68del)`). An unbalanced bracket inside the allele marks an uncertain position, so it is left alone and reported as an ERROR with the new fix code `UNBALANCED_BRACKETS` (`NM_000059.4:c.(100_200del` was silently turned into the certain `NM_000059.4:c.100_200del`). Balanced brackets are untouched as before
- #165 - HGVS cleaning: a genome build name before the reference (`hg19`, `hg38`, `GRCh37`, `GRCh38` with an optional patch suffix, `T2T`, `CHM13`, `T2T-CHM13v2.0`, any case) is dropped instead of being taken for a gene symbol (`hg19 NM_000059.4:c.68del` → `NM_000059.4:c.68del`, was `NM_000059.4(hg19):c.68del`; `GRCh38 BRCA2 NM_000059.4:c.68del` → `NM_000059.4(BRCA2):c.68del`, was `GRCh38BRCA2NM_000059.4:c.68del`). New op `HGVSCleanOp.STRIP_BUILD_NAME`, reported as a WARNING with the new fix code `STRIPPED_BUILD_NAME`; a `GRCh37`/`GRCh38` name directly before an accession (`GRCh38.p2 NM_000059.4:c.316+5G>A`) now reports this code instead of `STRIPPED_LEADING_JUNK`

## [0.2.33] - 2026-10-06

### Added

- #115 - `FastaSeqFetcher` (and `ExonsFromGenomeFastaSeqFetcher`) warn when building a transcript sequence that the data says differs from the genome (the build's `genome_mismatch`, RefSeq only), as that isn't the real transcript sequence. Set with the new `genome_mismatch` constructor argument: `'warn'` (default, a `GenomeMismatchWarning` once per transcript), `'raise'` (`GenomeMismatchError`, a subclass of `HGVSDataNotAvailableError`, so `ChainedSeqFetcher` falls through to the next seqfetcher) or `'off'`. `FastaSeqFetcher` now also passes its `cache` argument on, which it used to ignore
- #115 - `models.GenomeMismatch`: `transcript`, `protein` and `alignment` (data schema >= 0.2.36)

### Changed

- #113 - The biocommons HGVS data providers now use the `GRCh38.p14` assembly map for `GRCh38` (as `GRCh37` already used `GRCh37.p13`), so patch scaffold contigs (`NW_`) are recognised
- #113 - `FastaSeqFetcher` fetches `NW_` and `NT_` contigs from the genome FASTA (previously only `NC_`), so variants on patch scaffolds can be validated against a local FASTA

## [0.2.32] - 2026-09-27

### Added

- #131, #76 - `models`: `Transcript.translation` (new `Translation` struct), `GenomeBuild.genome_mismatch` (new `GenomeMismatch` struct) and `GenomeBuild.warnings`, all optional (data schema >= 0.2.35)
- #131, #76 - `cdot.hgvs.fix_c_to_p(variant_mapper, var_c)`: wraps a biocommons HGVS `VariantMapper`/`AssemblyMapper` `c_to_p()` and uses the translation data to return `(var_p, fixes)`, with problems as `HGVSFix` as in `fix_hgvs`. It picks the genetic code (mitochondrial `transl_table` 2, and the selenocysteine code for selenoproteins, which biocommons otherwise only detects when it can fetch the protein sequence, so eg with `FastaSeqFetcher` it read Sec as a stop and gave `p.?`), reports a new UGA in a selenoprotein as a stop rather than Sec (`NEW_UGA_READ_AS_STOP`), refuses transcripts with a ribosomal frameshift (`RIBOSOMAL_SLIPPAGE_UNSUPPORTED`), and checks translation exception codons against the transcript sequence used, so a genome derived sequence that differs from RefSeq (eg ACTN3 R577X on GRCh37) is an ERROR (`REFERENCE_CODON_MISMATCH`) instead of a silent `p.?`. RefSeq exceptions it can't apply and data conversion warnings are reported as WARNINGs

## [0.2.31] - 2026-08-19

### Added

- #112 - HGVS cleaning: `resolve_missing_accession_prefix(hgvs_string, data_provider, genome_build=None)` restores a RefSeq prefix the user dropped entirely, eg `000059.4:c.68del` → `NM_000059.4:c.68del`. The kind letter picks the candidate prefixes (`c.` → `NM_`/`XM_`, `n.` → `NR_`/`XR_`; a 7-8 digit field is zero-padded to the 9-digit form) and each candidate is checked against the data provider (via `get_tx_versions`), so the fix is applied only when exactly one candidate accession exists there - with zero or several matches the string is unchanged. Reported as a WARNING `HGVSFix` with the new code `ADDED_ACCESSION_PREFIX`, never silent. `fix_hgvs` runs it automatically whenever a `data_provider` is supplied (before gene resolution); without a provider the string is left alone
- #112 - HGVS cleaning: new string-level ops driven by the search-log corpus - a semicolon or underscore in place of the accession:allele colon (`NM_000059.4;c.68del` / `NM_001754.5_c.749G>A` → `…:c.…`, reported under `FIXED_SEPARATOR_TYPO`); a space in place of the reference:kind colon after a bare gene symbol (`BRCA2 c.68del` → `BRCA2:c.68del`, new op `HGVSCleanOp.FIX_SPACE_BEFORE_KIND`, which also repairs symbols that resemble an accession, eg `NR2E3 c.349G>A`); a dangling dot left by a dropped transcript version (`NM_001754.:c.749G>A` → `NM_001754:c.749G>A`, new op `HGVSCleanOp.DROP_EMPTY_VERSION` / fix code `DROPPED_EMPTY_VERSION`); and leading-junk stripping now also covers junk separated from the accession by whitespace or a stray dash (`": NM_000059.4(BRCA2):c.68del"`)
- #129 - Batch transcript-retrieval hook on the biocommons HGVS data providers: `_get_transcripts(tx_acs)` returns `{tx_ac: transcript_data}` and is now the single point every method that walks a list of transcripts goes through (`get_tx_for_gene`, `get_tx_ac_tags_for_gene`, `get_tx_for_region`, `_get_tags_by_tx_ac`). A provider backed by a database or a bulk REST call can override just this one method (one query instead of N) and speed all of them up. The default implementation loops `_get_transcript` per accession, which is what an in-memory provider such as `JSONDataProvider` already wants, so existing subclasses are unaffected; duplicate accessions collapse and the result follows the order of `tx_acs`. The batch tag hook `_get_tags_by_tx_ac` added in 0.2.28 keeps its signature and now retrieves through `_get_transcripts` too, so overriding `_get_transcripts` alone is enough to remove the N+1 from every one of these methods

### Changed

- #5 - Performance: the biocommons HGVS data providers now keep a small per-instance LRU cache (10,000 entries) of `get_tx_exons()` results, so repeated lookups of the same transcript/contig no longer rebuild the per-exon dict list (the per-variant hot path when resolving a stream of variants). Applies to `JSONDataProvider` and `RESTDataProvider`; results are cached only after a successful build, so lazily fetched providers still pick up transcripts that appear later (eg after a REST `prefetch()`)

## [0.2.30] - 2026-06-25

### Fixed

- Packaging: the wheel/sdist now ship only the `cdot` package. Previously setuptools namespace-package discovery swept in repo directories (`paper/`, `docs/`, `literature/`, `output/`) as top-level packages, so `pip install cdot` dumped `paper/scripts/*.py` and other non-library files into the environment. The package list is now an explicit `cdot*` allowlist

## [0.2.29] - 2026-06-25

### Added

- #120 - Coordinate-safety check for the opt-in transcript-version fallback: before substituting an adjacent version, cdot tests whether a coding `c.` position still maps to the same genomic coordinate and emits a richer `HGVSFix` (`USED_ADJACENT_VERSION_COORD_SAFE` when verified safe). Adds `is_version_substitution_safe(...)` on the data providers and `intrinsic_cds_structure()`/`cds_alignment_gaps()`/`utr_lengths()` in `cdot.hgvs.version_safety`.
- #28 - `get_tx_versions(accession, genome_build=None)` gained an optional `genome_build` argument that restricts the returned versions to those placeable in that build (the version-fallback candidate pool), so the fallback no longer offers a version that can't be mapped in the target build
- #112 - HGVS cleaning: new op that drops a genomic accession (`NC_`/`NG_`/`NW_`) wrongly placed in the gene-symbol parenthetical slot, eg `NM_000059.4(NC_000013.11):c.68del` → `NM_000059.4:c.68del`, so the string parses. Selectable as `HGVSCleanOp.DROP_GENOMIC_REF_IN_PARENS`; reported as an `HGVSFix` (`DROPPED_GENOMIC_REF_IN_PARENS`)
- #119 - HGVS cleaning now repairs the same formatting mistakes on all RefSeq accession types, not just `NM_`/`NC_`. The missing-kind step picks the kind from the reference-sequence type (mRNA `NM_`/`XM_` and Ensembl `ENST` → `c.`, non-coding RNA `NR_`/`XR_` → `n.`, genomic `NC_`/`NG_`/`NW_`/`NT_`/`AC_` → `g.`) and also handles 5'/3' UTR coordinates (eg `NM_000350.2:-52delC` → `:c.-52delC`); the missing-underscore step covers every RefSeq prefix (eg `NR003051.3` → `NR_003051.3`); the dropped-leading-`N` step covers `NM`/`NC`/`NR`/`NG`/`NT`/`NW`/`NS`/`NZ`/`NP` (eg `G_012337.3` → `NG_012337.3`); and the missing-reference-sequence validation flags a bare leading `n`/`m`/`p`/`r` kind as well as `c`/`g`

### Changed

- #119 - HGVS cleaning: renamed the colon/dot/underscore/kind collapse fix codes and ops to reflect that they collapse a run of 2+ characters, not only a pair - `HGVSFixCode.FIXED_DOUBLE_COLON`/`FIXED_DOUBLE_DOT`/`FIXED_DOUBLE_UNDERSCORE`/`FIXED_DOUBLE_KIND` are now `FIXED_MULTIPLE_COLON`/`FIXED_MULTIPLE_DOT`/`FIXED_MULTIPLE_UNDERSCORE`/`FIXED_MULTIPLE_KIND`, and the matching `HGVSCleanOp` members `FIX_DOUBLE_*` are now `FIX_MULTIPLE_*`. The colon, underscore and kind ops now also collapse a run of 3+ (eg `NM_000059.4:::c.123` → `NM_000059.4:c.123`, `NM___000059.4` → `NM_000059.4`, `c.c.c.123` → `c.123`), matching the dot op
- #28, #120 - Behaviour change to the opt-in version fallback: a substitution that is not structure-verified coordinate-safe is now refused by default (ERROR `HGVSFix` `REFUSED_UNSAFE_VERSION`, string left unchanged) rather than applied silently. Pass `on_unsafe_version=UnsafeVersionPolicy.SUBSTITUTE` to substitute anyway. Still opt-in (off unless `version_fallback=` is set). See #120 for the full description

### Fixed

- #119 - HGVS cleaning correctness fixes: reconstruction no longer mangles gene symbols that look like an accession (eg `NR3C1:c.1411A>G`), a kind letter glued to the version with a stray colon (eg `NM_000059.4c:123del`) now reconstructs correctly, and several fixes that were reported under the wrong `HGVSFix` code now use the right one. See #119 for the full list

## [0.2.28] - 2026-06-18

### Added

- #114 - Public HGVS ranking primitives + provider hooks (to let consumers stop importing private symbols / re-implementing ranking): `rank_transcript_versions` (the full version ordering behind `get_best_transcript_version`) and `rank_transcripts_for_gene` (the gene→transcript lookup/uppercase-retry/consortium-filter/tag-ranking core split out of `resolve_gene_hgvs`, returning the full ranked candidate list); `consortium_of` (promoted from `_consortium_of`); `HGVSFix.__str__` plus `messages`/`warning_messages`/`error_messages` helpers; an overridable batch tag hook `_get_tags_by_tx_ac(tx_acs, genome_build)` on data providers (default loops the per-transcript hook) so subclasses can answer in one query instead of N+1; and `cdot.hgvs` now also re-exports `Consortium`, `DEFAULT_CONSORTIUM`, `HGVSCleanOp`, and `ALL_CLEAN_OPS`
- #112 - HGVS cleaning: new ops driven by the search-log corpus - colon used in place of the kind dot (`:c:123` → `:c.123`), parens wrapping the whole accession (`(NM_000059.4):c.…` → `NM_000059.4:c.…`), and a gene-first wrapper missing the colon before the kind (`BRCA1(NM_…)c.…` → `NM_…(BRCA1):c.…`)

### Changed

- #114 - `_get_transcript_tags` now takes `tx_ac` explicitly (`_get_transcript_tags(tx_ac, transcript_data, genome_build)`) so overrides don't have to re-derive the accession from `transcript_data["id"]`
- #114 - Gene→transcript resolution now returns the canonical versioned accession: `LocalDataProvider.get_tx_ac_tags_for_gene` ranks, looks tags up by, and returns `transcript_data["id"]` (eg `NM_000059.4`) instead of whatever `_get_transcript_ids_for_gene` yielded - so providers whose gene→tx map is versionless (eg `NM_000059`) now resolve `BRCA2:c.36del` → `NM_000059.4:c.36del` (matching `resolve_gene_hgvs`'s documented output) rather than dropping the version; falls back to the id in hand when a record has no `"id"`. Providers whose ids are already versioned (eg `JSONDataProvider`) are unchanged
- #114 - `LocalDataProvider.get_tx_for_gene` now also returns the canonical versioned accession in `tx_ac` (`transcript_data["id"]`, falling back to the id in hand), so it agrees with `get_tx_ac_tags_for_gene` for versionless-id providers. biocommons `relevant_transcripts` maps via `get_tx_for_region` (not this method), so internal mapping is unaffected; already-versioned providers (eg `JSONDataProvider`) are unchanged

### Fixed

- #112 - HGVS cleaning no longer mangles LRG transcript references: the `t1`/`p1` transcript/protein suffix (eg `LRG_199t1(RUNX1):c.1415T>C`) is preserved through structure reconstruction instead of being mistaken for the gene symbol
- `JSONDataProvider` can now load 0.2.33+ data: the genome-build `source` field is accepted as either a string (early 0.2.32 data) or a list (eg `['BestRefSeq']`, 0.2.33 on). Previously `msgspec` rejected the list form with `Expected 'str | null', got 'array'`, so current production JSON.gz files failed to load. Mirrors the existing str-or-list handling of `biotype`
- #108, #109 - `RESTDataProvider.prefetch()` now respects server request-size limits: it POSTs the accession list to the batch `/transcripts` endpoint in chunks of `batch_size` (default `PREFETCH_BATCH_SIZE=1000`) instead of one giant request, and if a chunk is still rejected as too large (HTTP 413/414) it automatically halves the chunk and retries. Previously prefetching a large set (eg whole-ClinVar warming, ~tens of thousands of accessions) sent every id in a single POST and failed against the server's limit. Small inputs are unchanged (one request); the no-batch-endpoint fallback to concurrent singles is preserved

## [0.2.27] - 2026-06-16

### Added

- #86 Ensembl Tark Data Provider implementation 
- #27, #112 - HGVS cleaning: new `cdot.hgvs.clean` module (`clean_hgvs`, `get_best_transcript_version`) - cdot can now fix/clean common bad HGVS formatting (whitespace, quotes/punctuation, case, gene/transcript wrappers, separator typos, redundant del/dup counts, etc., driven by a real-world search-log corpus) and report warnings; cleaning ops are individually selectable via an `ops` set (`HGVSCleanOp` / `ALL_CLEAN_OPS`) with a `validate` flag, both forwarded by `fix_hgvs`
- #27, #36 - Gene HGVS: new `cdot.hgvs.gene_hgvs` module (`resolve_gene_hgvs`, `fix_hgvs`) - resolve gene-symbol HGVS (eg `BRCA2:c.36del`) via MANE/canonical transcript tags
- #28 - Transcript version fallback: `resolve_transcript_version` and a `version_fallback` arg on `fix_hgvs` (off by default) move to an adjacent transcript version (`VersionStrategy.UP_THEN_DOWN`/`CLOSEST`/`LATEST`) when the requested one isn't in the data, reporting a WARNING (`USED_ADJACENT_VERSION`); backed by a new `get_tx_versions(accession)` data-provider method on `JSONDataProvider` (enumerates in-memory) and `RESTDataProvider` (uses the versionless `/transcript/<ac>` lookup, which returns every stored version keyed by full accession, and warms the cache)
- #36 - Canonical transcript resolution: `get_tx_ac_tags_for_gene` (JSON + REST data providers) returns transcripts ranked by tag priority (MANE_Select > MANE_Plus_Clinical > RefSeq_Select > Ensembl_canonical > longest); `_get_transcript_tags` is overridable to supplement tags from an external source
- #36 - Added `Ensembl_canonical` tag to GRCh37 Ensembl data (pulled from the Ensembl REST API, see `generate_transcript_data/ensembl_grch37_canonical_transcripts.py`) so GRCh37 Ensembl transcripts can be resolved as canonical (only affects data not client code)
- #27 - `RESTDataProvider.get_tx_ac_tags_for_gene` - retrieve canonical/tagged transcripts over REST (needs cdot_rest endpoint `transcripts/gene/<gene>/tags/<genome_build>`, see SACGF/cdot_rest#12), enabling `resolve_gene_hgvs` against the REST API
- #37 - msgspec typed data models in `cdot.models`
- #109 - `RESTDataProvider.prefetch(tx_acs)` - read-ahead cache warming for bulk HGVS: populate the transcript cache up front so subsequent per-variant `c_to_g` calls are all cache hits. Plus `prefetch_from_hgvs(hgvs_strings)` which extracts transcript accessions (via `clean_hgvs`) and prefetches them
- #108 - `prefetch()` now warms the cache in a single round-trip via the batch `POST /transcripts` endpoint (needs cdot_rest, see SACGF/cdot_rest#9), expanding versionless accessions (eg `NM_000059`) to all versions server-side; falls back to a concurrent thread-pool of single `/transcript/<ac>` requests for servers without the batch endpoint

### Changed

- #27 - Gene HGVS resolution now matches RefSeq space-form tags (`MANE Select`) against the underscore-form tag priority, and retries a lowercase gene symbol (eg `brca2`) uppercased after a case-sensitive lookup miss
- #27 - `resolve_gene_hgvs` / `fix_hgvs` take a `prefer_consortium` arg (`Consortium.REFSEQ` default, or `ENSEMBL`/`None`) that **hard-filters** the resolved transcript's consortium - with a preference set you never get the other consortium back (eg `BRCA2` → `NM_000059.4` not `ENST00000380152.8`); errors if the preferred consortium has no transcript for the gene rather than crossing over; pass `None` to allow either
- #62 - Dropped Python 3.9 support (`requires-python` is now `>=3.10`); 3.9 reached end-of-life in Oct 2025 and the `build` extra's HTSeq>=2.1.2 requires 3.10+
- #111 - Gene `biotype` changed from comma-separated str (<= 0.2.19) to list (0.2.20) without a schema bump; `cdot.models` now accepts both forms and normalises legacy str biotype to `list[str]` on load
- #88 get_acs_for_protein_seq should return list not None
- #83 Ensembl files missing protein - breaking c_to_p (only affects data not client code)
- #17 RefSeq missing MT transcripts (only affects data not client code)
- #96 - All builds files now only contain 1 annotation consortium (RefSeq OR Ensembl not both)  
- #97 - Ensembl now has HGNC codes (used Gencode lookup) to match RefSeq
- #100 - Swapped cdot.cc domain references to cdotlib.org
- #106 - Collapse release nodes (only affects data not client code)
- Adjusted requirements - depend on hgvs core, bring in pysam for FASTA SeqFetcher
- Bug fixes + optimisations

### Fixed

- #27 - `clean_hgvs` no longer corrupts gene symbols that contain a mutation-type substring (eg `INSR`→`insR`, `INVS`→`invS`); mutation-type lowercasing is now confined to the allele (after the first `:`)
- #27 - `clean_hgvs` no longer strips brackets from valid HGVS with multiple *balanced* parentheses (gene symbol in parens plus uncertain-range notation, eg `NM_004006.2(DMD):c.(4071+1_4072-1)_(5154+1_5155-1)del`); only genuinely unbalanced brackets are stripped now
- #86 - `EnsemblTarkDataProvider` now reports coding transcripts that lack a 5' or 3' UTR correctly (CDS bounds derived from `cds_info`/UTR lengths), instead of returning `cds_start_i`/`cds_end_i` as `None` (which made them look non-coding) when either UTR was absent
- #36 - `JSONDataProvider.get_tx_ac_tags_for_gene` now ranks transcripts by spliced (exonic) length rather than genomic span (which included introns), matching the documented "decreasing transcript length" ordering used to pick the longest fallback transcript
- #36 - `resolve_gene_hgvs` now recognises a lowercase transcript accession (eg `enst00000617537.5:c.36del`, `nm_000059.4:c.36del`) as a transcript and passes it through unchanged, instead of treating it as a gene symbol and raising an error
- #86 - `EnsemblTarkDataProvider.get_tx_exons` now raises `HGVSDataNotAvailableError` (instead of crashing with a `TypeError`) when a transcript is not aligned to the requested (supported) contig
- #37 - `cdot.models.GenomeBuild` now retains the `source`, `ccds` and `transcript_support_level` build fields (data schema >= 0.2.32/0.2.33) instead of silently dropping them, preserving the dict drop-in contract
- #39 - `ExonsFromGenomeFastaSeqFetcher(cache=False)` now actually disables caching (the flag was being passed positionally and ignored)

## [0.2.26] 2024-08-15

Bumped version to 0.2.26 to catch up with data release. Only new client functionality is #81 'data_release' helper functions

All other changes in this release were for data (and contained in data_v0.2.26)

### Added

- #81 New 'data_release' code eg 'get_latest_combo_file_urls' that looks on GitHub to find latest data
- New GFFs: RefSeq RS_2023_10, Ensembl 111, 112
- #79 - RefSeq MT transcripts
- #66 - We now store 'Note' field (thanks holtgrewe for suggestion)
- client / JSON data schema version compatability check

### Changed

- #56 - Fix occasional UTA duplicated exons
- #57 - Correctly handle retrieving genomic position and dealing w/indels in GFF (thanks ltnetcase for reporting)
- #60 - Fix for missing protein IDs due to Genbank / GenBank (thanks holtgrewe)
- #64 - Split code/data versions. json.gz are now labelled according to data schema version (thanks holtgrewe)
- Renamed 'CHM13v2.0' to 'T2T-CHM13v2.0' so it could work with biocommons bioutils
- #72 - Correctly handle ncRNA_gene genes (thanks holtgrewe for reporting)
- #73 - HGNC ID was missing for some chrMT genes in Ensembl

## [0.2.21] - 2023-08-14

### Changed

- #45 - FastaSeqFetcher - fix alignment gaps properly
- #52 - Added transcripts from Ensembl 110 GRCh38 release
- #53 - UTA to cdot transcript start/end conversion issue  

## [0.2.20] - 2023-07-10

### Changed

- #50 - Biotype was missing in Ensembl transcripts 

## [0.2.19] - 2023-07-06

### Changed

- #49 - MT not converted to contigs correctly (GRCh37/Ensembl only) #49 
- Removed accidental logging

## [0.2.18] - 2023-07-05

### Added

- #44 - Support for mouse transcripts (Mus Musculus GRCm38 and GRCm39). Bash scripts only, never released as data, and removed in a later version (see #103)
- #47 - Implement HGVS DataProvider get_alignments_for_region

### Changed

- #45 - FastaSeqFetcher - handle deletions correctly (had swapped HGVS cigar projections around)
- #46 - HGVS DataProvider get_tx_info should properly handle alt_ac and alt_aln_method

## [0.2.17] - 2023-05-08

### Added

- #42 - Ensembl T2T CHM13v2.0

### Changed

- #43 - Contigs not converted to accession numbers properly (this was breaking local Biocommons HGVS conversion using 0.2.16 data)  

## [0.2.16] - 2023-04-12

### Added

- Added historical release 110 (2022-04-12) for T2T CHM13v2.0
- Added latest GRCh38.p14 release (2023-03-21)

## [0.2.15] - 2023-04-03

### Added

- Support for T2T CHM13v2.0

## [0.2.14] - 2023-03-21

### Added

- #39 - Fasta file SeqFetcher implementation
- Add Ensembl 109 GTF

### Changed

- #38 - Differing implementation of get_tx_for_region to hgvs one (reported by Manuel Holtgrewe) 
- #35 - Tags (ie MANE Select / RefSeq select etc) should be genome build specific 
- #34 - Stick to PyHGVS conventions, throw ValueError: transcript is required on missing transcript

## [0.2.13] - 2023-02-23

### Changed

- Fix for #25 - Pyhgvs data conversion - non-coding transcripts have bad cds start/end conversion
- Fix for #32 - Signature of get_pyhgvs_data consistent for all return statements

## [0.2.12] - 2022-12-08

### Added

- #30 - We now store "tag" attributes (eg "MANE Select", "RefSeq Select")
- Switch to using Ensembl GFF3 (so we can get tags out)
- Add Ensembl 108 GFF3

### Changed

- Fix for #25 -  GeneInfo currently fails for some records
- Fix for #27 -  Change URL for missing RefSeq GFFs

## [0.2.11] - 2022-09-27

### Added

- Now support all methods (get_gene_info, get_tx_for_gene, get_tx_for_region) for REST
- Add Ensembl 107 GTF

### Changed

- Ensembl gene info was missing "description"

## [0.2.10] - 2022-09-19

### Added

- [Implement get_gene_info](https://github.com/SACGF/cdot/issues/20) - For local JSON data only

### Changed

- Fixed issue [#23 UTA transcripts for PyHGVS](https://github.com/SACGF/cdot/issues/23)

## [0.2.9] - 2022-09-01

### Changed

- [BugFix for get_tx_for_region](https://github.com/SACGF/cdot/issues/22)


## [0.2.8] - 2022-08-29

### Added

- [Implemented get_pro_ac_for_tx_ac](https://github.com/SACGF/cdot/issues/14) (c_to_p can now generate p.HGVS)
- [Implemented get_tx_for_region](https://github.com/SACGF/cdot/issues/18) for local JSON data only

## [0.2.7] - 2022-05-19

### Added

- Add transcripts from latest RefSeq GRCh37 (105) and RefSeq GRCh38 (110)

### Changed

- Fixed default arguments bug where PyHGVS only worked on SACGF fork
- gtf_to_json now goes straight to cdot format (without intermediary PyReference format)
- UTA is not included in generation scripts by default, to enable, set environment variable UTA_TRANSCRIPTS=True
- Handle mismatches in UTA CIGAR alignments (convert to match (no indels) as GFF format has no support for mismatch)

## [0.2.6] - 2022-05-19

### Changed

- Fixed issue [Ensembl contigs g_to_c](https://github.com/SACGF/cdot/issues/9) - Ensembl JSON was using chrom names ie "17" instead of "NC_000017.11" for contig 

## [0.2.5] - 2022-04-14

### Changed

- PyHGVS conversion fix - non-coding cds_start/cds_end is set to start/end (not None)

## [0.2.4] - 2022-04-13

### Added

- Latest RefSeq (110) and Ensembl (106) transcripts

### Changed

- Fixed bug where all UTA transcripts were '-' strand
- Add "other_chroms" to combined historical file

## [0.2.3] - 2022-03-29

### Changed

- Fixed bug where HGNC not extracted properly from Ensembl GTFs
- Gene information is now included by default (only adds 5%)
- Clean artifacts from UTA data
- Support for [SACGF PyHGVS fork](https://github.com/SACGF/hgvs) (which adds alignment gap support)

## [0.2.2] - 2022-03-03

### Added

- Support for HTTPS (bought SSL certificate for REST server)

## [0.2.1] - 2022-03-03

### Added

- [Download/Convert UTA transcripts](https://github.com/SACGF/cdot/issues/1)
- [REST client](https://github.com/SACGF/cdot/issues/4) for [REST Service](https://github.com/SACGF/cdot_rest/)

### Changed

- [JSON format changed](https://github.com/SACGF/cdot/issues/2), separating common/build specific coordinates. This is so a transcript can contain data for multiple builds.
- [Use ijson to reduce RAM usage](https://github.com/SACGF/cdot/issues/7) - uses iterator vs loading all JSON into RAM

## [0.1.1] - 2022-01-19

### Added

- Initial commit

[unreleased]: https://github.com/SACGF/cdot/compare/v0.2.33...HEAD
[0.2.33]: https://github.com/SACGF/cdot/compare/v0.2.32...v0.2.33
[0.2.32]: https://github.com/SACGF/cdot/compare/v0.2.31...v0.2.32
[0.2.31]: https://github.com/SACGF/cdot/compare/v0.2.30...v0.2.31
[0.2.30]: https://github.com/SACGF/cdot/compare/v0.2.28...v0.2.30
[0.2.28]: https://github.com/SACGF/cdot/compare/v0.2.27...v0.2.28
[0.2.27]: https://github.com/SACGF/cdot/compare/v0.2.26...v0.2.27
[0.2.26]: https://github.com/SACGF/cdot/compare/v0.2.21...v0.2.26
[0.2.21]: https://github.com/SACGF/cdot/compare/v0.2.20...v0.2.21
[0.2.20]: https://github.com/SACGF/cdot/compare/v0.2.19...v0.2.20
[0.2.19]: https://github.com/SACGF/cdot/compare/v0.2.18...v0.2.19
[0.2.18]: https://github.com/SACGF/cdot/compare/v0.2.17...v0.2.18
[0.2.17]: https://github.com/SACGF/cdot/compare/v0.2.16...v0.2.17
[0.2.16]: https://github.com/SACGF/cdot/compare/v0.2.15...v0.2.16
[0.2.15]: https://github.com/SACGF/cdot/compare/v0.2.14...v0.2.15
[0.2.14]: https://github.com/SACGF/cdot/compare/v0.2.13...v0.2.14
[0.2.13]: https://github.com/SACGF/cdot/compare/v0.2.12...v0.2.13
[0.2.12]: https://github.com/SACGF/cdot/compare/v0.2.11...v0.2.12
[0.2.11]: https://github.com/SACGF/cdot/compare/v0.2.10...v0.2.11
[0.2.10]: https://github.com/SACGF/cdot/compare/v0.2.8...v0.2.10
[0.2.8]: https://github.com/SACGF/cdot/compare/v0.2.7...v0.2.8
[0.2.7]: https://github.com/SACGF/cdot/compare/v0.2.6...v0.2.7
[0.2.6]: https://github.com/SACGF/cdot/compare/v0.2.5...v0.2.6
[0.2.5]: https://github.com/SACGF/cdot/compare/v0.2.4...v0.2.5
[0.2.4]: https://github.com/SACGF/cdot/compare/v0.2.3...v0.2.4
[0.2.3]: https://github.com/SACGF/cdot/compare/v0.2.2...v0.2.3
[0.2.2]: https://github.com/SACGF/cdot/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/SACGF/cdot/compare/v0.1.1...v0.2.1
[0.1.1]: https://github.com/SACGF/cdot/releases/tag/v0.1.1
