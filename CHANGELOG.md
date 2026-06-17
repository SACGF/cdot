## [unreleased]

### Added

- #114 - Public HGVS ranking primitives + provider hooks (to let consumers stop importing private symbols / re-implementing ranking): `rank_transcript_versions` (the full version ordering behind `get_best_transcript_version`) and `rank_transcripts_for_gene` (the gene→transcript lookup/uppercase-retry/consortium-filter/tag-ranking core split out of `resolve_gene_hgvs`, returning the full ranked candidate list); `consortium_of` (promoted from `_consortium_of`); `HGVSFix.__str__` plus `messages`/`warning_messages`/`error_messages` helpers; an overridable batch tag hook `_get_tags_by_tx_ac(tx_acs, genome_build)` on data providers (default loops the per-transcript hook) so subclasses can answer in one query instead of N+1; and `cdot.hgvs` now also re-exports `Consortium`, `DEFAULT_CONSORTIUM`, `HGVSCleanOp`, and `ALL_CLEAN_OPS`
- #112 - HGVS cleaning: new ops driven by the search-log corpus - colon used in place of the kind dot (`:c:1811` → `:c.1811`), parens wrapping the whole accession (`(NM_000548.3):c.…` → `NM_000548.3:c.…`), and a gene-first wrapper missing the colon before the kind (`MYB(NM_…)c.…` → `NM_…(MYB):c.…`)

### Changed

- #114 - `_get_transcript_tags` now takes `tx_ac` explicitly (`_get_transcript_tags(tx_ac, transcript_data, genome_build)`) so overrides don't have to re-derive the accession from `transcript_data["id"]`
- #114 - Gene→transcript resolution now returns the canonical versioned accession: `LocalDataProvider.get_tx_ac_tags_for_gene` ranks, looks tags up by, and returns `transcript_data["id"]` (eg `NM_000059.4`) instead of whatever `_get_transcript_ids_for_gene` yielded - so providers whose gene→tx map is versionless (eg `NM_000059`) now resolve `BRCA2:c.36del` → `NM_000059.4:c.36del` (matching `resolve_gene_hgvs`'s documented output) rather than dropping the version; falls back to the id in hand when a record has no `"id"`. Providers whose ids are already versioned (eg `JSONDataProvider`) are unchanged
- #114 - `LocalDataProvider.get_tx_for_gene` now also returns the canonical versioned accession in `tx_ac` (`transcript_data["id"]`, falling back to the id in hand), so it agrees with `get_tx_ac_tags_for_gene` for versionless-id providers. biocommons `relevant_transcripts` maps via `get_tx_for_region` (not this method), so internal mapping is unaffected; already-versioned providers (eg `JSONDataProvider`) are unchanged

### Fixed

- #112 - HGVS cleaning no longer mangles LRG transcript references: the `t1`/`p1` transcript/protein suffix (eg `LRG_308t1(PALB2):c.3113G>A`) is preserved through structure reconstruction instead of being mistaken for the gene symbol

## [0.2.27] - 2026-06-16

### Added

- #62 - Run the test suite in GitHub Actions CI (`.github/workflows/tests.yml`) on push/PR across Python 3.10–3.14 (dev/infra only, no client code change)
- #86 Ensembl Tark Data Provider implementation 
- #70 Use Snakemake to build transcripts (only affects data not client code)
- #27, #112 - HGVS cleaning: new `cdot.hgvs.clean` module (`clean_hgvs`, `get_best_transcript_version`) - cdot can now fix/clean common bad HGVS formatting (whitespace, quotes/punctuation, case, gene/transcript wrappers, separator typos, redundant del/dup counts, etc., driven by a real-world search-log corpus) and report warnings; cleaning ops are individually selectable via an `ops` set (`HGVSCleanOp` / `ALL_CLEAN_OPS`) with a `validate` flag, both forwarded by `fix_hgvs`
- #27, #36 - Gene HGVS: new `cdot.hgvs.gene_hgvs` module (`resolve_gene_hgvs`, `fix_hgvs`) - resolve gene-symbol HGVS (eg `BRCA2:c.36del`) via MANE/canonical transcript tags
- #28 - Transcript version fallback: `resolve_transcript_version` and a `version_fallback` arg on `fix_hgvs` (off by default) move to an adjacent transcript version (`VersionStrategy.UP_THEN_DOWN`/`CLOSEST`/`LATEST`) when the requested one isn't in the data, reporting a WARNING (`USED_ADJACENT_VERSION`); backed by a new `get_tx_versions(accession)` data-provider method on `JSONDataProvider` (enumerates in-memory) and `RESTDataProvider` (uses the versionless `/transcript/<ac>` lookup, which returns every stored version keyed by full accession, and warms the cache)
- #36 - Canonical transcript resolution: `get_tx_ac_tags_for_gene` (JSON + REST data providers) returns transcripts ranked by tag priority (MANE_Select > MANE_Plus_Clinical > RefSeq_Select > Ensembl_canonical > longest); `_get_transcript_tags` is overridable to supplement tags from an external source
- #36 - Added `Ensembl_canonical` tag to GRCh37 Ensembl data (pulled from the Ensembl REST API, see `generate_transcript_data/ensembl_grch37_canonical_transcripts.py`) so GRCh37 Ensembl transcripts can be resolved as canonical (only affects data not client code)
- #27 - `RESTDataProvider.get_tx_ac_tags_for_gene` - retrieve canonical/tagged transcripts over REST (needs cdot_rest endpoint `transcripts/gene/<gene>/tags/<genome_build>`, see SACGF/cdot_rest#12), enabling `resolve_gene_hgvs` against the REST API
- #37 - msgspec typed data models in `cdot.models`
- #77 - Generate JSON schema docs (`generate_transcript_data/generate_json_docs.py`)
- #109 - `RESTDataProvider.prefetch(tx_acs)` - read-ahead cache warming for bulk HGVS: populate the transcript cache up front so subsequent per-variant `c_to_g` calls are all cache hits. Plus `prefetch_from_hgvs(hgvs_strings)` which extracts transcript accessions (via `clean_hgvs`) and prefetches them
- #108 - `prefetch()` now warms the cache in a single round-trip via the batch `POST /transcripts` endpoint (needs cdot_rest, see SACGF/cdot_rest#9), expanding versionless accessions (eg `NM_000059`) to all versions server-side; falls back to a concurrent thread-pool of single `/transcript/<ac>` requests for servers without the batch endpoint
- #102 - Docs: new `docs/coordinates_and_exons.md` (how exon coordinates and the alignment gap strings work) and `docs/advanced_usage.md` (fixing messy HGVS input via `fix_hgvs`/`clean_hgvs`, and bulk read-ahead retrieval via `prefetch`), plus a `docs/` index; linked from the README (docs only, no client code change)
- #102 - Docs: migrated the GitHub wiki into `docs/` (examples, FastaSeqFetcher, local JSON files, release file details, create-data-from-scratch, cdot vs UTA, design notes) so docs are versioned with the code; README and the generated `json_data_format.md` now link to `docs/` instead of the wiki (docs only, no client code change)

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
- Added requirements.txt for 'generate_transcript_data' sections
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

- #44 - Support for mouse transcripts (Mus Musculus GRCm38 and GRCm39)
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

[unreleased]: https://github.com/SACGF/cdot/compare/v0.2.26...HEAD
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
[0.2.10]: https://github.com/SACGF/cdot/compare/v0.2.9...v0.2.10
[0.2.9]: https://github.com/SACGF/cdot/compare/v0.2.8...v0.2.9
[0.2.8]: https://github.com/SACGF/cdot/compare/v0.2.7...v0.2.8
[0.2.7]: https://github.com/SACGF/cdot/compare/v0.2.6...v0.2.7
[0.2.6]: https://github.com/SACGF/cdot/compare/v0.2.5...v0.2.6
[0.2.5]: https://github.com/SACGF/cdot/compare/v0.2.4...v0.2.5
[0.2.4]: https://github.com/SACGF/cdot/compare/v0.2.3...v0.2.4
[0.2.3]: https://github.com/SACGF/cdot/compare/v0.2.2...v0.2.3
[0.2.2]: https://github.com/SACGF/cdot/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/SACGF/cdot/compare/v0.1.1...v0.2.1
[0.1.1]: https://github.com/SACGF/cdot/releases/tag/v0.1.1
