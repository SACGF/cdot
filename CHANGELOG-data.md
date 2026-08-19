# cdot data changelog

Changes to the transcript data published as GitHub releases (`data_v<version>`), which is
versioned independently of the pip package. The version here is
[`JSON_SCHEMA_VERSION`](generate_transcript_data/json_schema_version.py); it is embedded in every
output filename and becomes the release tag. Clients check compatibility on major/minor only, so a
patch bump is always readable by an existing client.

For changes to the Python library itself see [CHANGELOG.md](CHANGELOG.md). For how a data release
is built and published see [docs/data_release_workflow.md](docs/data_release_workflow.md).

## [unreleased]

## [0.2.34] - 2026-08-13

### Added

- #51 - RefSeq GRCh38 data now includes NCBI's historical transcript alignments (`RefSeq_historical_alignments`, RS_2024_08 set): alignments of replaced and suppressed NM_/NR_ versions, many of which never appeared in any annotation release cdot ingests. They are merged at the lowest priority above UTA, so a transcript version from any official annotation release is unchanged. Note some historical alignments are partial: a few transcript bases at the ends (eg an unaligned poly-A tail) or occasionally the first bases do not align to the genome, and the stored exon cDNA coordinates reflect that, so positions in an unaligned region cannot be projected

### Fixed

- #123 - `start_codon`/`stop_codon` are now positions along the whole transcript, so they agree with the exon `cds_start`/`cds_end`. A few RefSeq alignments leave a run of transcript bases unaligned between two exons (58 multi-exon GRCh38 records, eg `NM_033517.1`), and the generator collapsed those out of the codon positions but not out of the exon coordinates, so the CDS was reported short and its last codons were rejected as outside CDS bounds. Also corrects the documented convention: the codon fields are 0-based half-open (the biocommons `cds_start_i`/`cds_end_i` they are passed through as), not 1-based
- #95 - The genomic `cds_start`/`cds_end` of UTA-sourced minus-strand transcripts were both 1 too high (a 0-based/1-based conversion error in the UTA import; plus-strand transcripts were unaffected). These fields are only used by the PyHGVS integration, so biocommons HGVS resolution was never affected. The UTA import also no longer requires the SACGF PyHGVS fork: when the pipeline ran with standard pyhgvs installed instead, every coding UTA transcript failed conversion and was silently skipped, which is why release 0.2.33 has only 3,954 UTA-sourced GRCh37 records (all non-coding) against 46,659 in 0.2.28. Coding UTA transcripts are restored here, and the pipeline's UTA fetch now fails loudly instead of keeping a truncated download

## [0.2.33] - 2026-06-26

### Added

- Store `ccds` and `transcript_support_level`

## [0.2.32] - 2025-10-08

### Added

- `source` (GTF column #2) added to the genome build data

## [0.2.31] - 2025-08-18

### Added

- `metadata` added (generation method and source urls)

## [0.2.30] - 2025-07-29

### Added

- Ensembl GRCh37 has canonical transcripts added from outside the GTFs

## [0.2.29] - 2025-07-18

### Added

- Ensembl now has HGNC added from outside the GTFs

<!-- Data releases before 0.2.29 were not separately recorded; see the git history of
     generate_transcript_data/ for those. -->

[unreleased]: https://github.com/SACGF/cdot/compare/data_v0.2.34...HEAD
[0.2.34]: https://github.com/SACGF/cdot/compare/data_v0.2.33...data_v0.2.34
[0.2.33]: https://github.com/SACGF/cdot/compare/data_v0.2.32...data_v0.2.33
[0.2.32]: https://github.com/SACGF/cdot/compare/data_v0.2.31...data_v0.2.32
[0.2.31]: https://github.com/SACGF/cdot/compare/data_v0.2.30...data_v0.2.31
[0.2.30]: https://github.com/SACGF/cdot/compare/data_v0.2.29...data_v0.2.30
[0.2.29]: https://github.com/SACGF/cdot/releases/tag/data_v0.2.29
