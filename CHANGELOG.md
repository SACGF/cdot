## Unreleased [unreleased]

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

[unreleased]: https://github.com/SACGF/cdot/compare/v0.2.7...HEAD
[0.2.7]: https://github.com/SACGF/cdot/compare/v0.2.6...v0.2.7
[0.2.6]: https://github.com/SACGF/cdot/compare/v0.2.5...v0.2.6
[0.2.5]: https://github.com/SACGF/cdot/compare/v0.2.4...v0.2.5
[0.2.4]: https://github.com/SACGF/cdot/compare/v0.2.3...v0.2.4
[0.2.3]: https://github.com/SACGF/cdot/compare/v0.2.2...v0.2.3
[0.2.2]: https://github.com/SACGF/cdot/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/SACGF/cdot/compare/v0.1.1...v0.2.1
[0.1.1]: https://github.com/SACGF/cdot/releases/tag/v0.1.1
