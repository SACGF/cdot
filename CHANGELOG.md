## [unreleased]

### Added

- New GFFs: RefSeq RS_2023_10, Ensembl 110, 111
- #66 - We now store 'Note' field (thanks holtgrewe for suggestion)
- Added requirements.txt for 'generate_transcript_data' sections
- client / JSON data schema version compatability check

### Changed

- #56 - Fix occasional UTA duplicated exons
- #57 - Correctly handle retrieving genomic position and dealing w/indels in GFF (thanks ltnetcase for reporting)
- #60 - Fix for missing protein IDs due to Genbank / GenBank (thanks holtgrewe)
- #64 - Split code/data versions. json.gz are now labelled according to data schema version (thanks holtgrewe)
- Renamed 'CHM13v2.0' to 'T2T-CHM13v2.0' so it could work with biocommons bioutils

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

[unreleased]: https://github.com/SACGF/cdot/compare/v0.2.21...HEAD
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
