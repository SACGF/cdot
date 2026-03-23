# Discussion

*Target: ~800–1000 words. Interpret results, position vs related work, limitations, future.*

---

## 4.1 cdot in the biocommons ecosystem

cdot plays the same role for transcript coordinates that SeqRepo [Hart 2020] plays for biological sequences in the biocommons stack: both replace slow, infrastructure-heavy remote database access with local, high-performance flat-file storage. The combination of cdot (transcript coordinates) + SeqRepo (sequences) + biocommons/hgvs (parsing and conversion engine) constitutes a complete, offline-capable HGVS processing stack requiring no running services.

## 4.2 Comparison with related resources

**UTA**: cdot provides [9]× more versioned transcript alignments than UTA without requiring PostgreSQL. The key structural difference is that cdot regenerates data from primary annotation sources (RefSeq GFF3, Ensembl GTF) using a reproducible pipeline, while UTA is a manually curated snapshot. cdot uses UTA CSV dumps as a baseline for accessions predating its oldest GFF3/GTF files, preserving UTA's historical coverage while extending it.

**Ensembl TARK**: The Ensembl Transcript Archive (TARK) provides REST-based access to Ensembl transcript data. Unlike cdot, TARK is online-only, covers Ensembl transcripts only (no RefSeq), and does not store RefSeq alignment gap strings. cdot's local file approach supports air-gapped clinical environments where TARK is inaccessible.

**VariantValidator** [Freeman 2018]: Provides web-based HGVS validation built on biocommons/hgvs. VariantValidator is architecturally a monolithic service; cdot is an embeddable, offline data layer. The two are complementary: cdot provides the transcript data, VariantValidator provides validation logic. VariantValidator could in principle use cdot as its data provider.

**TransVar** [Zhou 2015]: Bundles its own annotation and provides multilevel variant conversion. Unlike cdot, TransVar does not separate the data generation pipeline from the conversion tool, limiting flexibility in transcript version management.

## 4.3 The shared transcript data layer argument

McCarthy et al. [2014] and Park et al. [2022] demonstrate that transcript set and tool choice can affect variant annotations for 15–56% of variants. Some of this disagreement stems from tools using different transcript versions for the same gene. cdot provides a single, versioned, source-attributed data source that any HGVS tool could use — not as a replacement for existing tools, but as a shared data layer that reduces one source of annotation disagreement. VEP-compatible cdot data releases (planned, issue #63) would enable cross-referencing between VEP annotations and HGVS-specific tools using matched transcript versions.

## 4.4 MANE and canonical transcript selection

MANE Select now covers >97% of protein-coding genes [Morales 2022], and clinical laboratory standards are increasingly mandating MANE transcript reporting [Wright 2023]. cdot's storage of MANE Select, MANE Plus Clinical, RefSeq Select, and Ensembl canonical tags, combined with the planned `CanonicalTranscriptSelector` API, makes cdot the first HGVS data provider to offer programmatic canonical transcript lookup. This transforms cdot from a passive data source into an active contributor to MANE adoption.

## 4.5 Pangenome and future reference assemblies

The Human Pangenome Reference [Liao 2023] represents 47 phased diploid assemblies, adding 119 million base pairs of polymorphic sequence. As variant calling against pangenome haplotypes becomes more common, HGVS tools will need transcript annotation for diverse reference sequences. cdot's multi-build architecture — treating each genome assembly as an independently indexed set of transcript alignments within a shared JSON structure — is designed to support additional assemblies without format changes.

## 4.6 Limitations

[List honest limitations:]
- **CDS phase / ribosomal frameshifts**: cdot does not currently model CDS phase at the exon level or ribosomal frameshifts, which can affect protein HGVS conversion for a small number of transcripts (issue #76).
- **Ensembl GFF3 protein versions**: Ensembl GFF3 files contain broken protein version information; cdot currently uses GTF files for Ensembl data as a workaround (issue #95).
- **Validation**: cdot provides coordinate data but does not validate HGVS descriptions — validation should be performed by a dedicated tool (VariantValidator, Mutalyzer).
- **Single-genome scope**: cdot does not yet support haplotype-aware transcript coordinates for pangenome references.
- **Non-HGVS use**: The JSON format has broader utility for general bioinformatics (gene annotation, splice analysis). The pyreference library (https://github.com/SACGF/pyreference) demonstrates this, but non-HGVS use cases are outside the scope of this paper.

## 4.7 Availability and maintenance

cdot data files are generated from publicly available annotation releases using a reproducible Snakemake pipeline. New data files will be released with each major Ensembl and RefSeq annotation release. The REST API at cdotlib.org is maintained by SACGF. All code is MIT-licensed and available at https://github.com/SACGF/cdot.

---

*Notes:*
- *Fill section 4.2 UTA comparison with exact transcript count once benchmarks are run*
- *Section 4.6 CDS phase limitation — check if this will be fixed before submission*
- *Consider whether to mention pyreference in more or less detail depending on word count*
