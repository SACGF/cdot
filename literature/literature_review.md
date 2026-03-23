# Literature Review: Context for cdot

## 1. The HGVS Nomenclature Standard

The Human Genome Variation Society (HGVS) nomenclature is the established standard for describing sequence variants in clinical reports, publications, and databases [den Dunnen 2016, PMID 26931183]. HGVS strings (e.g., `NM_000492.3(CFTR):c.1521_1523delCTT`) express variants in terms of a named transcript and a position within that transcript's coordinate system. This makes the nomenclature human-interpretable and clinically precise, but creates a fundamental computational challenge: converting between HGVS transcript coordinates and the genomic coordinates used by sequencing pipelines requires authoritative, versioned transcript data.

The scale of this problem is large. Mutalyzer, an HGVS syntax checker, processed over 133 million variant descriptions in five years, finding approximately 50% of submitted descriptions contained errors [Lefter 2021, PMID 33538839] — a rate consistent with the original Mutalyzer study that found 13–75% error rates across clinical mutation databases [Wildeman 2008, PMID 18000842]. ClinVar, the primary archive of clinically interpreted variants, now holds >3 million variants from >2,800 organisations [Landrum 2025, PMID 39578691], the vast majority described using HGVS transcript notation. The ACMG/AMP variant interpretation guidelines [Richards 2015, PMID 25741868], now the de-facto standard for clinical variant classification worldwide, assume reliable HGVS coordinate conversion. Together, these numbers quantify the clinical demand that HGVS infrastructure must serve.

---

## 2. Existing Python HGVS Tools

Two Python libraries dominate programmatic HGVS handling in research and clinical pipelines.

**biocommons/hgvs** [Hart 2015, PMID 25273102; Wang 2018, PMID 30129167] is the most comprehensive open-source HGVS library, supporting parsing, formatting, validating, normalising, and projecting variants between aligned sequences including those with gapped alignments. It obtains its transcript data from the Universal Transcript Archive (UTA), a PostgreSQL database. This creates two practical barriers: UTA requires either a local PostgreSQL installation or a remote connection to a public server, and firewalls frequently block outbound PostgreSQL connections in clinical environments. UTA also lacks Ensembl transcript version support, limiting the library to RefSeq transcripts despite Ensembl's widespread use in research genomics.

**PyHGVS** (Counsyl) is a simpler, more lightweight library with a smaller feature set. It lacks alignment gap support, which causes incorrect coordinate conversion for transcripts that do not align perfectly to the reference genome — a problem that affects a meaningful fraction of RefSeq transcripts that contain documented indels relative to the genome sequence.

Both libraries were designed around a pluggable data provider model, leaving open the question of what provides transcript data. cdot is designed to fill this role.

---

## 3. Other HGVS Validation and Conversion Tools

**VariantValidator** [Freeman 2018, PMID 28967166] provides web-based HGVS validation built on biocommons/hgvs, with additional error correction, claiming superior accuracy over competing tools. Like cdot-powered pipelines, it uses the HGVS standard and biocommons infrastructure, but is architecturally a monolithic service requiring its own deployment rather than a lightweight embeddable library.

**Mutalyzer** [Wildeman 2008, PMID 18000842; Lefter 2021, PMID 33538839] focuses on syntax checking and correction of HGVS descriptions. It does not provide HGVS-to-genomic coordinate conversion as its primary function.

**TransVar** [Zhou 2015, PMID 26513549] converts between genomic, transcript, and protein-level representations, supporting RefSeq, Ensembl, GENCODE, and CCDS. It bundles its own annotation rather than using an external data provider, which limits flexibility in transcript version management.

**CSN and CAVA** [Münz 2015, PMID 26315209] introduce Clinical Sequencing Nomenclature, a deterministic HGVS-aligned scheme, and a fast annotation tool achieving 100% concordance on BRCA1/BRCA2 variants versus other tools' 32% — illustrating the severity of annotation correctness problems in existing tools.

**SnpEff** [Cingolani 2012, PMID 22728672] and **ANNOVAR** [Wang 2010, PMID 20601685] are the dominant open-source variant consequence prediction tools. Both bundle their own transcript databases and do not use HGVS as a primary format. As shown by McCarthy et al. [2014, PMID 24944579], only 44% of putative loss-of-function variants agreed between RefSeq and Ensembl annotation sets when processed by ANNOVAR, and only 65% agreed between ANNOVAR and VEP for loss-of-function annotations. Park et al. [2022, PMID 34612497] found similar 85% agreement between ANNOVAR and SnpEff on clinical variants. These papers establish that transcript set and tool choice are not minor implementation details but have major consequences for research and clinical conclusions.

The **Ensembl Variant Effect Predictor (VEP)** [McLaren 2016, PMID 27268795] is the most widely used variant annotation framework globally, integrating annotation from Ensembl, GENCODE, and RefSeq gene sets. VEP-compatible cdot data releases (issue #63) would allow users to cross-reference VEP annotations with HGVS-specific tools using matched transcript versions.

---

## 4. Transcript Data Sources

**RefSeq** [O'Leary 2016, PMID 26553804; Goldfarb 2025, PMID 39526381] maintains curated, non-redundant reference sequences for genomes, transcripts, and proteins. Human transcript accessions (NM_, NR_, NP_) are the authoritative source for clinical HGVS nomenclature — transcript versions like NM_000492.3 appear in clinical reports, databases, and publications worldwide. RefSeq's 25-year history of curation and annotation has produced a rich set of historical transcript versions, each potentially representing variants from patients seen at that time. cdot consumes multiple historical RefSeq GFF3 annotation releases to support these historical HGVS strings.

**Ensembl** [Martin 2023, PMID 36318249] is the complementary major annotation source, widely used in research genomics (particularly in Europe). Ensembl and RefSeq have historically annotated different transcripts for the same genes, making multi-source coverage essential for comprehensive HGVS resolution. Ensembl transcript accessions (ENST) appear in research publications and some clinical contexts.

**GENCODE** [Frankish 2023, PMID 36420896] is the high-quality manually curated gene annotation used by Ensembl for human and mouse. Recent convergence between GENCODE, RefSeq, and UniProt protein-coding gene annotations reduces (but does not eliminate) differences between the two major sources. cdot already uses GENCODE HGNC files for gene symbol/HGNC ID mapping, and open issue #90 considers using GENCODE GTFs directly as a transcript source.

**Universal Transcript Archive (UTA)** is the biocommons project that cdot partially replaces. UTA stores transcript alignments in PostgreSQL; cdot replaces the PostgreSQL dependency with flat JSON files while providing 9× more transcript versions and adding Ensembl support. UTA remains a valuable source of historical transcript alignments, particularly for transcripts that predate the oldest GFF3 files cdot processes — cdot ingests UTA CSV dumps as a baseline source, overridden by official annotation where available.

---

## 5. Canonical Transcript Selection

The lack of a single "standard" transcript per gene has been a longstanding problem in genomics. For any given gene, RefSeq and Ensembl historically provided different representative transcripts, causing inconsistency in clinical reports and variant databases.

The **MANE** (Matched Annotation from NCBI and EMBL-EBI) initiative [Morales 2022, PMID 35388217] addresses this by jointly selecting one representative transcript per protein-coding gene (MANE Select) and additional clinically relevant transcripts (MANE Plus Clinical). MANE Select now covers 97% of protein-coding genes, including all ACMG secondary findings genes. The clinical importance of adopting MANE standards has been argued forcefully [Wright 2023, PMID 36441169], with calls for laboratory standards to mandate MANE transcript reporting.

Evidence supports MANE and **APPRIS** principal isoforms as biologically meaningful choices for canonical transcripts [Pozo 2022, PMID 36124785] — these approaches coincide with the main proteomics isoform for >98.2% of genes, while simpler heuristics (longest transcript) perform poorly.

cdot already stores MANE Select, MANE Plus Clinical, RefSeq Select, and Ensembl_canonical tags in build-specific transcript data (GRCh38 only for MANE; GRCh37 canonical Ensembl transcripts added in v0.2.30). Implementing a `CanonicalTranscriptSelector` API that surfaces these tags (issue #36) would directly enable gene-name HGVS lookup for the clinical community.

---

## 6. Reference Genome Evolution

The **T2T-CHM13** assembly [Nurk 2022, PMID 35357919] presents the first complete (telomere-to-telomere) human genome, resolving ~8% of the genome previously missing from GRCh38, including all centromeric satellite arrays and the short arms of acrocentric chromosomes. The assembly adds ~200 million base pairs and ~1,956 gene predictions. cdot is among the first HGVS tools to support T2T-CHM13v2.0 as a genome build, using both Ensembl and RefSeq annotations.

The **Human Pangenome Reference** [Liao 2023, PMID 37165242] extends this further, representing 47 phased diploid assemblies from genetically diverse individuals. Adding 119 million base pairs of polymorphic sequence and 1,115 gene duplications, the pangenome reduces small variant errors by 34% and increases structural variant detection by 104% compared to GRCh38. This creates new challenges for HGVS tooling — variants called against pangenome haplotypes will require transcript annotation for diverse reference sequences, motivating cdot's extensible multi-build architecture.

---

## 7. Sequence Data Infrastructure

**SeqRepo** [Hart 2020, PMID 33270643] provides local high-performance biological sequence storage, achieving up to 1300× speedup over remote sequence retrieval with a 50× throughput improvement for variant validation pipelines. SeqRepo is the sequence storage complement to biocommons/hgvs, providing actual nucleotide sequences that the HGVS library needs for variant normalisation and validation. cdot's `FastaSeqFetcher` provides an alternative local sequence source using genome FASTA files, avoiding SeqRepo's installation overhead while retaining offline capability.

---

## 8. Synthesis: Where cdot Fits

The cdot project addresses a specific and well-motivated gap in the HGVS tool ecosystem. The existing landscape has two dominant problems:

**1. Data access.** UTA (the standard biocommons/hgvs data backend) requires PostgreSQL, is blocked by firewalls in clinical environments, lacks Ensembl support, and provides only ~141,000 versioned transcript alignments. cdot provides the same data as a flat JSON file — loadable in any environment, offline-capable, and covering >1.3 million versioned transcript/genome alignments across both RefSeq and Ensembl.

**2. Coverage fragmentation.** The McCarthy 2014 and Park 2022 papers demonstrate that transcript set choice (RefSeq vs Ensembl) and tool choice dramatically affect annotation outcomes. cdot normalises this by providing a single, comprehensive, versioned data source that covers both consortia.

The broader ecosystem context is one of increasing standardisation pressure: MANE transcripts [Morales 2022; Wright 2023] are driving convergence on canonical transcript selection, ACMG/AMP guidelines [Richards 2015] have made HGVS accuracy a clinical requirement, and new reference genomes (T2T [Nurk 2022], pangenome [Liao 2023]) are expanding the genomic landscape that HGVS tools must cover. cdot's architecture — a simple, versioned JSON format served by a REST API or loaded locally, with Python integrations for both major HGVS libraries — is well-positioned to grow with these demands.

The key gap cdot fills that no existing tool addresses: **comprehensive, versioned, multi-source transcript data accessible without database infrastructure, at speeds suitable for production clinical pipelines.**

---

## References (by PMID)

| PMID | Citation |
|------|---------|
| 18000842 | Wildeman et al. (2008) Mutalyzer. Hum Mutat 29:6-13 |
| 20601685 | Wang et al. (2010) ANNOVAR. Nucleic Acids Res 38:e164 |
| 22728672 | Cingolani et al. (2012) SnpEff. Fly 6:80-92 |
| 14_mccarthy | McCarthy et al. (2014) Transcript choice effects. Genome Med 6:26 |
| 25273102 | Hart et al. (2015) biocommons/hgvs. Bioinformatics 31:268-70 |
| 25741868 | Richards et al. (2015) ACMG/AMP guidelines. Genet Med 17:405-24 |
| 26315209 | Münz et al. (2015) CSN/CAVA. Genome Med 7:76 |
| 26513549 | Zhou et al. (2015) TransVar. Nat Methods 12:1002-3 |
| 26553804 | O'Leary et al. (2016) RefSeq. Nucleic Acids Res 44:D733-45 |
| 26931183 | den Dunnen et al. (2016) HGVS 2016. Hum Mutat 37:564-9 |
| 27268795 | McLaren et al. (2016) VEP. Genome Biol 17:122 |
| 28967166 | Freeman et al. (2018) VariantValidator. Hum Mutat 39:61-68 |
| 29165669 | Landrum et al. (2018) ClinVar. Nucleic Acids Res 46:D1062-7 |
| 30129167 | Wang et al. (2018) biocommons/hgvs 2018. Hum Mutat 39:1803-13 |
| 33270643 | Hart & Prlić (2020) SeqRepo. PLoS ONE 15:e0239883 |
| 33538839 | Lefter et al. (2021) Mutalyzer 2. Bioinformatics 37:2811-7 |
| 34612497 | Park & Park (2022) Nomenclature variations. Lab Med 53:242-5 |
| 35357919 | Nurk et al. (2022) T2T genome. Science 376:44-53 |
| 35388217 | Morales et al. (2022) MANE. Nature 604:310-15 |
| 36124785 | Pozo et al. (2022) APPRIS/MANE. Bioinformatics 38:ii89-94 |
| 36318249 | Martin et al. (2023) Ensembl 2023. Nucleic Acids Res 51:D933-41 |
| 36420896 | Frankish et al. (2023) GENCODE 2023. Nucleic Acids Res 51:D942-9 |
| 36441169 | Wright et al. (2023) MANE clinical. Genet Med 25:100331 |
| 37165242 | Liao et al. (2023) Pangenome. Nature 617:312-24 |
| 39526381 | Goldfarb et al. (2025) RefSeq 25 years. Nucleic Acids Res 53:D243-57 |
| 39578691 | Landrum et al. (2025) ClinVar somatic. Nucleic Acids Res 53:D1313-21 |
