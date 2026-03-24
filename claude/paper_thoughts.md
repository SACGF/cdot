# cdot — Paper Framing Thoughts

## The Core Contribution

cdot solves a specific, painful infrastructure problem in clinical and research genomics: **getting transcript data into HGVS conversion tools, fast, without a database server, and with comprehensive historical coverage.**

The headline numbers frame the contribution cleanly:
- **Coverage:** 1.3M+ versioned transcript alignments vs ~141k in UTA uta_20210129 (9× increase)
- **Speed:** ~500–1000 transcripts/s local JSON vs ~1 transcript/s UTA remote (500–1000× faster)
- **Infrastructure:** zero-dependency (no PostgreSQL, works through firewalls, works offline)
- **Breadth:** GRCh37, GRCh38, T2T-CHM13v2.0; RefSeq + Ensembl; multi-species potential

### Why this matters — numbers from the literature

The clinical scale of the problem is large. ClinVar now holds >3 million variants from >2,800 organisations [Landrum 2025], virtually all expressed in HGVS transcript notation. Mutalyzer processed 133 million HGVS descriptions over five years, finding ~50% contained errors [Lefter 2021] — many attributable to missing or misidentified transcript data. McCarthy et al. [2014] showed only 44% of putative loss-of-function variants agreed between RefSeq and Ensembl annotation sets using ANNOVAR; Park et al. [2022] found ~85% agreement between ANNOVAR and SnpEff on clinical variants. Most dramatically, Münz et al. [2015] found that most annotation tools correctly annotated only 32% of BRCA2 variants — the rest failed due to improper handling of alignment gaps or strand orientation.

These numbers reframe the cdot contribution: this is not merely a *speed* problem (though cdot is 500–1000× faster than UTA remote access). It is fundamentally a **correctness and coverage** problem. Tools disagree, transcript versions are missing, and gapped alignments are handled incorrectly. cdot addresses all three.

---

## Project Ecosystem Map

```
cdot (JSON generation)          →  cdot_rest (REST delivery)
      ↓                                    ↓
  JSON.gz files            HTTP endpoints at cdot.cc / cdotlib.org
      ↓                                    ↓
biocommons/hgvs client ←——————————————————→  PyHGVS client
      ↓
  pyreference (non-HGVS bioinformatics)
```

### cdot (this repo)

Three separable components:
1. **Data generation** (`generate_transcript_data/`) — converts RefSeq/Ensembl GTF/GFF3 → JSON.gz
2. **HGVS client** (`cdot/hgvs/dataproviders/`) — biocommons HGVS data provider
3. **PyHGVS client** (`cdot/pyhgvs/`) — PyHGVS transcript factory

### cdot_rest (https://github.com/SACGF/cdot_rest)

REST API server that serves the JSON.gz data over HTTP. Users who don't want to download large files (72–122 MB each) can query transcripts on demand. Hosted at cdot.cc (migrating to cdotlib.org). This is the same data, different delivery mechanism — naturally the same paper.

### pyreference (https://github.com/SACGF/pyreference)

Uses cdot JSON for general bioinformatics (gene region queries, exon analysis) — not HGVS-specific. A separate use case that demonstrates the format's broader utility.

### biocommons/hgvs (https://github.com/biocommons/hgvs)

Community HGVS library that cdot integrates with as a data provider. cdot has made dozens of commits to this project. The relationship is: biocommons/hgvs is the engine; cdot is the fuel. Citing it extensively is appropriate.

---


---

## Question 3: What journal / scope to target?

**Target journal options:**
1. **Nucleic Acids Research** (Database/Web Server issue) — the JSON.gz files are a database, cdotlib.org is the web server, Python clients are the programmatic interface. Fits NAR's criteria exactly. Ensembl, RefSeq, GENCODE, and ClinVar all publish annual updates here. High visibility to the bioinformatics community.
2. **Genome Medicine** — McCarthy 2014 (transcript choice), Richards 2015 (ACMG), Münz 2015 (CAVA) all published here. Bridges the gap between computational tools and clinical genomics. Good fit for the clinical framing.
3. **Bioinformatics (Oxford)** — natural home for tool papers. Application Note (≤4 pages) is tight; full paper gives room for benchmarks.
4. **GigaScience** — open data focus; good fit given data release emphasis and T2T novelty.
5. **Human Mutation / Human Genetics and Genomics Advances** — Hart 2015 and Wang 2018 (biocommons/hgvs) and Freeman 2018 (VariantValidator) published here. Direct community relevance.

**Recommendation:** Target **Nucleic Acids Research Database issue** as first choice. If the scope feels too narrow for NAR reviewers (since cdot is software + data rather than purely a database), **Genome Medicine** is the strongest alternative — it has published the key papers in this space (McCarthy 2014, Münz 2015) and has strong clinical genomics readership.

**Human Mutation / HGGA** is worth considering specifically because it is where biocommons/hgvs and VariantValidator published — the paper would be read by exactly the right audience. The journal explicitly covers tools for variant handling.

---

## Paper Angle and Structure

### Title (candidates)

- *cdot: A fast, comprehensive transcript database for HGVS variant annotation*
- *cdot: Compact transcript annotation for high-performance HGVS coordinate conversion*
- *Complete Dict Of Transcripts: A scalable resource for versioned human transcript coordinate data*

### Abstract framing

Problem: HGVS variant nomenclature is the clinical standard [den Dunnen 2016] underpinning ACMG/AMP variant classification [Richards 2015] and ClinVar's >3 million variant archive [Landrum 2025]. Converting HGVS to/from genomic coordinates requires comprehensive, versioned transcript data, yet existing resources (UTA) require database infrastructure, have limited transcript coverage, and lack Ensembl support. Annotation tools drawing on different transcript sets disagree on up to 56% of loss-of-function variants [McCarthy 2014].

Solution: cdot generates compact JSON files from RefSeq and Ensembl GTF/GFF3 annotation sources, covering 1.3M+ versioned transcript alignments across GRCh37, GRCh38, and T2T-CHM13v2.0. Data is accessible via local files (500–1000 transcripts/s) or a REST API, with no database infrastructure required. Stores MANE Select, MANE Plus Clinical, and Ensembl canonical tags enabling programmatic canonical transcript selection.

Impact: Integration libraries for the two major Python HGVS tools (biocommons/hgvs, PyHGVS); resolves HGVS variants from ClinVar at X× the coverage and Y× the speed of UTA; first HGVS tool with T2T-CHM13v2.0 support.

### Key sections

1. **Introduction**
   - HGVS nomenclature as the clinical and research standard [den Dunnen 2016]
   - Scale of demand: ClinVar >3M variants [Landrum 2025], ACMG/AMP guidelines [Richards 2015], Mutalyzer 133M descriptions with 50% error rate [Lefter 2021]
   - The correctness problem: tool/transcript disagreement [McCarthy 2014, Park 2022]; 32% BRCA2 accuracy without gap support [Münz 2015]
   - The infrastructure problem: UTA requires PostgreSQL, is firewall-blocked, lacks Ensembl, covers only ~141k transcript versions
   - The MANE pressure: clinical labs now mandated to use MANE transcripts [Wright 2023] — tools must support canonical lookup
   - cdot's approach: unified, versioned, multi-source data accessible without infrastructure

2. **Data sources and generation**
   - Ensembl GTF (releases 81–115, GRCh37/38/T2T)
   - RefSeq GFF3 (multiple historical annotation releases — motivation: historical clinical HGVS strings)
   - UTA CSV dumps as baseline (fills pre-GFF gaps)
   - JSON format design: shared fields vs build-specific, exon 6-tuple, gap encoding (CIGAR)
   - Merge strategy and source priority order
   - MANE, RefSeq Select, Ensembl_canonical tag storage

3. **Access methods**
   - Local JSON.gz (primary use case)
   - REST API (cdot.cc / cdotlib.org) — analogous to SeqRepo REST [Hart 2020] for sequence data
   - Python clients: biocommons/hgvs integration, PyHGVS integration
   - FastaSeqFetcher as an alternative to SeqRepo for offline sequence access

4. **Benchmarks**
   - Transcript coverage comparison (cdot vs UTA vs TARK; include Ensembl-only count)
   - Speed: local JSON vs REST vs UTA PostgreSQL (informal: 500–1000× faster; formalise with script)
   - Accuracy: resolve ClinVar HGVS, compare to known VCF genomic coordinates
   - Coverage for historical strings: what fraction of ClinVar HGVS resolve with cdot vs UTA?
   - T2T-CHM13v2.0: show unique transcripts not in GRCh37/38

5. **Use cases**
   - Clinical genomics: resolving historical HGVS from old reports (VariantGrid/Shariant example)
   - Canonical transcript selection: `BRCA2:c.5946delT` resolved via MANE Select
   - Non-HGVS use: pyreference as a brief example that the JSON format is broadly useful

6. **Discussion**
   - cdot as the transcript coordinate layer in the biocommons stack (analogous to SeqRepo for sequences)
   - Comparison with TARK [Ensembl]: no RefSeq, online-only, no gap support for RefSeq transcripts
   - Comparison with VariantValidator [Freeman 2018]: monolithic service vs embeddable data layer
   - The pangenome future [Liao 2023]: cdot's multi-build architecture positions it for haplotype-level annotation
   - Format extensibility: CDS phase (#76), exon dicts (#78), VEP-compatible releases (#63)
   - Limitations: CDS phase/ribosomal frameshift not modelled; Ensembl GFF3 protein versions; no validation layer
   - Non-HGVS utility: pyreference demonstrates that the JSON format is a generally useful annotation resource

---

## Positioning vs Related Tools

| Tool | What it is | How cdot differs |
|------|-----------|-----------------|
| **UTA** | PostgreSQL transcript archive | No DB required; 9× more transcripts; includes Ensembl; works through firewalls |
| **Ensembl TARK** | Ensembl transcript REST archive | Multi-consortium (RefSeq + Ensembl); offline capable; RefSeq alignment gaps stored |
| **biocommons/hgvs** | HGVS parsing/conversion library | cdot is a *data provider* for it, not a replacement; analogous role to SeqRepo for sequences |
| **SeqRepo** | Local sequence storage | SeqRepo stores sequences; cdot stores transcript coordinates — complementary layers |
| **Counsyl PyHGVS** | HGVS parsing library | cdot adds alignment gap support and data; gap support fixes 32% → ~100% on BRCA2-class problems |
| **VariantValidator** | Web-based HGVS validator | cdot is embeddable/offline; no validation, just coordinate conversion; no server required |
| **Mutalyzer** | HGVS syntax checker | cdot provides data; Mutalyzer provides syntax validation — complementary |
| **TransVar** | Multilevel variant annotator | TransVar bundles own annotation; limited version management; cdot separates data from code |
| **VEP** | Variant annotation tool | VEP uses its own annotation; cdot is a transcript *source* not an annotation tool |
| **ANNOVAR / SnpEff** | Variant annotation tools | Bundle own transcripts; disagree on 15% of clinical variants [Park 2022]; cdot as shared layer would reduce this |
| **CSN / CAVA** | Clinical annotation tool | CAVA achieves 100% on BRCA2 by correct gap handling — same correctness goal cdot enables for biocommons/hgvs |
| **pyreference** | General annotation from cdot JSON | Non-HGVS use of same data format |

### The "shared data layer" angle

A secondary but powerful framing: if tools like VEP, ANNOVAR, TransVar, and biocommons/hgvs all drew transcript coordinate data from cdot, there would be one fewer source of annotation disagreement. The McCarthy 2014 and Park 2022 papers show tool disagreement is substantial; some of this stems from using different transcript versions. A shared, versioned, source-attributed data layer is a partial solution to this reproducibility problem.

---

## The T2T Angle

Supporting Telomere-to-Telomere (T2T-CHM13v2.0) assembly for HGVS is potentially a first. The T2T genome [Nurk 2022] adds ~200M bp of previously inaccessible sequence and corrects errors in GRCh38, including regions that contain real genes. Any HGVS string for a transcript annotated on T2T cannot be resolved by existing tools. If we can confirm no other HGVS coordinate conversion library supports T2T (check: VariantValidator, VEP, TARK), this is a strong novelty claim.

## The MANE Canonical Transcript Angle

This is more urgent than initially appreciated. The MANE Nature paper [Morales 2022] and the Wright 2023 editorial in Genetics in Medicine show that clinical labs are actively being called on to standardise on MANE transcripts. cdot already stores MANE tags (since v0.2.12). A `CanonicalTranscriptSelector` API (issue #36) that exposes these tags would make cdot the first HGVS data provider to offer programmatic canonical transcript lookup. This is a paper-ready feature that differentiates cdot from all other HGVS data sources and directly addresses a stated clinical need.

**This should be implemented before submission and featured prominently** — it turns cdot from a passive data source into an active contributor to the MANE adoption problem.

## The SeqRepo Analogy

The literature review clarifies a clean positioning statement: **cdot is to transcript coordinates what SeqRepo is to biological sequences**. Both projects provide local, high-performance, offline-capable access to data that was previously only available from slow remote databases. SeqRepo achieves 1300× speedup over remote access [Hart 2020]; cdot achieves 500–1000× speedup over UTA remote. Both are part of the biocommons stack. Framing cdot explicitly this way in the paper will resonate with the biocommons community and make the contribution concrete for reviewers unfamiliar with the HGVS pipeline.

---

## Scattered thoughts

JSON has become the standard for data transfer. Almost all libraries can rapidly load it. Useful not just for Python
Mention pyreference in the cdot paper as a one-paragraph demonstration that the JSON format has utility *beyond* HGVS (general bioinformatics, non-Python consumers, etc.). This strengthens the "format is generally useful" argument without making the paper unfocused. A separate pyreference bioRxiv preprint could follow.

## Author / Collaboration Strategy

- Core authors: Dave Lawrence, Shariant team who helped debug HGVS issues
- Consider inviting: a biocommons/hgvs maintainer (acknowledges the integration work and broadens the community reach)
- Consider inviting: a clinical genomics collaborator who can provide a real-world case study (VariantGrid/Shariant usage)
- Acknowledgements: NCBI RefSeq team, Ensembl team, biocommons community

---

## Pre-submission checklist

**Data / benchmarks**
- [ ] Finalise transcript count figures (exact, reproducible from Snakemake summary stats)
- [ ] Run formal ClinVar benchmark (#5) — resolve all ClinVar HGVS, report coverage and accuracy vs UTA
- [ ] Speed benchmark: local JSON vs REST vs UTA PostgreSQL (scripted, reproducible)
- [ ] T2T unique transcript count (not in GRCh37/GRCh38)

**Features that must exist before submission**
- [ ] Canonical transcript API working (#36) — `CanonicalTranscriptSelector` with MANE + RefSeq Select + Ensembl canonical; needed as paper use case
- [ ] Domain migrated to cdotlib.org (#100) — URL in paper must be stable
- [ ] Snakemake pipeline summary statistics output (transcript counts by source, genome build, biotype)

**Correctness / infrastructure**
- [ ] Get cdot into CI with test suite (#62)
- [ ] Python version bump + build dependencies (#104) — cannot publish with broken setup.cfg

**Verification**
- [ ] Confirm T2T novelty claim — check VariantValidator, VEP, TARK, TransVar for T2T support
- [ ] Confirm BRCA2 accuracy improvement is attributable to gap support (run with/without FastaSeqFetcher or gapped alignment handling)
- [ ] Check pyreference paper status — include as 1 paragraph in Discussion or cite as companion resource
