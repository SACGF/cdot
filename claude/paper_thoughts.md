# cdot — Paper Framing Thoughts

## The Mission Statement

> **cdot's goal is to resolve as many HGVS strings from the wild as possible.**

A HGVS string from a clinical report, a search box, or ClinVar can fail to resolve for several distinct reasons. cdot attacks each one:

| Failure mode | cdot's answer |
|---|---|
| String is malformatted | `clean_hgvs()` — fixes whitespace, casing, missing colons, swapped gene/transcript |
| Transcript version not in database | 1.3M+ versioned alignments across all historical RefSeq + Ensembl releases |
| Exact version missing, adjacent exists | `get_best_transcript_version()` — controlled fallback to nearest version |
| Transcript only in one consortium | Both RefSeq and Ensembl; clinical reports use both |
| Transcript on wrong build | GRCh37, GRCh38, T2T-CHM13v2.0 all supported |
| Alignment gap not handled | CIGAR gap encoding; fixes the 32% BRCA2 accuracy problem [Münz 2015] |
| Gene-only HGVS with no transcript | MANE/canonical tags enable programmatic transcript selection |
| No database server available | Zero-dependency local JSON.gz; works offline, through firewalls |
| Lookup too slow for bulk use | 500–1000× faster than UTA remote |

This framing unifies what might otherwise look like a miscellaneous collection of features. The cleaner, the transcript archive, the version fallback, the gap support, canonical lookup, and infrastructure choices are all expressions of the same goal.

---

## Title (candidates)

- *cdot: maximising HGVS resolution across clinical and research genomics*
- *cdot: a comprehensive transcript resource for resolving HGVS variant nomenclature*
- *cdot: Complete Dict Of Transcripts for high-coverage HGVS resolution*

The title should contain "HGVS" and signal the breadth/resolution angle, not just "fast" or "database".

---

## Abstract

**Problem:** HGVS variant nomenclature is the clinical standard [den Dunnen 2016] underpinning ACMG/AMP classification [Richards 2015] and ClinVar's >3 million variants [Landrum 2025]. Yet a large fraction of HGVS strings in the wild fail to resolve — due to formatting errors, missing transcript versions, absent Ensembl support, alignment gaps, or database infrastructure barriers. Mutalyzer found ~50% of 133 million submitted HGVS descriptions contained errors [Lefter 2021]; many tools correctly annotate only 32% of transcripts with alignment gaps [Münz 2015]; and the primary transcript archive (UTA) covers only ~141k versioned accessions and requires PostgreSQL.

**Solution:** cdot is a toolkit for maximising HGVS resolution. It provides: (1) `clean_hgvs()`, a pure-string preprocessor that fixes the most common real-world formatting mistakes derived from analysis of clinical search logs; (2) a compact JSON transcript archive covering 1.3M+ versioned alignments from RefSeq and Ensembl across GRCh37, GRCh38, and T2T-CHM13v2.0, with full CIGAR gap encoding; (3) `get_best_transcript_version()` for controlled version fallback; and (4) MANE Select / Ensembl canonical tags for resolving gene-only HGVS. Data is accessible via local JSON.gz (500–1000 transcripts/s, no infrastructure required) or a public REST API.

**Impact:** Resolves X% more ClinVar HGVS variants than UTA; [N]× faster for bulk processing; first HGVS tool with T2T-CHM13v2.0 support; integrated with the two major Python HGVS libraries (biocommons/hgvs, PyHGVS).

---

## Paper Structure

### 1. Introduction

Frame the problem as a **resolution rate** problem, not a speed or infrastructure problem:

- HGVS as the clinical standard [den Dunnen 2016]; scale: ClinVar >3M variants [Landrum 2025], ACMG/AMP [Richards 2015]
- The resolution problem has multiple independent causes:
  - Formatting errors in the wild (Mutalyzer: 50% error rate [Lefter 2021]; our search-log analysis shows the specific error types)
  - Missing transcript versions: UTA covers ~141k versioned accessions; clinical reports use historical versions
  - Annotation disagreement: tool/transcript divergence [McCarthy 2014, Park 2022]
  - Alignment gaps: 32% BRCA2 accuracy without gap support [Münz 2015]
  - Infrastructure: UTA requires PostgreSQL, is firewall-blocked, lacks Ensembl
  - MANE transition: labs must standardise on MANE transcripts [Wright 2023] — gene-only HGVS must be resolvable
- cdot's approach: address all failure modes in one toolkit

### 2. Methods — structured around the failure modes

#### 2.1 String cleaning (`cdot/hgvs/clean.py`)

`clean_hgvs()` is a pure-string preprocessing step requiring no data provider. Derived from analysis of real clinical search logs (Shariant, VariantGrid). Fix categories with real examples from the logs:

| Fix | Example |
|-----|---------|
| Whitespace stripping | `"NM_205768 c.44A>G"` → `"NM_205768:c.44A>G"` |
| Protein suffix removal | `"NM_000059.4:c.316G>A p.Arg106*"` → `"NM_000059.4:c.316G>A"` |
| Double colon | `"NM_004245: :c.337G>T"` → `"NM_004245:c.337G>T"` |
| Mutation type case | `"c.1126_1133DUP"` → `"c.1126_1133dup"` |
| Unbalanced brackets | `"NM_001754.5):c.557T>A"` → `"NM_001754.5:c.557T>A"` |
| Missing `c.`/`g.` | `"NM_001754.5:557T>A"` → `"NM_001754.5:c.557T>A"` |
| Swapped gene/transcript | `"BRCA1(NM_000059.4):c.316G>A"` → `"NM_000059.4(BRCA1):c.316G>A"` |
| Lowercase transcript | `"nm_000059.4:c.316G>A"` → `"NM_000059.4:c.316G>A"` |

Pull actual percentages by error type from `analysis/HGVS cleaning.ipynb` — this is the key figure.

#### 2.2 Transcript archive (coverage)

- **Sources:** Ensembl GTF releases 81–115 (GRCh37/38/T2T); RefSeq GFF3 multiple historical releases; UTA CSV as baseline
- **Why historical releases:** a clinical HGVS from 2015 may use a transcript version no longer in the current RefSeq annotation; cdot keeps every version from every release
- **Format:** compact JSON.gz; exon 6-tuple with CIGAR gap encoding; shared vs build-specific fields
- **Merge strategy:** UTA as baseline, annotation releases override chronologically (newer wins)
- **Coverage:** 1.3M+ versioned alignments vs ~141k in UTA uta_20210129 (9×)

#### 2.3 Transcript version fallback (`get_best_transcript_version()`)

When the exact requested version is absent, fall back to the nearest available version. Three strategies: UP_THEN_DOWN (prefer higher versions first), CLOSEST (alternate ±1), LATEST. This is only meaningful when cdot's archive has adjacent versions — which UTA often cannot provide given its 9× coverage gap.

#### 2.4 Alignment gap support

Transcripts with insertion/deletion gaps relative to the genome (RefSeq `cDNA_match` features) are encoded in GFF3 gap format (`M196 I1 M61`) and converted to HGVS CIGAR on read. Without this, tools fail on a substantial fraction of clinically relevant transcripts — Münz et al. [2015] showed 32% accuracy on BRCA2-class variants without gap support.

#### 2.5 Multi-build and multi-source coverage

- Both RefSeq and Ensembl: clinical labs use both; some transcripts exist only in one consortium
- GRCh37 + GRCh38: older clinical reports are on GRCh37
- T2T-CHM13v2.0: first HGVS tool to support this assembly (verify claim against VariantValidator, VEP, TARK)

#### 2.6 Canonical transcript selection (MANE)

Clinical labs are moving to MANE transcripts [Wright 2023, Morales 2022]. cdot stores MANE Select, MANE Plus Clinical, RefSeq Select, and Ensembl canonical tags. A `CanonicalTranscriptSelector` API (issue #36) enables resolution of gene-only HGVS (e.g. `BRCA2:c.5946delT`) by selecting the appropriate canonical transcript — a failure mode that no other HGVS data provider currently addresses programmatically.

**Must be implemented before submission** — it's a paper-ready feature.

#### 2.7 Access methods

- **Local JSON.gz:** primary use case; 500–1000 transcripts/s; no infrastructure; works offline
- **REST API (cdotlib.org):** for users who cannot download 72–122 MB files; analogous to SeqRepo REST [Hart 2020] for sequences
- **Client libraries:** biocommons/hgvs integration (`JSONDataProvider`); PyHGVS integration (`JSONPyHGVSTranscriptFactory`)
- **FastaSeqFetcher:** offline sequence access from genome FASTA, handling gapped alignments

### 3. Results / Benchmarks

All benchmarks framed as: **how many more HGVS strings can be resolved, and how fast?**

- **Resolution rate:** resolve all ClinVar HGVS strings; report % resolved with cdot vs UTA (the primary figure)
- **Coverage:** transcript count by source, build, biotype (from Snakemake summary stats)
- **Speed:** local JSON vs REST vs UTA PostgreSQL (scripted, reproducible; informal target: 500–1000×)
- **Cleaning impact:** from search-log analysis, % of failed searches that `clean_hgvs()` would fix
- **Version fallback impact:** % of "version not found" failures where an adjacent version exists in cdot
- **Gap support:** reproduce the BRCA2 accuracy comparison with/without gap support
- **T2T:** unique transcript count not in GRCh37/38; example HGVS resolvable only via T2T

### 4. Use Cases

Structured as "classes of HGVS that previously couldn't resolve":

1. **Malformatted clinical reports:** `clean_hgvs()` preprocessing pipeline; show before/after on real examples from the search log
2. **Historical HGVS:** old transcript version in a 2010 clinical report → resolved via cdot's historical archive
3. **Gene-only HGVS:** `BRCA2:c.5946delT` → resolved via MANE Select canonical transcript
4. **T2T variants:** transcript only annotated on T2T; unresolvable by any other tool
5. **Bulk ClinVar resolution:** speed benchmark in context of real workload

### 5. Discussion

- **The resolution rate lens:** cdot is best understood not as "a database" but as a resolution-maximising toolkit; each component addresses an independent failure mode
- **SeqRepo analogy:** cdot is to transcript coordinates what SeqRepo is to sequences — local, fast, offline-capable, part of the biocommons stack [Hart 2020]
- **The shared data layer argument:** if tools like VEP, ANNOVAR, SnpEff all drew from a shared versioned transcript source, tool disagreement [McCarthy 2014, Park 2022] would be reduced; cdot's format is a candidate for this role
- **Comparison with TARK:** no RefSeq, online-only, no gap support
- **Comparison with VariantValidator [Freeman 2018]:** monolithic service vs embeddable data layer; cdot has no validation, only coordinate conversion
- **Pangenome future [Liao 2023]:** multi-build architecture positions cdot for haplotype-level annotation
- **Limitations:** CDS phase/ribosomal frameshift not modelled (#76); Ensembl GFF3 protein versions broken (#101); no validation layer (by design — validation is Mutalyzer/VariantValidator's job)
- **Non-HGVS utility:** pyreference demonstrates the JSON format is broadly useful; one paragraph only

---

## Positioning vs Related Tools

| Tool | What it is | How cdot differs |
|------|-----------|-----------------|
| **UTA** | PostgreSQL transcript archive | No DB required; 9× more transcripts; includes Ensembl; gap support; works offline |
| **Ensembl TARK** | Ensembl transcript REST archive | RefSeq + Ensembl; offline; RefSeq gaps stored; T2T support |
| **biocommons/hgvs** | HGVS parsing/conversion library | cdot is a *data provider* for it — the fuel, not the engine |
| **SeqRepo** | Local sequence storage | Complementary: SeqRepo = sequences; cdot = transcript coordinates |
| **Counsyl PyHGVS** | HGVS parsing library | cdot adds gap support + 9× more data; fixes 32% → ~100% on BRCA2-class |
| **VariantValidator** | Web-based HGVS validator | cdot is offline/embeddable; no validation, just coordinate conversion |
| **Mutalyzer** | HGVS syntax checker | Complementary: Mutalyzer validates syntax; cdot provides data |
| **TransVar** | Multilevel variant annotator | TransVar bundles own annotation; no version management; cdot separates data from code |
| **VEP / ANNOVAR / SnpEff** | Variant annotation tools | Bundle own transcripts; disagree on 15% of clinical variants [Park 2022]; cdot as shared layer would reduce this |
| **CAVA / CSN** | Clinical annotation tool | CAVA achieves 100% on BRCA2 by correct gap handling — same correctness goal cdot enables for biocommons/hgvs |

---

## Journal Target

The primary contribution of cdot is a **data generation pipeline + library integration toolkit**, not a web server or database. The REST service is a convenience wrapper around the JSON files. This matters for journal fit.

**First choice: Bioinformatics (Oxford)** — cdot's users are developers integrating it into pipelines and tools; this is exactly the Bioinformatics readership. Full article gives room for benchmarks; Application Note (≤4 pages) is too tight given the breadth of the contribution.

**Second choice: Human Mutation / Human Genetics and Genomics Advances** — where biocommons/hgvs [Wang 2018] and VariantValidator [Freeman 2018] published. Strong domain fit; reviewers will understand the problem without background. More clinical audience than cdot's primary users.

**Third choice: Genome Medicine** — McCarthy 2014, Münz 2015 published here; bridges computational and clinical. Good fit for the resolution-rate framing.

**NAR is a weak fit.** NAR Database issue is for curated biological databases (Ensembl, ClinVar, UniProt); the JSON.gz files are a derived resource, not an independently curated database. NAR Web Server issue requires the web server to be central — but cdot_rest is a minor convenience component, not the contribution. Submitting to NAR risks rejection on scope grounds.

---

## Project Ecosystem Map

```
cdot (JSON generation)          →  cdot_rest (REST delivery)
      ↓                                    ↓
  JSON.gz files            HTTP endpoints at cdotlib.org
      ↓                                    ↓
biocommons/hgvs client ←——————————————————→  PyHGVS client
      ↓
  pyreference (non-HGVS bioinformatics)
```

- `cdot/hgvs/clean.py` — string preprocessing; no data provider needed; ships with cdot
- `cdot/hgvs/dataproviders/` — biocommons HGVS integration (`JSONDataProvider`, `RESTDataProvider`)
- `cdot/pyhgvs/` — PyHGVS integration (legacy; PyHGVS abandoned)
- `generate_transcript_data/` — data pipeline (not pip-distributed)

---

## Author / Collaboration Strategy

- Core authors: Dave Lawrence, Shariant team who helped debug HGVS issues
- Consider inviting: a biocommons/hgvs maintainer (acknowledges integration work; broadens community reach)
- Consider inviting: a clinical genomics collaborator for the real-world case study
- Acknowledgements: NCBI RefSeq team, Ensembl team, biocommons community

---

## Pre-submission Checklist

**Data / benchmarks**
- [ ] Finalise transcript count figures (exact, reproducible from Snakemake summary stats)
- [ ] Run formal ClinVar benchmark — resolve all ClinVar HGVS, report % coverage vs UTA (the primary figure)
- [ ] Speed benchmark: local JSON vs REST vs UTA PostgreSQL (scripted, reproducible)
- [ ] T2T unique transcript count (not in GRCh37/GRCh38)
- [ ] Search-log analysis: pull per-error-type percentages from `analysis/HGVS cleaning.ipynb`
- [ ] Version fallback analysis: what % of "version not found" cases have an adjacent version in cdot?

**Features that must exist before submission**
- [ ] Canonical transcript API (#36) — `CanonicalTranscriptSelector` with MANE + RefSeq Select + Ensembl canonical
- [ ] Domain migrated to cdotlib.org (#100) — URL in paper must be stable
- [ ] Snakemake pipeline summary statistics (transcript counts by source, build, biotype)

**Correctness / infrastructure**
- [ ] cdot into CI with test suite (#62)
- [ ] Python version bump + build dependencies (#104)

**Verification**
- [ ] Confirm T2T novelty claim — check VariantValidator, VEP, TARK, TransVar for T2T support
- [ ] Confirm BRCA2 accuracy improvement is attributable to gap support
- [ ] Check pyreference paper status — cite as companion or include as one Discussion paragraph
