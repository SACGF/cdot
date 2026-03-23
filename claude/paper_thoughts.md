# cdot — Paper Framing Thoughts

## The Core Contribution

cdot solves a specific, painful infrastructure problem in clinical and research genomics: **getting transcript data into HGVS conversion tools, fast, without a database server, and with comprehensive historical coverage.**

The headline numbers frame the contribution cleanly:
- **Coverage:** 1.3M+ versioned transcript alignments vs ~141k in UTA uta_20210129 (9× increase)
- **Speed:** ~500–1000 transcripts/s local JSON vs ~1 transcript/s UTA remote (500–1000× faster)
- **Infrastructure:** zero-dependency (no PostgreSQL, works through firewalls, works offline)
- **Breadth:** GRCh37, GRCh38, T2T-CHM13v2.0; RefSeq + Ensembl; multi-species potential

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

## Question 1: Should cdot client code be contributed into biocommons/hgvs?

**Arguments for:**
- Makes cdot the official/canonical data provider for biocommons/hgvs
- Higher adoption — users of biocommons/hgvs get cdot "for free"
- Better long-term maintenance (community instead of SACGF alone)
- biocommons/hgvs#634 and #747 (canonical transcripts) suggest they want this

**Arguments against:**
- Loss of independent release cadence — cdot data schema and client would be coupled to biocommons/hgvs release cycle
- cdot serves PyHGVS too — the client is not biocommons-specific; a split would be needed
- The paper becomes a co-authorship / PR paper rather than a standalone contribution
- cdot's broader scope (data generation pipeline, REST service, T2T support) wouldn't fit inside biocommons/hgvs cleanly

**Recommended position:** Keep the client code in cdot, but explicitly position cdot as *the* recommended data provider for biocommons/hgvs, and make the integration seamless enough that biocommons could reference cdot in their own documentation. The paper can describe this as "a data provider for biocommons/hgvs," cite our contributions, and potentially include a biocommons author as a collaborator — without merging the repos.

After publication: contribute the data provider to biocommons/hgvs as an optional extra, keeping cdot as the source of truth for the data generation and REST service.

---

## Question 2: Should pyreference be in the same paper?

**Probably not, for these reasons:**
- Different audience: pyreference users are doing general bioinformatics (gene expression, splice analysis), not variant annotation
- Different story: the cdot paper is about *solving an HGVS infrastructure problem*; pyreference is about *using a well-structured annotation format*
- Dilutes the message of the cdot paper — reviewers would ask "what's the contribution, HGVS or general annotation?"
- pyreference likely has simpler functionality and wouldn't carry a separate paper on its own; a bioRxiv preprint is the right vehicle

**Recommended position:** Mention pyreference in the cdot paper as a one-paragraph demonstration that the JSON format has utility *beyond* HGVS (general bioinformatics, non-Python consumers, etc.). This strengthens the "format is generally useful" argument without making the paper unfocused. A separate pyreference bioRxiv preprint could follow.

---

## Question 3: What journal / scope to target?

**Target journal options:**
1. **Bioinformatics (Oxford)** — the natural home for tool papers. Application Notes (≤4 pages) or full papers. Application Note would be tight but doable; full paper is easier with benchmarks + case studies.
2. **Nucleic Acids Research** (Web Server / Database issue) — if framing as a database/resource paper. The Database issue specifically wants resources like this. Good fit.
3. **Genome Biology** — higher impact, more competitive; possible if we have strong benchmarks and case studies.
4. **GigaScience** — open data focus; good fit given the data release emphasis.
5. **F1000Research** — fast, open access; lower prestige but reaches the community.

**Recommendation:** Target **Nucleic Acids Research Database issue** primarily. The framing is: cdot is a database/resource — the JSON.gz files *are* the database, the REST API is the web server interface, and the Python clients are the programmatic interface. This fits NAR's Database issue criteria exactly. If rejected, Bioinformatics Applications Note is the fallback.

---

## Paper Angle and Structure

### Title (candidates)

- *cdot: A fast, comprehensive transcript database for HGVS variant annotation*
- *cdot: Compact transcript annotation for high-performance HGVS coordinate conversion*
- *Complete Dict Of Transcripts: A scalable resource for versioned human transcript coordinate data*

### Abstract framing

Problem: HGVS variant nomenclature is the clinical standard, but converting HGVS to/from genomic coordinates requires comprehensive, versioned transcript data. Existing resources (UTA) require database infrastructure, have limited transcript coverage, and lack Ensembl support.

Solution: cdot generates compact JSON files from RefSeq and Ensembl GTF/GFF3 annotation sources, covering 1.3M+ versioned transcript alignments across GRCh37, GRCh38, and T2T-CHM13v2.0. Data is accessible via local files (500–1000 transcripts/s) or a REST API, with no database infrastructure required.

Impact: Integration libraries for the two major Python HGVS tools (biocommons/hgvs, PyHGVS); resolves HGVS variants from ClinVar at X× the coverage and Y× the speed of UTA.

### Key sections

1. **Introduction**
   - HGVS nomenclature as clinical standard
   - The coordinate conversion problem
   - Limitations of UTA (PostgreSQL, firewall, no Ensembl, limited historical coverage)
   - Limitations of Counsyl PyHGVS (no alignment gap support, bugs)
   - cdot's approach: simple format, high coverage, no infrastructure

2. **Data sources and generation**
   - Ensembl GTF (releases 81–115, GRCh37/38/T2T)
   - RefSeq GFF3 (multiple annotation releases)
   - UTA CSV dumps (fills gaps in historical coverage)
   - JSON format design: shared fields vs build-specific, exon 6-tuple, gap encoding
   - Merge strategy and source priority

3. **Access methods**
   - Local JSON.gz (primary use case)
   - REST API (cdot.cc / cdotlib.org)
   - Python clients: biocommons/hgvs integration, PyHGVS integration

4. **Benchmarks**
   - Transcript coverage comparison (cdot vs UTA vs TARK)
   - Speed comparison (local JSON vs REST vs UTA PostgreSQL)
   - Accuracy validation (resolve ClinVar HGVS, compare to known genomic coordinates)
   - T2T-CHM13v2.0 novelty (first resource to support this for HGVS?)

5. **Case studies / use cases**
   - Clinical genomics: resolving historical HGVS from old clinical reports
   - MANE Select transcript lookup by gene name
   - Canonical transcript selection for gene-name HGVS

6. **Discussion**
   - Format extensibility (future: CDS phase, phasing array, exon dicts)
   - Non-HGVS uses (mention pyreference)
   - Comparison with TARK (Ensembl's own transcript archive)
   - Comparison with VariantValidator
   - Limitations: no protein sequence validation, CDS phase not yet modelled, Ensembl GFF3 protein version parsing

---

## Positioning vs Related Tools

| Tool | What it is | How cdot differs |
|------|-----------|-----------------|
| **UTA** | PostgreSQL transcript archive | No DB required; 9× more transcripts; includes Ensembl |
| **Ensembl TARK** | Ensembl transcript REST archive | Multi-consortium (RefSeq + Ensembl); offline capable; gap support |
| **biocommons/hgvs** | HGVS parsing/conversion library | cdot is a *data provider* for it, not a replacement |
| **Counsyl PyHGVS** | HGVS parsing library | cdot adds alignment gap support and data |
| **VariantValidator** | Web-based HGVS validator | cdot is embeddable/offline; no validation, just coordinate conversion |
| **VEP** | Variant annotation tool | VEP uses its own annotation; cdot is a transcript *source* not an annotation tool |
| **ANNOVAR** | Variant annotation tool | Same as VEP |
| **pyreference** | General annotation from cdot JSON | Non-HGVS use of same data format |

---

## The T2T Angle

Supporting Telomere-to-Telomere (T2T-CHM13v2.0) assembly for HGVS is potentially a first. If we can confirm no other tool supports HGVS resolution on T2T, this is a strong novelty claim. Worth checking against VariantValidator, VEP, and TARK.

---

## Author / Collaboration Strategy

- Core authors: Dave Lawrence, SACGF team
- Consider inviting: a biocommons/hgvs maintainer (acknowledges the integration work and broadens the community reach)
- Consider inviting: a clinical genomics collaborator who can provide a real-world case study (VariantGrid/Shariant usage)
- Acknowledgements: NCBI RefSeq team, Ensembl team, biocommons community

---

## Pre-submission checklist

- [ ] Finalise transcript count figures (exact, reproducible)
- [ ] Run formal ClinVar benchmark (#5)
- [ ] Write Snakemake pipeline summary statistics output
- [ ] Get cdot into CI with test suite (#62)
- [ ] Canonical transcript API working (#36) — strong use case for paper
- [ ] Domain migrated to cdotlib.org (#100)
- [ ] Confirm T2T novelty claim (check VariantValidator, TARK)
- [ ] Check pyreference paper status (include as 1 paragraph or separate ref)
