# cdot Paper Plan

Working plan for the cdot manuscript. Supersedes the scope decisions in `paper/notes.md`
(which assumed an Application Note). Sits alongside `claude/paper_thoughts.md` (framing) and
`claude/benchmark_plan.md` (benchmark methodology). Section drafts live in `paper/*.md`.

Last updated: 2026-06-17.

---

## 1. Decisions locked

| Decision | Choice | Consequence |
|---|---|---|
| **Article type** | **Full Original Article** (Bioinformatics, or Human Mutation / GHGA Advances) | Room for a real Results section: coverage + cleaning + taxonomy + speed + T2T. Resolves the AppNote-vs-full tension in `notes.md`. |
| **Private-corpus use** | **Hybrid** — reproducible injection benchmark is primary; real-corpus stats are aggregate, non-reproducible *production validation* | Two-tier fact model (§3). Referee can rerun the headline cleaning claim. |
| **Cleaning emphasis** | **Co-equal contribution** — "repair, don't just validate" + the real-world HGVS entry-error taxonomy as headline novelty | `clean_hgvs()` gets its own Methods + Results treatment (it is absent from the current `implementation.md` draft — must be added). Appears in title/abstract. |

## 2. Decisions resolved (2026-06-17 chat)

- **Cleaning webapp → deferred to a SECOND paper (NAR Web Server).** This paper covers the cleaning
  *library* (`clean_hgvs()`) + the reproducible injection benchmark + the taxonomy. A dedicated
  interactive cdot **cleaner web service** ("paste messy HGVS → repaired + fix-codes + resolved")
  becomes paper #2, targeting **NAR Web Server** — a strong fit for a real server, where cdot-the-
  resource is a weak fit. Removes the webapp from this paper's critical path entirely. When built,
  spec it **stateless / no query logging** (or aggregate counts only) so it does not recreate the
  private-corpus governance problem.
- **Governance — CONFIRMED.** Aggregate statistics from the production query logs may be published.
  Needs a data-availability + ethics sentence stating aggregate-only, underlying queries not
  shareable. Corpus described by source (see below).
- **Journal — DECIDED: Bioinformatics (Oxford), Original Paper.** (Not the 4-page Application Note —
  full article for the coverage + cleaning + taxonomy + benchmark breadth.) Framing stays
  "developer tool + data resource". See §8 for rationale and runners-up.

### Corpus sources (for the Methods corpus description — aggregate use only)

Production search-bar HGVS queries pooled from SA Pathology / Australian clinical & research
variant platforms (CSV sources in `cdot_private/csv/`):

| CSV | Platform |
|---|---|
| `shariant_search_hgvs.csv` | Shariant — national variant-classification sharing platform |
| `vgaws_search_hgvs.csv` | VariantGrid (AWS-hosted instance) |
| `vg3upgrade2_search_hgvs.csv` | **SA Pathology VariantGrid clinical server** |
| `runx1db4_search_hgvs.csv` | RUNX1 database (research/curation) |

Combined: N = 32,752 queries (`combined_search_hgvs.csv`). Describe generically as "production
clinical and research variant-curation platforms" in the paper; name them in acknowledgements/data
statement per governance. Report **aggregates only** — no individual query strings.

### Still open
- [ ] **Author list / corresponding author** (carried from `notes.md`).
- [ ] Optional: second-rater pass on the taxonomy for a Cohen's κ reliability figure (§6).

## 3. The two-tier fact model (the key privacy mechanism)

The paper's fact system (`output/facts/*.json` → `{{ }}` templating) is for **reviewer-reproducible**
numbers. Private-corpus numbers must NOT masquerade as reproducible facts. Split:

**Tier 1 — public, reproducible (`output/facts/`, normal `{{ }}` facts):**
- ClinVar resolution (correct / no_data / incorrect), cdot vs UTA — `benchmark_resolution.py`.
- **Injection-cleaning recovery** — perturb public ClinVar c.HGVS with each `clean.py` fix
  category at frequencies *informed by* the aggregate real distribution; report recovery % per
  class + **0 regressions**. New script (§5). This is the headline cleaning claim.
- Version-bump recovery + false-rescue rate (`get_best_transcript_version()`).
- Speed (load, throughput, prefetch) — already in `output/facts/resolution_*.json`.
- Coverage counts (Snakemake summary stats).

**Tier 2 — private corpus, aggregate-only, NON-reproducible (frozen constants, flagged):**
- Real rescue rate: **91.5% → 96.4% (+4.9%), 0 regressions**, N = 32,752 production queries.
- Rescue-op distribution (which fixes do the rescuing — whitespace, case, swapped gene/tx, …).
- Residual **8-class error taxonomy** distribution (TRUNCATED 24.2%, MISSING_REFERENCE 23.6%, …).
- Stored in `cdot_private/output/` (committed there: `residual_classification_metadata.json`,
  `residual_classification_summary.txt`, `METHODS_residual_classification.md`). **Never** copy the
  strings into this repo; examples in the paper are *synthesised* from public ClinVar transcripts
  (NM_000059.4 BRCA2, NM_001754.5 RUNX1).

Mechanism: Tier 2 numbers go into the manuscript as literal constants with an explicit
"production validation; underlying queries are not shareable for privacy reasons" caveat — they do
**not** get a regenerable `output/facts/` entry. Consider a `paper/private_facts.md` (or a clearly
named constants block) so it is obvious which numbers a referee cannot reproduce.

## 4. Narrative spine and section map

Thesis (from `paper_thoughts.md`): **cdot maximises resolution of real-world HGVS** by attacking
each independent failure mode. Cleaning is now a first-class pillar, not a footnote.

**File structure (LOCKED + WIRED 2026-06-17).** `paper.toml`, `paper/Snakefile`, `paper/README.md`
updated to the Original-Paper layout. `implementation.md` → renamed `methods.md`; new `results.md`
stub added. `paper/notes.md` carries a SUPERSEDED banner.

| Section | File | Current state | Plan |
|---|---|---|---|
| **Title** | (abstract) | candidates in `paper_thoughts.md` | Pick one signalling resolution + repair (e.g. "cdot: maximising real-world HGVS resolution"). |
| **Abstract** | `abstract.md` | drafted (resolution framing present) | Add one cleaning sentence (repair + 0 regressions) and the taxonomy hook. |
| **Introduction** | `introduction.md` | drafted, 3 paras | Keep. Strengthen the "errors in the wild" para with the taxonomy angle (Lefter reported a *rate*; we give a repair-oriented *taxonomy* + a fix). |
| **Methods** | `methods.md` | renamed from `implementation.md`; has restructure-TODO banner | **Add String-cleaning subsection** (`clean_hgvs()` — pure-string, ordered fix categories, no data provider). **Move benchmark/coverage prose OUT to results.md.** Keep data/format/access/canonical. |
| **Results** | `results.md` | **STUB created** (R1–R6 skeleton) | Write prose: (R1) coverage + ClinVar resolution [Fig 1]; (R2) cleaning — injection recovery [Tier 1] + real-corpus validation [Tier 2]; (R3) residual taxonomy [Fig 2]; (R4) version-bump recovery; (R5) speed; (R6) gap support (Münz) + T2T uniqueness. |
| **Discussion** | `discussion.md` | drafted | Add: repair-vs-validate positioning vs Mutalyzer/VariantValidator; the taxonomy as a reusable artifact; the *ceiling of cleaning* (residual ~48% irreducible user input, 6.9% biocommons grammar gap, 6.9% non-HGVS); point to the planned NAR Web Server cleaner (Paper 2). |
| **Figures** | `figures.md` | Fig 1 (coverage+ClinVar), Fig 2 (arch, optional) | Repurpose **Fig 2 = cleaning**: bar of injection recovery per class + the residual taxonomy donut/stacked bar. Arch diagram → supplementary. |
| **Supplementary** | `supplementary.md` | tables S1–S4 | Add S5 injection-recovery-per-class table; S6 residual taxonomy table (aggregate); S7 rescue-op distribution. |

## 5. New / changed work items

> Done & removed: `analysis/inject_and_clean.py` (Tier-1 injection benchmark script — exists);
> String-cleaning Methods subsection (written in `paper/methods.md`).

**Scripts & data (this repo, `analysis/`):**
1. [ ] Wire `inject_and_clean.py` output into `output/facts/` so Fig 2 / Table S5 regenerate from facts.
2. [ ] Version-bump sub-benchmark (drop requested version → `get_best_transcript_version()` →
   rescue % + false-rescue %). (`benchmark_plan.md` Table 3.2.)
3. [ ] Larger ClinVar resolution run (1k/10k via VEP recipe) so resolution numbers are not just the
   "easy" 100-set (`benchmark_plan.md` §6a caveat). Keep seed + accession list committed.
4. [ ] Taxonomy figure data: aggregate counts already in
   `cdot_private/output/residual_classification_summary.txt` — copy as *numbers only* into a
   Tier-2 constants block; synthesise example strings from public transcripts.

**Manuscript (`paper/`):**
5. [ ] Promote Results to its own section; weave Tier 1 vs Tier 2 with explicit reproducibility flags.
6. [ ] `paper/private_facts.md` (or constants block) for the Tier-2 numbers + caveat text.
7. [ ] Ethics/data-availability statement covering the aggregate search-log use.
8. [ ] Update `figures.md` (Fig 2 = cleaning) and `supplementary.md` (S5–S7).

**Optional (deferred, post-benchmark-lock):**
9. [ ] `cdot_rest` `POST /clean` endpoint + minimal demo UI (stateless, no query logging).

**Carried from existing checklists (`notes.md`, `paper_thoughts.md`):** #55 FastaSeqFetcher fix
(gap-correctness claim), Zenodo DOI, confirm T2T-novelty + Münz gap attribution.
*(Done: #36 gene/canonical resolution, #100 cdotlib.org live, #62 CI.)*

## 6. Methodology notes to state explicitly (referee-proofing)

- Cleaning recovery headline is **reproducible** (injection on public ClinVar); real-corpus rate is
  **production validation** (aggregate, not reproducible) — state which is which, do not blur.
- Error taxonomy was **LLM single-label, single-rater** (8-class decision-tree; see
  `cdot_private/output/METHODS_residual_classification.md`). Disclose the single-rater limit; if a
  reliability figure is wanted, run a second independent pass and report Cohen's κ before submission.
- "incorrect" buckets need a manual audit (UTA-vs-annotation alignment diffs, #95) — some are a
  *paragraph*, not a bug. Normalise before string-compare.
- Log every silent cap (sample size, seed, skipped buckets).

## 7. Immediate next actions (proposed order)

1. Confirm the open decisions in §2 (webapp scope, governance, journal).
2. Add the String-cleaning Methods subsection + restructure Results into its own section (drafts only).
3. Build `analysis/inject_and_clean.py` and produce the Tier-1 cleaning facts.
4. Draft Fig 2 (cleaning recovery + residual taxonomy) and the Tier-2 constants block.
5. Scale the ClinVar resolution run (1k) for non-trivial resolution numbers.

## 8. Journal recommendation (full Original Article)

Two-paper strategy: **Paper 1 (this one)** = resource + library + cleaning + taxonomy + benchmarks;
**Paper 2 (later)** = interactive cdot cleaner web service → **NAR Web Server**.

Ranked for Paper 1:

1. **Bioinformatics (Oxford) — Original Paper. [DECIDED — 2026-06-17.]** cdot's users are developers
   integrating it into pipelines — exactly the readership. biocommons/hgvs (Hart 2015), the engine
   cdot fuels, published here. Full Original Paper (not the 4-page Application Note) gives room for
   the Results section. Strong methods/reproducibility expectations — our Tier-1/Tier-2 split and
   committed scripts fit well.
2. **Human Genetics and Genomics Advances (HGGA, ASHG/Cell Press, OA).** The clinical-variant
   domain home: VariantValidator (Freeman 2018) and biocommons/hgvs (Wang 2018) appeared in
   *Human Mutation*, which Wiley discontinued (end 2022); HGGA is its de-facto successor venue.
   Reviewers feel the HGVS-resolution pain directly. Best if we lead with clinical impact.
   *(Verify Human Mutation status / current best successor before committing.)*
3. **GigaScience.** Best if we foreground the *reproducible data resource + Snakemake pipeline +
   Zenodo deposition + benchmark* angle; values data availability and reproducibility highly.
4. **Genome Medicine.** McCarthy 2014 and Münz 2015 (both cited) published here; bridges
   computational and clinical. Favours a clinical-impact narrative over a tool/resource one.

Reserved: **NAR Web Server** → Paper 2 (cleaner web service). **NAR Database** remains a weak fit
(derived resource, not a curated database) — do not target.

Decision driver: if the framing stays "developer tool + data resource" → Bioinformatics; if it
shifts to "clinical-variant resolution impact" → HGGA / Genome Medicine; if "reproducible resource"
→ GigaScience.
