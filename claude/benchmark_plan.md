# cdot Benchmark Plan

Benchmarks for the paper. The framing (from `claude/paper_thoughts.md`) is **resolution
rate**, not just speed: *how many HGVS strings from the wild can cdot resolve, how many can it
clean/fix, and how fast?* Speed is necessary but secondary — the headline figure is coverage.

## 0. The three questions to answer

| Axis | Question | Primary metric |
|------|----------|----------------|
| **Resolution** | Of real ClinVar c.HGVS, how many can we convert to the correct genomic coordinate? | % correct vs ground-truth g.HGVS, cdot vs UTA |
| **Recovery (clean/fix)** | Of strings that *fail* (malformed, wrong version), how many does `clean_hgvs()` + version-bumping rescue? | % of failures recovered |
| **Speed** | How fast for bulk processing? | transcripts/s and HGVS/s; local JSON vs REST vs UTA |

Each axis maps to a paper Results subsection and to an open issue: resolution → #5, recovery →
#27/#28 (version bumping) + `clean.py`, speed → #5.

## 1. Ground truth & datasets

**Source of truth:** ClinVar VCF. `CLNHGVS` gives the normalised genomic g.HGVS
(`NC_000001.11:g.12007126G>A`). That is the gold standard a c→g conversion must reproduce.

**Getting the c.HGVS** (ClinVar VCF has no transcript HGVS):
- **RefSeq:** VEP annotate with `--pick --hgvs --refseq` → `HGVSc` from the CSQ field. (Recipe
  already at the bottom of `tests/benchmark_hgvs.py`.)
- **Ensembl:** import the ClinVar subset into VariantGrid, populate ClinGen Alleles, pull the
  MANE Ensembl HGVS. (Recipe also in `tests/benchmark_hgvs.py`.)

**Already extracted (use for dev / small runs — no VEP needed):**
`tests/test_data/clinvar_hgvs/*.tsv` — tab-separated `g.HGVS<TAB>c.HGVS`, RefSeq and Ensembl,
10/50/100/500 rows. Good enough to validate the harness before scaling.

**Scale tiers** (per issue #5): 100 → 1k → 10k → full ClinVar. Keep the random seed fixed and
commit the sampled accession list so the numbers are reproducible. **Do not** silently sample —
log N and the seed in the output.

**Sampling recipe (full):**
```
zgrep "^#" clinvar.vcf.gz > header.txt
zgrep -v "^#" clinvar.vcf.gz | shuf --random-source=<(yes 42) -n 10000 > sample.vcf
cat header.txt sample.vcf | bgzip > clinvar_10k.vcf.gz
# then VEP --pick --hgvs --refseq  →  extract (g,c) pairs
```

## 2. Systems under test

cdot is a *data provider* for biocommons/hgvs; the engine (parse, c_to_g, normalise) is held
constant and we swap the provider. That isolates the contribution of the data layer.

| System | How | Needs |
|--------|-----|-------|
| **cdot local JSON.gz** | `JSONDataProvider([...json.gz])` | download release json.gz (§4) |
| **cdot REST** | `RESTDataProvider()` → `cdotlib.org` | network only (tests new domain too) |
| **UTA remote** | `hgvs.dataproviders.uta.connect()` | network; PostgreSQL at uta.biocommons.org |
| **UTA local** | `uta.connect(db_url=...)` | local PostgreSQL + UTA dump (docker) |
| **Ensembl TARK** | `EnsemblTarkDataProvider()` | network; Ensembl-only, online-only |
| **PyHGVS + cdot** | `JSONPyHGVSTranscriptFactory` | the legacy path; report but don't over-invest |

**Sequence fetching** (orthogonal to the transcript provider — needed for `replace_reference`
and ref-base validation): `SeqRepo` (local, fast), `FastaSeqFetcher` (genome FASTA, handles
gaps), or `EnsemblTarkSeqFetcher`. For a pure *coordinate-projection* benchmark we can run
`AssemblyMapper(replace_reference=False)` and skip sequence fetching entirely — the c→g position
math doesn't need it. For a *correctness* benchmark that also checks the reference base, install
SeqRepo locally. Report both modes; note which one each number used.

## 3. Metrics & result tables

For every (system, dataset) cell, classify each input into exactly one bucket:

- **correct** — c→g equals ground-truth g.HGVS (after normalisation)
- **incorrect** — converted but disagrees (alignment/normalisation difference — investigate, see #95)
- **no_data** — transcript/version not in provider (`HGVSDataNotAvailableError`)
- **error** — parse/other failure

**Table 1 — Resolution (headline).** % correct + % no_data, per system, per build, per source
(RefSeq / Ensembl). The cdot-vs-UTA delta in `no_data` is the coverage story (1.3M vs ~141k).

**Table 2 — Speed.** load time (s), throughput (HGVS/s and tx/s), cold vs warm cache, per
system. Local JSON vs REST vs UTA. Target framing: local JSON ≫ REST > UTA remote (~1/s).

**Table 3 — Recovery.** Two sub-experiments:
1. *Cleaning* — take valid c.HGVS, inject each `clean.py` fix category (whitespace, missing
   colon, lowercase transcript, swapped gene/transcript, missing `c.`, case of dup/del, protein
   suffix, unbalanced brackets), measure: does `clean_hgvs()` restore the string so it resolves?
   Report recovery % per category. Cross-check against real failures from the search-log
   analysis (`analysis/HGVS cleaning.ipynb`) for the *real-world* distribution.
2. *Version bumping* — for each input, drop the exact requested version from the provider's
   available set, then `get_best_transcript_version()` (UP_THEN_DOWN / CLOSEST / LATEST); does
   an adjacent version resolve, and does it still convert *correctly*? This is the #28 / #27
   feature and a unique cdot capability (only meaningful given cdot's version depth). Report %
   of "version not found" that an adjacent version rescues, **and** the false-rescue rate (cases
   where the bumped version converts to a *wrong* coordinate — important caveat).

**Table 4 — Gap support (targeted).** BRCA2-class transcripts with `cDNA_match` gaps; reproduce
the Münz 2015 32%→~100% comparison with gap encoding on vs off.

**Table 5 — T2T uniqueness.** Count transcripts/HGVS resolvable only on CHM13v2.0; show one
example unresolvable by any other tool (verify against VariantValidator / VEP / TARK).

## 4. Getting cdot release data (local JSON benchmark)

Programmatic (respects schema compatibility): `cdot.data_release.get_latest_combo_file_urls(
["RefSeq","Ensembl"], ["GRCh38"])` returns the `browser_download_url`s for the matching
json.gz. Or list assets from the latest `data_v*` GitHub release. Download GRCh38 RefSeq +
Ensembl at minimum; add GRCh37 for the historical-version story. Files are large (tens–hundreds
of MB) — download once, cache, **never regenerate** from GTF/GFF for a benchmark.

## 5. Connecting to other systems

Beyond UTA, these strengthen the comparison / discussion:

- **UTA local (PostgreSQL):** the *fair* speed baseline — strips network latency so the
  comparison is data-structure vs data-structure. Run `uta` docker image, load `uta_20210129`.
- **SeqRepo (local):** for reference-base validation; biocommons stack. Install + `seqrepo pull`.
- **Ensembl TARK REST:** Ensembl-only online archive — the closest comparator for the REST mode;
  shows cdot adds RefSeq + offline + gaps.
- **VariantValidator (REST):** end-to-end validator; useful as an upper bound on "what a heavy
  online service resolves", and to confirm the T2T-novelty claim. Not a like-for-like (it
  validates; cdot only converts) — frame carefully.
- **Mutalyzer:** syntax checker — relevant to the *cleaning* axis (what fraction of our cleaned
  strings Mutalyzer then accepts).
- **ClinGen Allele Registry:** cross-reference for ground-truth g.HGVS on the Ensembl side
  (already used in the VariantGrid extraction recipe).

## 6. Harness

Extend the existing `tests/benchmark_hgvs.py` (already does c→g vs g.HGVS string compare and
tracks correct/incorrect/no_data/errors/timing) rather than starting fresh. Additions needed:
- Pluggable provider already supported (`--uta/--rest/--json/--ensembl-tark`); add `--rest`
  default now points at cdotlib.org. ✔
- Add **normalised** comparison (not just string equality) to avoid counting normalisation-only
  diffs as incorrect — use `hgvs.normalizer` on both sides before compare.
- Emit a **machine-readable** results row (CSV/JSON) per (system, dataset), not just prints, so
  the paper tables regenerate from `output/facts/`.
- Add the **cleaning** and **version-bump** sub-benchmarks (new `analysis/benchmark_resolution.py`).
- `analysis/compute_benchmark.py` already does throughput/load-time → fold its CSV output in.

Keep everything scripted + seeded so a reviewer can reproduce. Log every silent cap (sample
size, no-retry, skipped buckets).

## 6a. Results so far (first run — 2026-06-12)

Harness: `analysis/benchmark_resolution.py`, replace_reference=False, GRCh38, over the committed
`tests/test_data/clinvar_hgvs/` RefSeq pairs (100). Data: release 0.2.32 RefSeq GRCh38 (482,519
transcripts). Sequence fetcher: remote NCBI vs local `FastaSeqFetcher`
(`/data/annotation/fasta/GCF_000001405.39_GRCh38.p13_genomic.fna.gz`).

| Provider | seqfetcher | N | correct | no_data | load | throughput | total wall |
|----------|-----------|---|---------|---------|------|------------|-----------|
| cdot local JSON | NCBI (remote) | 100 | **100%** | 0 | 21 s | 0.4 HGVS/s | ~282 s |
| cdot local JSON | **local FASTA** | 100 | **100%** | 0 | 22 s | **46.0 HGVS/s** | ~25 s |
| cdot REST (cdotlib.org) | **local FASTA** | 100 | **100%** | 0 | 0.01 s | **30.9 HGVS/s** | ~3.2 s |

Findings:
- **Resolution: 100% correct, 0 no_data** — REST and local. Clean baseline; this set is "easy"
  (current MANE-ish RefSeq). Need a larger/harder set (old versions, indels, Ensembl) to surface
  the interesting failures.
- **The new cdotlib.org domain works end-to-end** through the biocommons stack (live #100 test).
- **Sequence fetching was the whole bottleneck.** Swapping remote NCBI → local FASTA took local
  JSON from 0.4 → **46 HGVS/s (≈120×)**, same correctness. The default bioutils `SeqFetcher` hits
  NCBI E-utilities at ~0.9 s/fetch, several per `c_to_g`. **Always pass a local seqfetcher
  (FASTA/SeqRepo) before quoting end-to-end throughput.** cdot's raw provider lookup
  (`get_tx_info`) is ~0 ms regardless — report provider-lookup vs end-to-end as separate numbers.
- **REST vs local JSON, once sequence is local:** REST is only ~33 % slower per-variant (30.9 vs
  46 HGVS/s) **but has no ~22 s load cost**, so for small/medium batches REST is actually faster
  end-to-end. Crossover ≈ **~2,000 variants** (load 22 s ÷ per-variant delta ~0.011 s); below
  that, REST wins; above, local JSON wins. Good framing for the paper: REST = convenience +
  small batches, JSON.gz = bulk.

### REST batching / read-ahead (answering "do we need an interface change?")

**No biocommons interface change is needed.** `RESTDataProvider` already caches per transcript in
`self.transcripts`, and a single `c_to_g` for transcript X makes **one** `/transcript/X` call
(all of `get_tx_info`/`get_tx_exons`/`get_tx_mapping_options` reuse the cache). So REST cost ≈ one
round-trip per *unique* transcript, issued **sequentially**. Over 100 variants (~88 unique tx)
that added only ~1 s here — modest, but it dominates at 10k-variant scale.

Options, cheapest-first (all cdot-specific; biocommons `Interface` untouched):
1. **Concurrent read-ahead prefetch** — add `RESTDataProvider.prefetch(tx_acs)` that warms
   `self.transcripts` using a thread pool / async before the AssemblyMapper runs. Turns N
   sequential round-trips into a few concurrent waves. **Zero server change, biggest bang/buck.**
2. **Batch REST endpoint** in cdot_rest (`POST /transcripts {ids:[...]}` or `?ids=A,B,C`) so the
   prefetch is a *single* round-trip. Small server addition; pairs with (1).
3. Leave the interface alone for the engine itself — the prefetch only warms the existing cache,
   so per-variant calls become cache hits. A true streaming/batch *Interface* change would have
   to go upstream into biocommons and isn't worth it.

Recommendation: implement (1) now (it's a few lines on `RESTDataProvider` + a `prefetch` call in
the benchmark), add (2) in cdot_rest when convenient. For genuine bulk, still steer users to
JSON.gz — REST is the convenience tier.

**Prototype result (option 1, implemented).** Added `RESTDataProvider.prefetch(tx_acs,
max_workers=10)` — concurrent (ThreadPoolExecutor) read-ahead that warms `self.transcripts`,
skips already-cached, stores `None` on failure (matches `_get_transcript`). Benchmark gains a
`--prefetch` flag. Same 100 RefSeq pairs, REST + local FASTA, GRCh38:

| REST mode | prefetch | resolution wall | resolution throughput | total wall |
|-----------|----------|-----------------|-----------------------|-----------|
| naive (sequential per-variant) | — | 3.24 s | 30.9 HGVS/s | ~3.24 s |
| **read-ahead** | 88 tx in **0.47 s** (10 workers) | **0.22 s** | **451.7 HGVS/s** | **~0.69 s** |

The 88 sequential transcript round-trips collapse into one 0.47 s concurrent wave; the resolution
loop is then a pure cache-hit pass (network removed). End-to-end **~4.7× faster** (3.24 → 0.69 s),
and the win grows with batch size and network latency. No biocommons interface change; full test
suite still green (104 passed, 5 skipped). Next: a batch REST endpoint (option 2) would shrink
the 0.47 s prefetch to a single round-trip; and `prefetch` should be exposed as public API + doc'd.

**Important correction — the flat 2-pass only holds for clean, exact-version HGVS.** cdot's job
includes fixing and version-bumping, which the simple "extract accession → prefetch → resolve"
model ignores:
- **Cleaning** folds into pass 1 for free — `clean_hgvs()` is pure-string and *produces* the
  reliable accession from messy input; run it as pass 1 (you run it anyway).
- **Version bumping** breaks the flat model: the string says `.4`, but if `.4` is absent the
  bumped target (`.3`/`.5`) is **not knowable from the string** — it depends on which versions the
  provider holds. **Gene-only HGVS** is the same (canonical depends on provider tags).
- So the real shape is a **bounded warm / resolve-misses loop**: clean → warm requested versions +
  gene lists → resolve → for misses, discover (versions/canonical) → warm bumped targets → resolve
  misses. This iterates over *cache-warming waves*, NOT over suspended mid-parse HGVS state (the
  mapper still runs per-variant), so a `c_to_g` coroutine is still unnecessary.
- This makes the batch/discovery endpoint **load-bearing**: version discovery needs a
  "list versions" / "best version" capability, else it's N probe round-trips and the pre-warm gain
  evaporates on the messy-HGVS path. Tracked: cdot#108, cdot#109, cdot_rest#9.

## 7. Risks / caveats

- **Never run the full pipeline** on a whim — full ClinVar × all providers × UTA-remote (~1/s)
  is hours. Tier up: 100 → 1k → 10k. Quote wall-clock before going bigger.
- **"incorrect" is not always cdot's fault** — UTA-vs-annotation alignment differences (#95)
  show up as disagreements. Manually audit a sample of incorrects before reporting; some are a
  paper *paragraph*, not a bug.
- **Normalisation** differences (3'-shift, ref-base) inflate "incorrect" if you string-compare —
  normalise first.
- **Sequence dependence dominates end-to-end speed** (see §6a) — even with `replace_reference=False`,
  biocommons `hgvs` normalisation fetches sequence; the default remote NCBI `SeqFetcher` is
  ~1 s/fetch and swamps everything. Install local SeqRepo before quoting any end-to-end
  throughput, and keep provider-lookup throughput as a separate number.

## 8. Concrete next steps

1. [x] Validate harness on `tests/test_data/clinvar_hgvs/*.tsv` over cdot REST (see
   `analysis/benchmark_resolution.py`, first run logged below).
2. [ ] Download GRCh38 RefSeq + Ensembl json.gz; rerun local vs REST vs UTA for speed (Table 2).
3. [ ] Install SeqRepo; rerun resolution with `replace_reference=True` + normalised compare (Table 1).
4. [ ] VEP-annotate a 1k / 10k ClinVar sample for a fresh, larger resolution set.
5. [ ] Cleaning sub-benchmark: inject fix categories, measure recovery (Table 3.1).
6. [ ] Version-bump sub-benchmark: drop version, measure rescue + false-rescue (Table 3.2).
7. [ ] Stand up UTA-local for the fair speed baseline.
8. [ ] BRCA2 gap experiment (Table 4) and T2T-uniqueness count (Table 5).
