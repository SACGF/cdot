# Paper rewrite: outline

Working file for rewriting the main text in your own words. This outline is the *target*
structure: paragraphs have already been split, merged and moved where the current text was
badly built (see "Structural changes" below). Write from the beats, not from the current
text.

Each entry gives:

- **Job**: what the paragraph does for the argument.
- **Opens with**: the topic. The reader should know what the paragraph is about by the end
  of sentence one or two.
- **Lands on**: the point, in the last sentence.
- **Beats**: the content in order, as notes rather than sentences.
- **Must keep**: `{{ fact }}` references and `[@citations]`. Wording a claims sidecar checks is
  marked *checked*. You can reword it, but the claim must stay true. The sidecar gets updated
  afterwards.
- **Your version**: blank until written. Your dictation is recorded as-is. The accepted text
  only has minimal corrections (dictation slips, factual fixes, facts in your own phrasing). A
  heavier "tidy" by Claude scored 100% AI on Pangram, so that isn't done.

Main text only. Supplementary and the figure caption are not included. Table captions are
listed but are optional to rewrite.

## Recurring tells to avoid

- **"X, not Y" / "X rather than Y" endings.** The current text has more than a dozen. Say what
  something is, and leave out what it isn't unless a reader would really assume the opposite.
- **Long colon and semicolon chains.** Prefer two or three plain sentences.
- **Signpost openers**: "cdot is, first,", "Comparing like with like:", "Caution is
  warranted", "The roles are converging, though.", "the real test of".
- **Saying it twice.** Each point now has one home in the outline. If you want to refer back,
  a short clause is enough.

## Structural changes from the current text

Listed so you can veto any of them.

1. **Intro.** I1 is split into "needs versioned data" (I1a) and "strings are broken" (I1b).
   I2 is split into the libraries (I2a) and UTA's limits (I2b). I3 now ends on the
   three-part approach. The coverage figure and formats move to a short new I4, which says
   what the paper provides and tests.
2. **Methods, cleaning.** The motivation paragraph (copy/paste damage list) is merged into
   the `clean_hgvs()` paragraph as one opening sentence. The damage list itself stays in
   Results R4.1, where the data shows it.
3. **Methods, clients.** The biocommons and PyHGVS paragraphs are merged into one.
4. **Methods, benchmarking.** Reordered to follow the Results: private-data note first, then
   resolution corpora, throughput, cleaning evaluation, version stability. The VV/Mutalyzer
   method joins the LOVD method. Version stability, previously tacked onto the VV paragraph,
   gets its own.
5. **R1.** The VEP gap-handling point is split out of the coverage paragraph into R1b.
6. **R2.** The one-sentence lead-in merges into R2.1. The method detail in the historical
   corpus paragraph is cut to a reference to Methods. "Cleaning rescued N on the submitted
   corpus" moves to R4, where the two corpora are contrasted. In the Shariant paragraph, the
   "no ground truth" caveat moves up so the paragraph ends on the breakdown.
7. **R3.** It now opens on "faster than UTA locally and remotely". Load time moves to the
   full-scale paragraph.
8. **R4.** The production result is split into headline (R4.1) and fix breakdown (R4.2). The
   ceiling (R4.3) follows directly. The two-corpora contrast becomes its own paragraph (R4.4).
   **The LOVD and VV/Mutalyzer numbers move from the Discussion into Results (R4.6)**. The
   Discussion keeps only what they mean.
9. **R5.** The one-sentence lead-in is merged. The very long drift paragraph is split into
   three: how often it's safe (R5.1), what kind of drift (R5.2), and where the drift is
   concentrated (R5.3).
10. **Discussion.** Reordered. It opens on what the results showed (new), then adoption as
    evidence of general use, related repair tools, how they compare, and limitations, and it
    ends on the design principle (safe cleaning automatic, heuristics opt-in) rather than on a
    caveat. The Ensembl REST and TARK asides are dropped, since TARK is already in Methods.
11. **Abstract.** Add a cleaning result.

---

# Introduction

**Flow:** what resolution needs (I1a) → the strings are broken (I1b) → the libraries and
their usual backend (I2a) → why that backend falls short (I2b) → cdot's origin and its
three-part answer (I3) → what the paper provides and tests (I4).

### I1a. Resolving HGVS needs versioned transcript data

- **Job:** Define the task the paper is about.
- **Opens with:** HGVS is the standard way to describe variants.
- **Lands on:** resolving it needs the right transcript version and its alignment to the build.
  This is the requirement that I2b says UTA fails.
- **Beats:**
  1. HGVS is the international standard (reports, databases, papers).
  2. Scale: ClinVar alone has more than N variants with HGVS.
  3. c./n. to genome needs versioned transcript data: the accession names a sequence at a
     version (`NM_000492.3`), and its alignment to a build places the exons.
- **Must keep:** `literature.clinvar_variants`; `[@DenDunnen2016]`, `[@Landrum2025]`.
- **Your version (as dictated):** HTVS is the clinical standard for describing variants.
  People often like using C and P.dot forms because they find it easier to reason about. But
  in genomics, we have to convert it to genomic coordinates. So in order to do this, you need
  to have the reference transcript to have a version. And these versions of the transcript
  are aligned to the genome by RefSeq and Ensemble, and they provide these. And then you can
  take those alignments and do the calculations to convert them to genomic coordinates.
- **Accepted (minimal corrections; Pangram 100% human):** HGVS is the clinical standard for
  describing variants [@DenDunnen2016]. People often like using c. forms because they find them
  easier to reason about. But in genomics, we have to convert them to genomic coordinates. So
  in order to do this, you need to have the reference transcript at a specific version. These
  versions of the transcript are aligned to the genome by RefSeq and Ensembl. You can take
  those alignments and do the calculations to convert them to genomic coordinates.

### I1b. Human-entered HGVS is often broken

- **Job:** The second problem, the strings themselves. This sets up cleaning.
- **Opens with:** HGVS written by people is often wrong.
- **Lands on:** existing tools fix only a few of those strings.
- **Beats:**
  1. Mutalyzer logs, about 26M unique descriptions over 5 years [@Lefter2021].
  2. Only ~X% correct as written, ~Y% had a syntax or semantic error.
  3. Mutalyzer could auto-correct only ~Z%.
  4. Optional, from you: where broken strings come from in practice (pasted reports, lab
     systems), or what that meant for Shariant.
- **Must keep:** `literature.mutalyzer_correct_pct`, `literature.mutalyzer_error_pct`,
  `literature.mutalyzer_autocorrect_pct`; `[@Lefter2021]`. Quoted literature (*checked*).
- **Your version (as dictated):** One of the surprising things we've seen is that HGVS
  strings are often broken in a way that suggests that humans have either just typed them out
  by hand or copy-pasted them incorrectly. Looking through the literature, we find other tools
  report the x% errors. and could fix only another certain %.
- **Accepted (minimal corrections; Pangram 100% human):** One of the surprising things we've
  seen is that HGVS strings are often broken in a way that suggests that humans have either
  just typed them out by hand or copy-pasted them incorrectly. Mutalyzer looked at 26 million
  descriptions submitted over five years and found only half were correct,
  {{ literature.mutalyzer_error_pct | pct(0) }} had errors, and it could only fix
  {{ literature.mutalyzer_autocorrect_pct | pct(0) }} [@Lefter2021].

### I2a. The Python HGVS libraries and their data backends

- **Job:** Introduce the libraries and the idea of a swappable data backend.
- **Opens with:** two dominant Python HGVS libraries.
- **Lands on:** the standard backend for biocommons/hgvs is UTA, a PostgreSQL database.
- **Beats:**
  1. biocommons/hgvs [@Wang2018] and PyHGVS.
  2. biocommons/hgvs reads transcript data through a pluggable provider interface. Other tools
     (hgvs-weaver [@HgvsWeaver]) implement it too, so one backend serves several clients.
     PyHGVS uses its own structures.
  3. The standard backend is UTA [@Hart2020], which stores alignments in PostgreSQL.
- **Must keep:** `[@Wang2018]`, `[@HgvsWeaver]`, `[@Hart2020]`.
- **Your version (as dictated):** The two main Python HGVS libraries are Biocommons HGVS and
  PyHGVS. Biocommons gets its transcript genome data from the data provider interface, which
  has become a bit of a standard, having been adopted by a number of other tools. PyHVS
  provides its own simple data structures.

  The default and only implementation of the Biocommons HTTPS data provider library is UTA,
  the Universal Transcript Archive. And this works by performing its own transcript to genome
  alignments and storing it in a Postgres database.

### I2b. UTA's limits: access and coverage

- **Job:** The gap cdot fills.
- **Opens with:** UTA has two problems, access and coverage.
- **Lands on:** no Ensembl at all.
- **Beats:**
  1. Access: install PostgreSQL locally, or connect remotely over the PostgreSQL protocol.
     Hospital and corporate networks often block that, and it's slow (~1 transcript/s at best,
     N HGVS/s in our measurement, Table 1).
  2. Coverage: ~N versioned alignments in the 2021 release, RefSeq only, limited history.
  3. Ensembl is absent entirely, though it's widely used in research and the main source in
     Europe.
- **Must keep:** `benchmark.uta_remote_tps`, `literature.uta_count`; `[@Hart2020]`.
  *Checked:* "about one transcript per second at best"; "RefSeq only ... absent entirely".
- **Your version (as dictated):** This is a great resource, but has a number of drawbacks.
  Firstly, doing its own alignment means that it's quite an effort to do this, and for
  instance, Ensemble is not supported. The other drawback of a Postgres SQL database is that
  either if you set it up locally, it requires quite a bit of configuration work, or you
  connect over the internet, many workplaces or clinical hospitals basically only allow HTTPS
  out and block postgres conneections
- **Accepted (your words, cut and joined only; Pangram 100% human):**

  The two main Python HGVS libraries are biocommons HGVS [@Wang2018] and PyHGVS. Biocommons
  gets its transcript genome data from the data provider interface, which has become a bit of
  a standard, having been adopted by a number of other tools, including VariantValidator
  [@Freeman2018] and MaveDB [@Rubin2025]. We found
  {{ adoption.n_interface_implementers | number_word }} by searching public GitHub
  repositories, so this is a lower bound. PyHGVS provides its own simple data structures.

  The default implementation of the biocommons HGVS data provider is UTA, the Universal
  Transcript Archive [@Hart2020]. This works by performing its own transcript to genome
  alignments and storing them in a Postgres database.

  This is a great resource, but has a number of drawbacks. Firstly, doing its own alignment
  means that it's a large amount of effort to add new transcripts, and for instance, Ensembl
  is not supported. The total number of transcript alignments is
  {{ coverage.uta_count | commas }}, RefSeq only, and it provides a limited amount of
  historical support. The other drawback of a PostgreSQL database is that if you set it up
  locally, it requires quite a bit of configuration work. Many workplaces, clinical labs and
  hospitals firewall off all connections other than HTTPS, which means that they are unable
  to connect to the public Postgres instance.
  The protocol for communication with UTA requires several round trips to perform the
  conversion, which means the performance is {{ benchmark.uta_remote_tps | dp(2) }} HGVS/s
  (Table 1).

### I3. Shariant, and cdot's three-part approach

- **Job:** Why cdot exists and how it attacks the problem.
- **Opens with:** cdot was built for Shariant.
- **Lands on:** the three parts: coverage, cleaning, safe substitution.
- **Beats:**
  1. Shariant [@Tudini2022]: a national platform pooling lab classifications.
  2. Classifications build up over years, each against whatever transcript version the lab
     used.
  3. To resolve as many as possible, cdot does three things: broad versioned coverage (the
     cited version is there), string cleaning (a malformed string still parses), and opt-in,
     coordinate-safe version substitution (for when the version really is gone).
- **Must keep:** `[@Tudini2022]`.
- **Your version (as dictated):** CDOT was split out of the Australian Genomic Shariant
  project. The goal of Shariant is to collect HGV, or to collect variant classifications from
  the different labs across Australia. We found the most common form that these variants were
  stored in was HGVS, but they came from dozens of different systems. So there were multiple
  commercial systems, custom built software solutions, and also spreadsheets. And these often
  went back many, many years. So one of the great problems in Shariant was how to resolve all
  of these HGVSs including transcript versions that are many years old.

  We also found that one of the most common features in the web interface was people entering
  HGVSs in the search box, looking for existing classifications against this. And we found
  that people were using all kinds of invalid HGVSs, and we basically tried to do all we can to
  resolve these to what the people want, the people actually intended.

  there are a number of problems that we had to solve. So there was the cleaning. There was
  also the just collecting as many historical transcript versions as we could. And the other
  problem was a lot of, well, some labs were on GRCH 38 and some on 37, and we had to convert
  between them. And there were often no official releases of a transcript version against a
  particular genome build. So what we did then was try and find the nearest version,
  transcript version available for that build that still had the same exonic structure and
  thus was safe to substitute.
- **Accepted (your words, cut and joined only; Pangram 100% human):**

  cdot was split out of the Australian Genomics Shariant project [@Tudini2022]. The goal of
  Shariant is to collect variant classifications from the different labs across Australia. We
  found the most common form that these variants were stored in was HGVS, but they came from
  dozens of different systems. There were multiple commercial systems, custom built software
  solutions, and also spreadsheets. These often went back many, many years. So one of the
  great problems in Shariant was how to resolve all of these HGVSs including transcript
  versions that are many years old.

  We also found people entering HGVSs in the search box, looking for existing
  classifications. People were using all kinds of invalid HGVSs, and we tried to resolve these
  to what the people actually intended.

  We had to solve a number of problems. There was the cleaning. There was also
  collecting as many historical transcript versions as we could. And the other problem was
  some labs were on GRCh38 and some on GRCh37, and we had to convert between them. There were
  often no official releases of a transcript version against a particular genome build. So
  what we did then was try and find the nearest transcript version available for that build
  that still had the same exonic structure and thus was safe to substitute (opt-in and
  flagged).


### I4. What cdot provides and what the paper tests

- **Job:** The concrete product and a preview of the evaluation.
- **Opens with:** cdot builds JSON from every RefSeq and Ensembl release.
- **Lands on:** what the paper evaluates (on ClinVar against UTA: coverage, accuracy on
  historical strings, speed, cleaning, substitution safety).
- **Beats:**
  1. Compact JSON from all RefSeq GFF3 and Ensembl GTF releases: N alignments across GRCh37,
     GRCh38, T2T-CHM13v2.0 [@Nurk2022].
  2. Loads locally in seconds, also served by REST (cdotlib.org).
  3. New: one sentence on what the paper measures.
- **Must keep:** `coverage.total_count`; `[@Nurk2022]`. *Checked:* "from all available RefSeq
  GFF3".
- **Your version (as dictated):** So what CDOT does is it goes and reads the official GFF or
  GTF releases from RefSeq and Ensemble, and then converts them into a JSON format that allows
  us to basically collect all of the 1.9 million version transcript alignments across 37, 38,
  and T2T genome builds. We conver that into JSON a widely supported data format. Almost all
  languages have libraries that have been highly optimized so that you can load it very
  rapidly. The same data can also be transferred across the internet via a REST server, which
  means that workplaces with web access can retrieve them without local copies.

  We provide our cleaning library, which goes and fixes a lot of common mistakes, giving a
  record of what substitutions were made.

  We validate and benchmark this by comparing C. and versus UTA on the historical ClinVar
  submissions, the VCV XML, which records the original transcript version that was submitted
  at the time, giving a historical test record. We also benchmark, we also use our private
  search history of invalid HGVSs and compare our fix rate versus other HGVS fix tools,
  including X, whatever they are.
- **Accepted (your words, cut and joined only):**

  cdot reads the official GFF or GTF releases from RefSeq and Ensembl, and then converts them
  into a widely supported JSON format that allows us to collect all of the
  {{ coverage.total_count | commas }} versioned transcript alignments across GRCh37, GRCh38,
  and T2T genome builds [@Nurk2022]. Almost all languages have libraries that have been highly
  optimized so that you can load it very rapidly. The same data can also be transferred across
  the internet via a REST server, which means that workplaces with web access can retrieve them
  without local copies. We provide our cleaning library, which fixes a lot of common mistakes,
  giving a record of what fixes were made.

  We validate and benchmark this by comparing cdot versus UTA on the historical ClinVar
  submissions, the VCV XML, which records the original transcript version that was submitted
  at the time, giving a historical test record. We measure performance in HGVS resolved per
  second, as well as the number resolved. To test the HGVS cleaning and fixes, we take public
  ClinVar data and then synthetically modify it to match the problems we found in the invalid
  HGVSs from our own private search data, and we compare our fixes to other tools, including
  LOVD, VariantValidator and Mutalyzer.

---

# Methods

**Status: the move below has been applied to the files.** `methods.md` holds cut-down old
text as placeholders for these six paragraphs, and the full old Methods is in
`supplementary.md` ("... (detail)" subsections).

**First pass: most of Methods moves to Supplementary Methods.** The current Methods is about
1,500 words. The plan below keeps six short paragraphs in the main text (about 500 words) and
moves the detail to the supplement, saving about 1,000 words. That roughly cancels the
overrun. The moved text stays as it is for now. It's supplementary, so it can be lightly
edited later rather than rewritten.

**Flow:** where the data comes from and how it is built (M1) → the file format (M2) →
cleaning (M3) → version substitution (M4) → how clients reach the data (M5) → how the
results were measured (M6).

## Moving to Supplementary Methods

| Current text | Goes to |
|---|---|
| Per-source release counts, RS_2025_08 / RS_2024_08 detail, why GTF not GFF3 | Tables S1/S2 notes (new "Data sources" subsection) |
| Pipeline steps (HTSeq, bioutils, CDS from codons, gene metadata) | new "Data generation pipeline" subsection |
| JSON per-build fields, gap-to-CIGAR conversion, unaligned bases between exons, canonical tags | Table S3 (already the schema) |
| `HGVSFix` fields, unfixable cases, the accession-prefix repair in full | Table S7 (already the operation catalogue) |
| Substitution strategies, the three safety conditions in full | "Version-substitution safety gate detail" (exists) |
| Per-provider detail (lazy interval trees, `prefetch()`, TARK), code example, `FastaSeqFetcher` RefSeq caveat, `ChainedSeqFetcher`, PyHGVS | new "Data providers and clients" subsection (the RefSeq caveat subsection exists) |
| `resolve_gene_hgvs()` / MANE selection | same new "Data providers and clients" subsection |
| Script names, ClinVar pair construction, submitted-corpus build | "Submitted-string corpus construction" (exists) |
| LOVD and VariantValidator/Mutalyzer scoring rules, service dates | "LOVD head-to-head scoring detail" (exists, widened) |
| Throughput protocol (N passes, IQR, shared SeqRepo, remote UTA subset) | "Throughput measurement detail" (exists) |
| Version-stability sample and positional bins | "Version-substitution safety gate detail" (exists) |

Facts and claims checks travel with their text (`sources.*`, `benchmark.n_*`,
`lovd_comparison.lovd_version`, `vv_mutalyzer_comparison.service_date`,
`version_stability.sample_n`, and the `methods.claims.toml` quotes).

### M1. Data sources and generation

- **Job:** Where the data comes from and how it becomes one file.
- **Opens with:** cdot merges every historical RefSeq GFF3 and Ensembl GTF release, plus UTA,
  per genome build (Figure 1A).
- **Lands on:** newest wins, which also backfills fields missing from early files.
- **Beats:**
  1. RefSeq GFF3 [@OLeary2016] and Ensembl GTF [@Martin2023], every archived release, for
     GRCh37, GRCh38 and T2T (Tables S1, S2).
  2. UTA alignments fill gaps from before the oldest releases, and official annotation
     overrides them.
  3. Merged oldest to newest, newer wins. This also backfills fields that early files lack.
  4. A Snakemake pipeline does the parsing and writes gzipped JSON (Supplementary Methods).
- **Must keep:** `[@OLeary2016]`, `[@Martin2023]`. (The `sources.*` counts move to the
  supplement.)
- **Your version (as dictated):** CDOT data collection works via a snake make pipeline that
  downloading all of the official RefSeq GFF3 and Ensemble GTF releases for 37, 38, and T2T.
  We also use download UTA, which, as it computes its own alignments, has a few accessions
  that were never in any of the official releases.

  We convert these the the JSON format

  We then also provide a release which takes all of the historical releases for such as
  RefSeq for 38 and merges them chronologically, with older data being replaced by new, which
  sometimes contains new information such as cigar strings for alignment gaps.

  *Corrections (as dictated):* The historical versions are loaded with UTA first and then
  chronologically. So the official releases will override any custom UTA, and older versions
  will be overwritten by newer ones which contain new information such as gaps. Say any
  custom alignments performed by UTA.

  *Corrections after acceptance (as dictated):* We also provide a cdot release which includes all of the historical releases RefSeq or Ensemble have made for a particular genome build.
- **Accepted (your words, cut and joined only):**

  cdot data collection works via a Snakemake pipeline downloading all of the official RefSeq
  GFF3 [@OLeary2016] and Ensembl GTF [@Martin2023] releases for GRCh37, GRCh38, and T2T
  (Figure 1A; Supplementary Tables S1 and S2). We also download UTA, which, as it computes its
  own alignments, has a few accessions that were never in any of the official releases. We
  also provide a cdot release which includes all of the historical releases RefSeq or Ensembl
  have made for a particular genome build. The historical versions are loaded with UTA first
  and then chronologically, so the official releases will override any custom alignments
  performed by UTA, and older versions will be overwritten by newer ones which contain new
  information such as gaps.


### M2. JSON format

- **Job:** What a file holds.
- **Opens with:** each transcript entry holds its metadata plus its alignment on each build.
- **Lands on:** `schema_version`, so clients can reject incompatible files (full schema in
  Table S3).
- **Beats:**
  1. One entry per transcript version, with a `genome_builds` dict, so one entry carries all
     builds. A single-build file can stand in for a GTF/GFF.
  2. Per-exon alignment gaps (transcript vs genome indels), applied during conversion.
  3. Canonical tags (MANE, RefSeq Select, Ensembl canonical).
  4. `schema_version`.
- **Must keep:** nothing fact-bound. (I4 already says why JSON, so don't repeat it.)
- **Your version (as dictated):** The JSON format has one entry per transcript version with
  build specific alignment data in a genome_builds dictionary kedyed by build. This gives the
  same structure for a JSON result, whether there is information to align it to 1, or multiple
  builds.

  The Exxon data contains a gap field that records the indels between found during alignment,
  and can be applied for covering gaps during conversion. We can tag transcripts marking them
  as MANE select and other information

  We include a top-level schema version, which can be used to check what information it
  contains and whether it can be parsed by older clients.

  *Corrections (as dictated):* Covering was used to describe handling. And for the, I guess,
  we copy tags from the original GFF, such as main plus clinical, ref select, and sonical
  canonical.

  For the JSON, I was describing the layout of the JSON data is the same whether genome builds
  contain one or multiple genome builds, which means loadder we can handle it the same from
  various sources

  *Corrections (as dictated):* add in main select and GFF or GTF.

  *Corrections after acceptance (as dictated):* For M2, I think we can just say that it's consistent and then cut everything else.
- **Accepted (your words, cut and joined only):**

  The JSON format has one entry per transcript version with build-specific alignment data in a
  `genome_builds` dictionary keyed by build. The layout of the JSON data is the same whether
  `genome_builds` contains one or multiple genome builds. The exon data contains a `gap` field
  that records the indels found during alignment, and can be applied for handling gaps during
  conversion. We copy tags from the original GFF or GTF, such as MANE Select, MANE Plus
  Clinical, RefSeq Select, and Ensembl canonical. We include a top-level `schema_version`,
  which can be used to check what information it contains and whether it can be parsed by
  older clients (Supplementary Table S3).


### M3. String cleaning

- **Job:** The cleaning API: `fix_hgvs()` (the recommended entry point) and `clean_hgvs()`
  (the pure string function inside it, for callers who want only that).
- **Opens with:** `fix_hgvs()` is the one call we recommend: HGVS string in, repaired string
  plus a record of each fix out.
- **Lands on:** `clean_hgvs()` is available on its own as a pure function, for callers who want
  string repair with no data.
- **Beats:**
  1. `fix_hgvs()` always cleans the string first (a fixed order of small steps, Table S7).
  2. With a data provider it also: restores a dropped accession prefix (only when exactly one
     candidate exists in the loaded data) and, given a genome build, resolves a gene symbol to
     its MANE transcript [@Morales2022; @Wright2023].
  3. Each fix is recorded (`HGVSFix`), so a caller can audit what changed. It never raises an
     error by default, and it flags what it cannot fix.
  4. The string step is exposed as `clean_hgvs()`: needs no genome, sequence, parser or data
     provider, so it runs anywhere, before any library.
- **Must keep:** `[@Morales2022; @Wright2023]` (moved here from M5, since gene resolution is
  part of `fix_hgvs()`).
- **Your version (as dictated):** For cleaning HTTPS strings, we provide a method fixHTTPS. So
  you pass in an HTTPS string and get back the repaired string plus a list of the fixes
  applied. Some examples include whitespace quoting, double dots, the uppercase, lowercase,
  and or the swapped gene and transcript examples. If you pass a data provider to fix HTTPS, it
  can also look up, fix some transcript errors, such as a dropped accession prefix. It only
  does this when exactly one candidate is possible, and give the example you gave. If you also
  pass a genome build, it resolves a gene-only string, such as BRCA2 C68del, by looking up the
  main transcript. And so every change is recorded as an HTTPS fix enumeration, which has
  severity, code, and before after, so a caller can audit what was changed or present them to
  the user. What it can't fix is flagged as an error level fix. By default, it returns errors
  instead of raising them, but you can pass a flag to throw an exception. Internally, we also
  have another function called clean HTTPS, which is called by fix HTTPS, but it is a pure
  string conversion function and doesn't load any data

  *Corrections (as dictated):* For the examples, just sub in the examples that we used
  previously, because there's no way to tell the difference between just including examples,
  AI and humans. Yeah, change method fix HVS to function. And yeah, feel free to trim it down
  if there's some superfluous words there.
- **Accepted (your words, cut and joined only):**

  For cleaning HGVS strings, we provide a function `fix_hgvs()`. You pass in an HGVS string
  and get back the repaired string plus a list of the fixes applied. Some examples include
  whitespace, quoting, double dots (`NM_000059..4`), uppercase, lowercase (`c.123delg`), and
  the swapped gene and transcript (`BRCA2(NM_000059.4)` → `NM_000059.4(BRCA2)`)
  (Supplementary Table S7). If you pass a data provider to `fix_hgvs()`, it can also fix some
  transcript errors, such as a dropped accession prefix (`000059.4:c.68del` →
  `NM_000059.4:c.68del`). It only does this when exactly one candidate is possible. If you
  also pass a genome build, it resolves a gene-only string, such as `BRCA2:c.68_69del`, by
  looking up the MANE transcript [@Morales2022; @Wright2023]. Every change is recorded as an
  `HGVSFix`, which has severity, code, and before/after, so a caller can audit what was
  changed or present them to the user. What it can't fix is flagged as an error-level fix. We
  also have `clean_hgvs()`, which is called by `fix_hgvs()`, but it is a pure string
  conversion function and doesn't load any data.

### M4. Version substitution

- **Job:** The opt-in helper and its safety check.
- **Opens with:** when the cited version is missing, `get_best_transcript_version()` offers the
  nearest available version, opt-in and reported.
- **Lands on:** it is only offered as safe when the coding structure matches (validated in
  R5).
- **Beats:**
  1. Opt-in, returned as an `HGVSFix`. The caller decides.
  2. Safety check: a version's intrinsic CDS structure (CDS length and coding-exon lengths in
     transcript coordinates) is the same on every build, so two versions can be compared even
     where one was never aligned to that build.
  3. Safe only if the structures match (plus gap and UTR conditions, Supplementary Methods).
- **Must keep:** nothing fact-bound.
- **Your version (as dictated):** When the requested transcript version isn't available, Fix
  HTTPS can substitute in another version. This is off by default because it can be wrong. You
  can pass in a strategy, version fallback, which you can set to either try higher versions
  than lower, closest, or the latest. It's reported back as an HTTPS fix, so the user can see
  what happened. We have a check that each version has, transcript version has a CDS
  structure. The CDF length plus the length of coding exons. This is the same on every genome
  build, so the two versions can be compared even when one of them is not aligned or not
  available for the genome build you request. For instance, if you can read the old version
  structure from 37 and compare it to the current version in 38. It's safe to swap only if the
  structures match. There are two extra conditions which are gone into in the supplementary
  methods. The alignment gaps must also match, and for a variant in the UTR length must match
  as well.

  *Corrections (as dictated):* We have a check that each transcript version has an identical
  CDS structure

  For a UTR variant, the UTR lengths must be the same.

  If the version substitution doesn't match, by default we don't just leave the string
  unchanged and it returns an error HTTPS fix. You can override this to make it only a
  warning.

  You can override this to change the string and return a warning that the coordinates may
  have moves

  *Edits after acceptance (author: "Feel free to chop words when you think it makes sense"):* cut "either", the GRCh37/GRCh38 example sentence, and the "two extra conditions" sentence (Supplementary Methods reference moved onto the next sentence).
- **Accepted (your words, cut and joined only):**

  When the requested transcript version isn't available, `fix_hgvs()` can substitute in
  another version. This is off by default because it can be wrong. You can pass in a strategy,
  `version_fallback`, which you can set to try higher versions then lower, closest, or the
  latest. It's reported back as an `HGVSFix`, so the user can see what happened. We have a
  check that each transcript version has an identical CDS structure: the CDS length plus the
  length of coding exons. This is the same on every genome build, so the two versions can be
  compared even when one of them is not aligned or not available for the genome build you
  request. It's safe to swap only if the structures match. The alignment gaps must also match,
  and for a UTR variant, the UTR lengths must be the same (Supplementary Methods). If the
  version substitution doesn't match, by default we just leave the string unchanged and it
  returns an error `HGVSFix`. You can override this to change the string and return a warning
  that the coordinates may have moved.


### M5. Data providers and client libraries

- **Job:** How software reaches the data.
- **Opens with:** cdot implements the biocommons/hgvs data provider interface, so it drops in
  for UTA.
- **Lands on:** sequence can come from SeqRepo or a local genome FASTA, so it can run fully
  offline.
- **Beats:**
  1. Three providers: local JSON (in memory), REST (cdotlib.org, with batch prefetch), and
     Ensembl TARK.
  2. PyHGVS supported but unmaintained.
  3. Sequence from SeqRepo [@Hart2020] or `FastaSeqFetcher`, tried in order.
- **Must keep:** `[@Hart2020]`. (The MANE lookup and its citations moved to M3.)
- **Your version (as dictated):** C.Dot implements a Biocommons HTTPS data provider, so it
  can be swapped in to load the C.Dot transcript versions. We provide 3 providers: JSON data
  provider, which loads the gzip JSON into memory; REST data provider, which makes network
  calls to c.lib.org. It also has a prefetch, which batches things into more efficient calls;
  and Ensemble Tark data provider, which reads from Ensemble's official transcript archive,
  Tark. We also make a PyHVS data loader, but PyHVS is no longer maintained. In order to load
  we also need to be able to load the sequences. So the Biocommons HTTPS provide seq repo, but
  we also make a seq fetcher that can load a local FASTA file, which pulls out the exon
  sequences from the FASTA file and pastes them together to make a transcript. This works
  perfectly for Ensemble transcripts, which exactly match the exon, the genome, but the might
  have some troubles with refseq, where the transcript reference and the genome reference can
  have gaps and mismatches. In this case, we recommend that we've also made a chained seq
  fetcher, which can try several data sources and then fall back on others. And you can use
  this to take the seq repo sequences when available and then fall back on the other sources.

  *Corrections (as dictated):* Just cut what you think we can.

  *Edits after acceptance:* cut "We also need to be able to load the sequences."
- **Accepted (your words, cut and joined only):**

  cdot implements a biocommons HGVS data provider, so it can be swapped in to load the cdot
  transcript versions (Figure 1B). We provide three providers: `JSONDataProvider`, which loads
  the gzipped JSON into memory; `RESTDataProvider`, which makes network calls to cdotlib.org
  (it also has a prefetch, which batches things into more efficient calls); and
  `EnsemblTarkDataProvider`, which reads from Ensembl's official transcript archive, TARK.
  biocommons provides SeqRepo [@Hart2020], but we also make a seq fetcher (`FastaSeqFetcher`)
  that can load a local FASTA file, which pulls out the exon sequences and pastes them
  together to make a transcript. This works perfectly for Ensembl transcripts, which exactly
  match the genome, but might have some troubles with RefSeq, where the transcript reference
  and the genome reference can have gaps and mismatches. We've also made a chained seq fetcher
  (`ChainedSeqFetcher`), which can try several data sources and then fall back on others. You
  can use this to take the SeqRepo sequences when available and then fall back on the other
  sources.


### M6. Evaluation

- **Job:** One overview of how each result was measured. The detail is in the supplement.
- **Opens with:** every measurement comes from scripts committed in `paper/scripts/`.
- **Lands on:** what **[private data]** means (frozen aggregates, everything else reproducible
  from public data).
- **Beats:**
  1. Resolution (R2): ClinVar pairs, current-version and submitted strings (VCV XML 2026-10),
     scored against ClinVar's VCF coordinate.
  2. Throughput (R3): the same ClinVar pairs through every backend, with engine and sequence
     layer held fixed, so only the transcript data varies.
  3. Cleaning (R4): the production search corpus, plus corrupted public ClinVar strings run
     through cdot, LOVD [@LovdHgvsChecker], VariantValidator [@Freeman2018] and Mutalyzer
     [@Lefter2021].
  4. Substitution safety (R5): consecutive versions compared on GRCh38.
  5. [private data] marker.
- **Must keep:** `[@LovdHgvsChecker]`, `[@Freeman2018]`, `[@Lefter2021]`.
- **Your version (as dictated):** a reproducible script to calculate all of the benchmarks and
  stats in this paper is available in the repository under paper/scripts and described in the
  SuP.

  For HGVS resolution, we use ClinVar version 2026-10, and we score converting c.hgvs to
  genomic coords against ClinVar's own VCF coordinate, CROMPOS ref alt.

  We should say all data in this paper can be produced from a reproducible pipeline in the
  repository.

  For a performance benchmark, we run the 500 of the same ClinVar pairs through every backend.
  C.local JSON, REST, UTA, local and remote. The HTTPS engine and sequence source are held
  consistent betweenr runs

  we use two data sets. Our own production search bar collection, which is 32,000 queries and
  private data, and the 998 clean public ClinVar strings with synthetically injected errors.
  3.419 thousand cases in total. This is run through CDOT and the LOVD checker, variant
  validator, and mutilator.

  Transcript version substitution safety is run on consecutive transcript versions from GRCH38
  from a sample of 12,000.

  *Corrections (as dictated):* Results from this paper that involve public data can be
  reproduced from a pipeline in this repository. The invalid HGVs from production searches
  used as examples of real world data cleaning are private because they contain clinical
  information. We have converted this to aggregate counts which are included in the
  repository.

  The HCVSs sent from clinical labs into ShareAnt are also used as a benchmark. Again, we
  cannot provide this as it is private clinical data.

  We use 2 data sets to test cleaning functionality.

  And the synthetic data was also run on third-party cleaners such as LOVD, Variant Validator,
  and Mutilizer.

  remove the such as, and I thought I mentioned the clinvar pairs with VCV XML, just copy
  paste that from when I said it before. [Inserted from I4: "the historical ClinVar
  submissions, the VCV XML, which records the original transcript version that was submitted
  at the time"]

  Also uses the VCV XML as historical data, which contain the original C.HVS with the
  transcript versions provided by labs at submission time

  *Corrections after acceptance (as dictated):* For a performance benchmark, we ran 500 clinvar pairs through the c.local JSON REST and the UTA local. For the UTA remote, we ran only 25 pairs because it was slow.

  *Corrections after acceptance (as dictated):* We check transcript version substitution safety by comparing 12,000 consecutive transcript versions each from both refseq and ensemble GCHR38.

  *Corrections after acceptance (as dictated):* For one, delete the versions. Feel free to change the tense to be consistent. It's a very mechanical change. Feel free to chop words when you think it makes sense.
- **Accepted (your words, cut and joined only):**

  Results from this paper that involve public data can be reproduced from a pipeline in this
  repository (`paper/scripts`, described in the Supplementary Methods). For HGVS resolution,
  we use ClinVar version {{ clinvar_vcf._prov.tools.clinvar }} and the VCV XML as historical
  data, which contain the original c.HGVS with the transcript versions provided by labs at
  submission time. We score converting c.HGVS to genomic coordinates against ClinVar's own VCF
  coordinate (CHROM/POS/REF/ALT). For a performance benchmark, we run
  {{ benchmark.n_pairs | int }} ClinVar pairs through the cdot local JSON, REST and the UTA
  local. For the UTA remote, we run only {{ benchmark.uta_remote_n | int }} pairs because it
  is slow. The HGVS engine and sequence source are held consistent between runs. We use two
  data sets to test cleaning functionality: our own production search bar collection, which is
  {{ cleaning_corpus.corpus_n | commas }} queries and private data, and the
  {{ cleaning.inject_sample_size | commas }} clean public ClinVar strings with synthetically
  injected errors, {{ lovd_comparison.n_cases | commas }} cases in total. The synthetic data
  is also run on third-party cleaners LOVD [@LovdHgvsChecker], VariantValidator
  [@Freeman2018], and Mutalyzer [@Lefter2021]. We check transcript version substitution safety
  by comparing consecutive versions of {{ version_stability.sample_n | commas }} transcripts
  each from both RefSeq and Ensembl GRCh38. The invalid HGVS from production searches, used as
  examples of real-world data cleaning, are private because they contain clinical information.
  The HGVS sent from clinical labs into Shariant [@Tudini2022] are also used as a benchmark.
  We cannot provide this as it is private clinical data. We have converted this to aggregate
  counts which are included in the repository.

---

# Results

**Flow:** coverage (R1) → accuracy on current and historical ClinVar strings, plus Shariant
(R2) → speed (R3) → cleaning, and how other tools compare (R4) → is version substitution safe
(R5).

### R1a. Transcript coverage

- **Job:** Headline coverage figure and where the gain comes from.
- **Opens with:** the unit (one transcript version aligned to one build).
- **Lands on:** T2T, a first for the Python HGVS libraries.
- **Beats:**
  1. Unit definition.
  2. N alignments vs ~N in UTA 2021 [@Hart2020], an N-fold increase.
  3. Gain 1, history: an old NM_ version still resolves after NCBI retires it.
  4. Gain 2, Ensembl: N on GRCh38 alone, none in UTA.
  5. T2T adds N, the first source to bring T2T to the Python HGVS libraries.
- **Must keep:** `coverage.total_count`, `literature.uta_count`, `coverage.improvement_fold`,
  `coverage.ensembl_unique_count`, `coverage.t2t_unique_count`; `[@Hart2020]`. *Checked:* "cdot
  covers … such alignments across all builds and sources", "Ensembl itself, absent from UTA
  entirely".
- **Your version:**

### R1b. Alignment gaps

- **Job:** Coverage includes correct gap information, which another major tool ignores.
- **Opens with:** the JSON stores per-exon alignment gaps, which libraries apply.
- **Lands on:** VEP doesn't apply them, so it can misplace variants after an indel in a gapped
  RefSeq alignment.
- **Must keep:** `[@VepHgvsGaps]`.
- **Your version:**

### R2.1. Current-version ClinVar: parity on RefSeq, Ensembl only in cdot

- **Job:** On current versions, cdot equals UTA for RefSeq, and only cdot does Ensembl.
- **Opens with:** compared against a local UTA (release N), current versions first.
- **Lands on:** UTA can't resolve Ensembl at all (X% vs X%). Note that the aggregate isn't
  reported because the sets are Ensembl-enriched; put that mid-paragraph or as a clause.
- **Beats:**
  1. Local UTA release N. Current-version check first, historical next.
  2. N NM_ + N ENST pairs [@Landrum2025], FASTA sequence so no Ensembl pair is lost to
     SeqRepo.
  3. UTA rate = membership of the exact accession.version in its GRCh38 set (the live database
     agrees).
  4. RefSeq: parity, X% vs X%.
  5. Ensembl: UTA X% vs cdot X%.
- **Must keep:** `clinvar_submitted._prov.tools.uta`, `clinvar.n_refseq`, `clinvar.n_ensembl`,
  `clinvar.cdot_refseq_pct`, `clinvar.uta_refseq_pct`, `clinvar.uta_ensembl_pct`,
  `clinvar.cdot_ensembl_pct`; `[@Landrum2025]`. *Checked:* "no Ensembl pair is dropped", "On
  RefSeq the two backends are at parity", "UTA stores no Ensembl alignments at all".
- **Your version:**

### R2.2. Full-scale ClinVar: almost every resolved variant matches

- **Job:** Correctness at scale.
- **Opens with:** every c.HGVS in variant_summary (N pairs, essentially all RefSeq).
- **Lands on:** what the small residual is (Table S4).
- **Beats:** X% resolve, and X% of those match ClinVar's VCF coordinate. Residual X%: mostly
  indel representation and paralog/copy-number transcripts. Plus N transcripts with unaligned
  bases between exons (Table S3), where cdot usually agrees with UTA and ClinVar differs.
- **Must keep:** `clinvar_vcf.n_pairs`, `clinvar_vcf.resolved_pct`,
  `clinvar_vcf.matched_of_resolved_pct`, `clinvar_vcf.incorrect_pct`,
  `clinvar_vcf_residual.unaligned_between_exons_transcripts`. *Checked:* "resolving every
  c.HGVS in ClinVar's variant_summary … effectively all RefSeq", "is dominated by
  indel-representation differences and by paralog and copy-number transcripts", "map to more
  than one genomic locus".
- **Your version:**

### R2.3. What labs actually submitted is mostly old

- **Job:** Those checks only test current versions. Submitted strings are the real record, and
  they're historical.
- **Opens with:** both corpora so far use ClinVar's recomputed names at the current version,
  so they are a ceiling.
- **Lands on:** cdot holds X% of the superseded versions, and only X% of cited versions are
  missing.
- **Beats:**
  1. Ceiling problem.
  2. Second corpus: submitted (SCV) strings joined to the VCF coordinate (N pairs from N
     expressions; Methods).
  3. All RefSeq: no submitted string cites Ensembl.
  4. X% cite a no-longer-current version (X% weighted by submission). X% cite a transcript no
     longer annotated at any version.
  5. cdot holds X% of those. Only X% of cited versions are absent.
- **Must keep:** `clinvar_submitted.n_unique_pairs`, `clinvar_submitted.n_scv_tx_strings`,
  `clinvar_submitted.version_not_current_pct`,
  `clinvar_submitted.scv_weighted_not_current_pct`, `clinvar_submitted.base_retired_pct`,
  `clinvar_submitted.not_current_in_cdot_pct`, `clinvar_submitted.absent_cdot_pct`.
  *Checked:* "not one submitted string cites an Ensembl transcript", "only … of cited versions
  are absent from its GRCh38 data".
- **Your version:**

### R2.4. The missing versions come from a few labs

- **Job:** The missing versions were never publicly released, which points to substitution as
  the fix.
- **Opens with:** cdot ingests every release, so a missing version should be one never
  released.
- **Lands on:** these are what version substitution (R5) is for.
- **Beats:** only N of N submitters ever cite one. One lab has X% of the absent citations. X%
  of the N absent versions are cited by one submitter alone. This looks like labs aligning
  their own transcripts.
- **Must keep:** the five `submitter_attribution.*` facts. *Checked:* "Because cdot ingests
  every RefSeq and Ensembl annotation release, … submitters only".
- **Your version:**

### R2.5. On submitted strings, cdot beats UTA, mostly on old submissions

- **Job:** The headline accuracy result.
- **Opens with:** cdot X% vs UTA X% on all N submitted pairs.
- **Lands on:** the era pattern: the gap is largest on the oldest submissions.
- **Beats:**
  1. cdot X% vs UTA X%.
  2. The RefSeq gap that was zero on current versions opens to N points: N cdot-only, N
     UTA-only, because UTA lacks X% of the cited versions.
  3. By era (2008-15 / 2016-20 / 2021-26): cdot X/X/X, UTA X/X/X. UTA is lowest on the oldest,
     and lower on the newest than the middle.
- **Must keep:** `clinvar_submitted.full_n`, `clinvar_submitted.cdot_resolved_pct`,
  `clinvar_submitted.uta_resolved_pct`, `clinvar.cdot_refseq_pct`,
  `clinvar_submitted.cdot_only_pct`, `clinvar_submitted.cdot_only`,
  `clinvar_submitted.uta_only`, `clinvar_submitted.uta_no_data_pct`, the six
  `clinvar_submitted.era_*`. *Checked:* "reproducing ClinVar's VCF coordinate for … versus …
  for the same locally loaded UTA", "The RefSeq gap invisible at the current-version ceiling
  (both backends ...)", "because UTA holds no GRCh38 alignment for", "lowest on the oldest
  submissions ... and lower on the newest than on 2016-2020".
- **Note:** six era numbers in one sentence is heavy. Consider giving only the oldest and
  newest in the text.
- **Your version:**

### R2.6. What is left on the submitted corpus

- **Job:** The residual is about version history.
- **Opens with:** X% (N of N) still fail after cleaning and substitution.
- **Lands on:** the rest (bad position or base, grammar, accession absent; Table S8).
- **Beats:** N cite a missing version where substitution declines (safety can't be verified).
  N resolve via the cited version to a coordinate that differs from ClinVar's current one.
  Then the rest.
- **Must keep:** `clinvar_submitted.residual_pct`, `clinvar_submitted.residual_n`,
  `clinvar_submitted.full_n`, `clinvar_submitted_residual.version_refused`,
  `clinvar_submitted_residual.coordinate_drift`. *Checked:* "dominated by historical-version
  effects: … or an accession absent entirely".
- **Your version:**

### R2.7. Shariant [private data]

- **Job:** The same gap on the clinical data that motivated cdot.
- **Opens with:** the same gap holds on Shariant.
- **Lands on:** what the cdot-only strings are (X% old RefSeq versions, X% Ensembl).
- **Beats:**
  1. N unique strings imported into Shariant [@Tudini2022] (N of N transcripts cited at more
     than one version).
  2. Only the transcript layer swapped: cdot X% vs UTA X%.
  3. No ground truth, so this is resolution rate only.
  4. Cdot-only X%: X% RefSeq versions UTA lacks, X% Ensembl.
- **Must keep:** the eight `historical.*` facts; `[@Tudini2022]`. *Checked:* "With only the
  transcript-data layer swapped … versus … for UTA", "over many years ( … of the … distinct
  transcripts are cited at more than", "UTA holds no GRCh38 alignment and".
- **Your version:**

### R3.0. Table 1 caption (optional)

- **Must keep:** all `benchmark.*` in table and caption. *Checked:* "HGVS/s over … timed
  passes".

### R3.1. cdot is faster than UTA, locally and remotely

- **Job:** Interpret Table 1.
- **Opens with:** in both like-for-like comparisons cdot is faster.
- **Lands on:** after one prefetch, REST is as fast as local JSON, because both are then
  limited by the engine and sequence layer.
- **Beats:**
  1. Local: UTA N HGVS/s, about a quarter of local JSON (several SQL queries per lookup vs one
     JSON object).
  2. Remote: cdot REST (no prefetch) is orders of magnitude faster than public UTA, which pays
     wide-area round trips on every lookup.
  3. One `prefetch()` batch (under a second, untimed) brings REST within X% of local JSON.
     Both are then limited by the shared engine and sequence layer.
- **Must keep:** `benchmark.uta_local_tps`, `benchmark.cdot_rest_vs_uta_remote_ratio`,
  `benchmark.rest_prefetch_vs_local_diff_pct`. *Checked:* "about a quarter of local-JSON
  throughput".
- **Your version:**

### R3.2. All of ClinVar in minutes

- **Job:** It scales.
- **Opens with:** a GRCh38 file loads in ~N s, and one process then resolves all of ClinVar.
- **Lands on:** remote UTA would take about N years.
- **Beats:** N pairs in ~N min, X% resolved. N hot repeats: N HGVS/s (range), within X% of
  cold (N), so the cache isn't the limit. Remote UTA extrapolates to ~N years.
- **Must keep:** `benchmark.grch38_load_time_s`, all `benchmark_fullscale.*`.
- **Your version:**

### R4.1. Cleaning rescues a large share of real queries [private data]

- **Job:** The main cleaning result.
- **Opens with:** production search-box queries (N) from VariantGrid-based platforms.
- **Lands on:** parseable rose from X% to X%, N rescued (~X% of failures), zero regressions.
- **Beats:**
  1. N real queries [@VariantGrid] matching a broad HGVS regex.
  2. What people pasted or typed, with damage from the route: Word/PDF whitespace and odd
     characters, casing, punctuation, trailing protein.
  3. The headline numbers.
- **Must keep:** `cleaning_corpus.corpus_n`, `cleaning_corpus.as_submitted_pct`,
  `cleaning_corpus.after_pct`, `cleaning_corpus.rescued`,
  `cleaning_corpus.rescued_share_pct`, `cleaning_corpus.regressions`; `[@VariantGrid]`.
  *Checked:* "failed as-submitted) with … regressions", "strings rescued, about … of those
  that failed as-submitted".
- **Your version:**

### R4.2. Which fixes did the work (Table 2)

- **Job:** Breakdown by fix type.
- **Opens with:** punctuation and whitespace dominate (Table 2).
- **Lands on:** the rarest fix, restoring a dropped prefix, needs local transcript data, which
  is something a string-only checker can't do.
- **Beats:** then gene/transcript wrapper and structure rebuilds, then the prefix point.
- **Must keep:** nothing fact-bound.
- **Your version:**

### R4.2t. Table 2 caption (optional)

- **Must keep:** `cleaning_corpus.corpus_n`, `cleaning_corpus.rescued`, the include-csv block.

### R4.3. The ceiling of string repair [private data]

- **Job:** What cleaning can't fix.
- **Opens with:** X% (N) of genuine HGVS-shaped queries still fail, which is the ceiling.
- **Lands on:** ~X% is valid HGVS the biocommons grammar rejects.
- **Beats:** N non-HGVS queries (URLs, templates, prose) excluded, as a clause. Of the rest
  (Table S6): well over half incomplete or missing a reference (can't be invented), about a
  quarter fixable in principle (future rules), then the grammar gap.
- **Must keep:** `cleaning_corpus.nonhgvs_n`, `cleaning_corpus.residual_pct`,
  `cleaning_corpus.residual_n`, `cleaning_corpus.residual_grammar_gap_pct`. *Checked:*
  "Classified under a fixed decision-tree taxonomy … about a quarter is in principle
  fixable".
- **Your version:**

### R4.4. Two corpora, two ways to fail

- **Job:** Tie R2 and R4 together.
- **Opens with:** the public submitted corpus behaves quite differently.
- **Lands on:** database submissions fail by version age, typed input fails by formatting.
- **Beats:** On submitted strings, `fix_hgvs()` repaired only N (N regressions), because they
  are mostly well-formed. Their failures are old versions (R2.5, R2.6).
- **Must keep:** `clinvar_submitted.rescued_by_fix`, `clinvar_submitted.regressions`.
- **Your version:**

### R4.5. Injection control

- **Job:** A reproducible no-regression check.
- **Opens with:** `inject_and_clean.py` injects each fix category into clean ClinVar strings.
- **Lands on:** its purpose is the no-regression guarantee. Recovery is in the supplement only,
  since the test injects the very errors it fixes.
- **Must keep:** `cleaning.inject_regressions`.
- **Your version:**

### R4.6. How other tools do on the same cases

- **Job:** Comparison numbers (moved from the Discussion).
- **Opens with:** the same injected cases were run through LOVD, VariantValidator and
  Mutalyzer (Methods).
- **Lands on:** no tool falsely corrected a valid input.
- **Beats:**
  1. Caveat first: the corpus is `clean_hgvs()`'s own targets, so it recovers all of them by
     construction. The test is how much of that ground each tool also covers.
  2. LOVD top-1 (N cases): X% weighted by production error mix (X% unweighted). It matches on
     single-token errors, not free-text damage or cross-family repairs.
  3. VariantValidator X%, Mutalyzer X% weighted (Table S5). Each repairs what it can interpret
     and rejects the rest.
  4. No tool altered a valid input.
- **Must keep:** `lovd_comparison.n_cases`, `lovd_comparison.lovd_version`,
  `lovd_comparison.lovd_top1_weighted_pct`, `lovd_comparison.lovd_top1_pct`,
  `vv_mutalyzer_comparison.vv_weighted_pct`, `vv_mutalyzer_comparison.mut_weighted_pct`.
  *Checked:* "neither tool altered any valid input", "neither falsely corrected a valid
  input" (these can merge into one claim).
- **Your version:**

### R5.1. How often version substitution is safe

- **Job:** Quantify the risk.
- **Opens with:** substituting a version can move a variant's coordinate, and biocommons/hgvs
  has no substitution with any provider.
- **Lands on:** per coding base, a random variant is safe X% (RefSeq) / X% (Ensembl) of the
  time.
- **Beats:** consecutive pairs: RefSeq N, X% preserve every coding coordinate. Ensembl N, X%.
  Then the per-base figures.
- **Must keep:** `version_stability.refseq_pairs`, `version_stability.refseq_preserving_pct`,
  `version_stability.ensembl_pairs`, `version_stability.ensembl_preserving_pct`,
  `version_stability.refseq_pervariant_safety`,
  `version_stability.ensembl_pervariant_safety`. *Checked:* "preserved every coding".
- **Your version:**

### R5.2. What kind of drift

- **Job:** When it moves, it's usually all or nothing. The partial case is rare but
  positional.
- **Opens with:** when a coordinate moves, it's usually the whole CDS.
- **Lands on:** in partial drift, risk grows toward the 3' end (Fig S1).
- **Beats:** whole-CDS drift X%/X% (re-annotation of the coding region). Partial drift, the
  dangerous case, is rare (X%/X%). It keeps a 5' prefix up to the first change, so preservation
  falls from X% (5' decile) to X% (3') for RefSeq (X% to X% Ensembl).
- **Must keep:** `version_stability.refseq_full_drift_pct`,
  `version_stability.ensembl_full_drift_pct`,
  `version_stability.refseq_partial_drift_pct`,
  `version_stability.ensembl_partial_drift_pct`, the four `positional_drift.*_decile*_pct`.
  *Checked:* "the most dangerous case, a partial substitution … is rare", "Within that tail the
  risk is positional: … so preservation falls".
- **Your version:**

### R5.3. Drift is concentrated

- **Job:** The risky transcripts are a short list.
- **Opens with:** few accessions ever drift.
- **Lands on:** short and identifiable enough to flag or blocklist.
- **Beats:** only X% of RefSeq (X% Ensembl) accessions ever drift. The top 5% hold X% (X%) of
  drifting pairs.
- **Must keep:** `version_stability.refseq_accessions_drift_pct`,
  `version_stability.ensembl_accessions_drift_pct`,
  `version_stability.refseq_top5_drift_share_pct`,
  `version_stability.ensembl_top5_drift_share_pct`. *Checked:* "Drift is also concentrated:
  only … so the risky transcripts form a short, identifiable list", "the 5% of accessions with
  the most drifting pairs".
- **Your version:**

### R5.4. The safety check works across builds and makes no errors on ClinVar

- **Job:** Validate the check from M4.
- **Opens with:** the check works even when the version was never aligned to the target build.
- **Lands on:** every accepted substitution that resolves gives the same coordinate (0
  changes).
- **Beats:** structure is identical across builds for X% RefSeq, X% Ensembl versions. On
  ClinVar: N accepted substitutions, N resolve, N coordinate changes.
- **Must keep:** `version_stability.refseq_struct_portable_pct`,
  `version_stability.ensembl_struct_portable_pct`,
  `version_safety_validation.n_safe_substitutions`,
  `version_safety_validation.n_resolved_under_w`,
  `version_safety_validation.n_coordinate_changes`. *Checked:* "was never aligned", "every one
  that resolves under the substitute".
- **Your version:**

---

# Discussion

**Flow:** what we showed (D1) → others already use it as a general source (D2) → other tools
that repair HGVS (D3) → how they compare and why they differ (D4) → limitations (D5) → the
design principle to take away (D6).

### D1. What the results show

- **Job:** Open the Discussion on findings rather than a product description.
- **Opens with:** the main result: on what labs actually submitted, cdot resolves far more than
  UTA, because of history depth.
- **Lands on:** it is a general source of versioned alignments, not tied to Python or
  biocommons/hgvs. This leads into D2.
- **Beats:**
  1. Historical strings: the gain is in old versions (R2).
  2. Speed: fast enough locally and over REST that no database is needed (R3).
  3. Cleaning: rescues typed input without breaking valid input (R4).
  4. The data itself: RefSeq + Ensembl, three builds, full history, JSON readable from any
     language. The per-release JSON can replace a GTF/GFF, and REST serves thin clients.
- **Must keep:** no facts needed. Refer back to results rather than repeating numbers.
- **Your version:**

### D2. Use beyond Python

- **Job:** Evidence for "general source".
- **Opens with:** others already use it this way.
- **Lands on:** the Rust projects read cdot JSON directly with only a JSON parser.
- **Beats:** GitHub search (date N): N projects outside our group (a lower bound). N use the
  Python providers, including MaveDB [@Rubin2025]. N are Rust: mehari [@Hartmann2026], hgvs-rs
  [@HgvsRs], ferro-hgvs [@FerroHgvs].
- **Must keep:** `adoption.census_date`, `adoption.n_projects`,
  `adoption.n_python_projects`, `adoption.n_rust_projects`; the four citations. *Checked:* "The
  other … are written in Rust".
- **Your version:**

### D3. Other tools that repair HGVS

- **Job:** Related work, and where `clean_hgvs()` fits.
- **Opens with:** cdot isn't the first to repair HGVS.
- **Lands on:** `clean_hgvs()` runs offline before parsing, reports every change, and never
  breaks a string that already parsed.
- **Beats:** VariantValidator [@Freeman2018] (biocommons + own UTA) corrects what it can
  interpret. Mutalyzer (I1b already gives its number) [@Lefter2021]. LOVD [@LovdHgvsChecker]
  suggests ranked corrections without a reference. ClinGen Allele Registry [@Pawliczek2018]
  normalises on registration. In all of these, repair happens inside validation or
  registration, usually behind a web service.
- **Must keep:** `[@Freeman2018]`, `[@Lefter2021]`, `[@LovdHgvsChecker]`, `[@Pawliczek2018]`.
  `literature.mutalyzer_autocorrect_pct` can be dropped here, since it is in I1b.
- **Your version:**

### D4. Why the tools score differently, and where they're heading

- **Job:** Interpret R4.6.
- **Opens with:** VariantValidator and Mutalyzer validate against the reference, so an
  unparseable string is an error for them to report.
- **Lands on:** the LOVD figure also describes the syntax-repair layer of current
  VariantValidator.
- **Beats:** LOVD is the closest comparator (offline, no reference). Since VV 4.0.0 (June 2026)
  it embeds the LOVD checker [@VariantValidator400] and returns its suggestions with the
  error, but doesn't apply them. cdot supplies the coordinate layer these services sit on, so
  it complements them.
- **Must keep:** `[@VariantValidator400]`. *Checked:* "offline, no reference sequence".
- **Your version:**

### D5. Limitations

- **Job:** Honest limits.
- **Opens with:** the main limit is currency.
- **Lands on:** ClinVar's submitter bias (Table S4), so the scale check is not an unbiased
  sample.
- **Beats:** a version released after the last ingest isn't covered until regeneration
  (automating that is planned). Some residue is valid HGVS biocommons doesn't parse (repeats,
  alleles; Tables S6, S8). ClinVar submitter bias toward current RefSeq.
- **Must keep:** nothing fact-bound.
- **Your version:**

### D6. Closing: safe by default, heuristics opt-in

- **Job:** The takeaway.
- **Opens with:** cdot separates unambiguous cleaning (automatic) from heuristics that can be
  wrong.
- **Lands on:** your call. Possibly what this makes possible, eg pooled clinical histories
  like Shariant become usable.
- **Beats:** adjacent version and canonical transcript are opt-in and always reported as an
  `HGVSFix` the caller can inspect or reject.
- **Must keep:** nothing.
- **Your version:**

---

# Abstract (do last)

### A1. Motivation

- **Beats:** Shariant pools lab classifications against whatever version each lab used.
  Making it usable means resolving as many real-world strings as possible, including
  malformed ones and old versions. UTA means either a local PostgreSQL or a slow,
  often-firewalled public server, has ~N alignments (2021), limited history, no Ensembl.
- **Must keep:** `literature.uta_count`; `[@Tudini2022]`.
- **Your version:**

### A2. Results

- **Beats:** N alignments, RefSeq + Ensembl, three builds (first T2T for the Python
  libraries), JSON.gz or REST, no database. ClinVar submitted strings (X% non-current): cdot X%
  vs UTA X%. ~N HGVS/s in memory, Nx local UTA. REST orders of magnitude faster than public
  UTA. Cleaning: add the R4.1 rescue figure. Opt-in substitution only when the safety check
  passes.
- **Must keep:** `coverage.total_count`, `clinvar_submitted.version_not_current_pct`,
  `clinvar_submitted.cdot_resolved_pct`, `clinvar_submitted.uta_resolved_pct`,
  `benchmark.cdot_local_tps`, `benchmark.cdot_local_vs_uta_local_ratio`,
  `benchmark.cdot_rest_vs_uta_remote_ratio`. *Checked:* "with no database required", "cdot
  reproduced ClinVar's genomic coordinate for … versus … for UTA".
- **Your version:**
