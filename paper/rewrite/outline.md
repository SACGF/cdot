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
12. **Results restructure (2026-11-26).** Results cut from 20 paragraphs to 7, one central
    result about cdot each. R1b cut. Details not about cdot move to a new Supplementary
    Results section. See "Restructure" under Results.

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

  *Informality pass (as dictated, 2026-11-24):* The C. form is typically preferred in clinical
  genetics because medical scientists find it easier to reason about

  In genomics, we need to convert them to genomic coordinates.

  To do this, we require the reference transcript to be at a specific version.

  The exon coordinates from the transcript genome alignments can be used to convert them to
  genomic coordinates.
- **Accepted after informality pass:** HGVS is the clinical standard for describing variants [@DenDunnen2016]. The c. form is typically preferred in clinical genetics because medical scientists find it easier to reason about. In genomics, we need to convert them to genomic coordinates. To do this, we require the reference transcript to be at a specific version. These versions of the transcript are aligned to the genome by RefSeq and Ensembl. The exon coordinates from the transcript genome alignments can be used to convert them to genomic coordinates.

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

  *Informality pass (as dictated, 2026-11-24):* Biocommons abstracts transcript data retrieval
  inside the data provider interface, which has become a standard being adopted by a number of
  other tools, including Variant Validator and MaveDB.

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

  *Informality pass (as dictated, 2026-11-24):* UTA's popularity demonstrates it is an excellent
  resource, though the Postgres implementation has a number of drawbacks.

  Performing their own alignments means that the process is lengthy and complex, and  Ensemble
  is not supported, and many historical refseq transcripts are missing.

  *Mechanical fixes he approved:* deleted "Postgres" from "the Postgres implementation",
  "their" → "its", dropped the first "and" before "Ensembl", and cut "RefSeq only, and it
  provides a limited amount of historical support" from the UTA count sentence (repeated).

  *Informality pass (as dictated):* A drawback of Postgres database is that setting up a local
  instance requires systems administration capability.

  *Informality pass, I3 (as dictated):* These often went back many years.

   One of the most critical problems in Sharint was how to resolve historical Hgvs that were
  often many years old.

  *Informality pass, I3 search box and three problems (as dictated):* People were searching
  using HTVSs with many subtle and different problems, and we tried to resolve these to what
  the users actually intended.

  We had to solve a number of problems. Cleaning invalid HTTPSs, collecting as many historical
  transcript versions as possible, and conversion between different genome builds.

  As not all transcript versions were available in all builds, we attempted to find the nearest
  transcript version available for that build that had the same exonic structure and was safe
  to substitute. Brackets opt-in and flagged.

  *Cut:* "There were often no official releases of a transcript version against a particular
  genome build." (now said by "As not all transcript versions were available in all builds").

  *Informality pass, I4 (as dictated):* Almost all languages have JSON libraries that have been
  highly optimized.

  JSON is also the standard protocol for transferring data across the internet via a REST
  server, which means it was easy to provide a web server which can be accessed through most
  workplace Firewalls.

  We provide our cleaning library, which fixes many common mistakes and returns a record of
  what fixes were made.

  *Approved fix:* "standard protocol" → "standard format" (JSON is a format; HTTP is the protocol).

  *Informality pass, I4 last paragraph (as dictated):* We validate and benchmark CDOT by
  comparing it against UTA on ClinVar VCV XML, which records the original transcript version
  that was submitted by the lab at the time, giving a historical test record.

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

  *Informality pass, M3 (as dictated):* You call the function with an HGVS string and it returns
  the repaired string plus a list of fixes applied.

  For the last one, change if you Pass to passing a [applied to both "If you pass" and "If you
  also pass"]

  Change you Call the function to calling the function. [Cut "and it" to make it grammatical;
  cut ", it" from both "Passing a ..., it" sentences to fix the dangling modifier.]

  What cannot be fixed is flagged as an error

  The pure string conversion function is exposed as clean HTTPS, which doesn't load any data.
  It is called by fixHTTPS. ["doesn't" → "does not"]

  When the requested transcript version is not available, fix HTTPS can substitute another
  version.

  The function accepts a parameter version fallback, where the caller sets the resolution
  strategy, such as trying higher versions then lower, closest or the latest.

  It is returned as an HGVS fix, so the user has visibility on what operations were performed.

  Checks are made to ensure the substituted transcript version has an identical CDS structure,
  the CDs length plus the length of coding exons.

  All string modifications return an accompanying HGVS fix, so the user has visibility on what
  operations were performed.

  [Cut from M4 at his request, as it repeats M3's "Every change is recorded as an `HGVSFix`".]
  [Then, at his request, M3's sentence was replaced by joining the two: his newer wording plus
  M3's "which has severity, code, and before/after".]

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

  *Informality pass, M4 (as dictated, 2026-11-26):* for the requested genome build [replaces
  "for the genome build you request"]

  The string is left unchanged if the version substitution is unsafe, and a HGVSFix error is
  returned

  It is possible to override this to return a changed string and a warning that the
  coordinates may have moved. Feel free to cut repetition if you think it is unecessary

  [Cut, at his go-ahead: "It's safe to swap only if the structures match." (repeats the CDS
  structure check sentence).]
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

  *Informality pass, M5 (as dictated, 2026-11-26):* "batches them into" [replaces "batches
  things into"]

  biocommons uses the SeqFetcher interface for sequence retrieval and provides a SeqRepo
  implementation [@Hart2020]. cdot implement a SeqFetcher (FastaSeqFetcher) that can load a
  local FASTA file, which pulls out the exon sequences and pastes them together to make a
  transcript. ["implement" → "implements"]

  genome, but is not guaranteed to be safe for RefSeq [replaces "genome, but might have some
  troubles with RefSeq"]

  We also provide ChainedSeqFetcher which can configured with other SeqFetchers and take the
  first that succeeds. This can be used to prefer SeqRepo sequences when available, and fall
  back on the Fasta file. ["can configured" → "can be configured"; replaces the last two
  sentences]

  *Informality pass, M6 (as dictated, 2026-11-26):* We ran only { benchmark } pairs on the
  public UTA instance for performance reasons [replaces "For the UTA remote, we run only N
  pairs because it is slow."]
  [At his request, "We ran" → "We run" to match the present tense of the rest of M6.]
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

**Flow:** coverage (R1) → correctness on ClinVar (R2) → old transcript versions in submitted
strings and Shariant (R3) → speed (R4) → cleaning (R5) → other tools (R6) → version
substitution, risk and check (R7).

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
- **Your version (as dictated):** The following counts refer to a transcript version aligned to
  a particular genome build. C.Dot has 1.9 million alignments compared to 300,000 in the 2024
  UTA release, 6.5 times as many.

  This gain comes from ingesting every archived refseq release, and also from bringing in
  ensemble transcripts which UTA lacks. Bringing in Ensemble GRCH30 alone accounts for 850,000
  transcript versions. adding T to T brings in 186,000.

  *Corrections (as dictated):* The following counts are of a transcript version aligned to a
  particular genome build.
- **Accepted (your words, cut and joined only):**

  The following counts are of a transcript version aligned to a particular genome build.
  cdot has {{ coverage.total_count | commas }} alignments compared to
  {{ coverage.uta_count | commas }} in the 2024 UTA release [@Hart2020],
  {{ coverage.improvement_fold | dp(1) }} times as many. This gain comes from ingesting every
  archived RefSeq release, and also from bringing in Ensembl transcripts which UTA lacks.
  Ensembl GRCh38 alone accounts for {{ coverage.ensembl_unique_count | commas }} transcript
  versions. Adding T2T brings in {{ coverage.t2t_unique_count | commas }}.

### R1b. Alignment gaps

- **Job:** Coverage includes correct gap information, which another major tool ignores.
- **Opens with:** the JSON stores per-exon alignment gaps, which libraries apply.
- **Lands on:** VEP doesn't apply them, so it can misplace variants after an indel in a gapped
  RefSeq alignment.
- **Must keep:** `[@VepHgvsGaps]`.
- **Cut (author's decision, 2026-11-26).** Not a measured result, and M2 already says the
  exon data has a `gap` field applied during conversion. When Results goes in, drop the VEP
  sentence from R1 in `results.md`. The `VepHgvsGaps` bib entry can stay (unused entries do
  not render).

### Restructure (agreed 2026-11-26)

Each paragraph has one central result about cdot (see the runbook). Results goes from 20
paragraphs to 7: R1 (coverage, accepted above as R1a), then R2 to R7 below. The "interesting
things we found" that are not about cdot move to a new **Supplementary Results** section in
`supplementary.md`, and the main text points to it. Make each move when Results goes in.

**Moving to Supplementary Results** (old entry → what moves; the old wording is in
`results.md` until then):

- Old R2.1: the current-version check method (UTA rate is membership of the exact
  accession.version in its GRCh38 set, the live database agrees; the pair sets are
  Ensembl-enriched so no aggregate is reported) and RefSeq parity on 500 pairs (100% vs
  100%). Its Ensembl number stays in R2.
- Old R2.2: the residual breakdown (indel representation, paralog/copy-number transcripts,
  five transcripts with bases unaligned between exons). Already Table S4, so just a pointer.
- Old R2.3: weighted-by-submission (69.3%), retired transcripts (0.7%), absent from cdot
  (0.6%), and "no submitted string cites Ensembl".
- Old R2.4: submitter attribution (88 of 3,028 submitters; one lab has 94.8% of absent
  citations; 82.9% of absent versions cited by one submitter only; looks like in-house
  alignments).
- Old R2.5: the era split (2008-15 / 2016-20 / 2021-26) and the cdot-only / UTA-only counts
  (562,173 / 87). The era split tracks submission date, not transcript-version age, and most
  of the gain is in 2021-26 submissions that cite old versions.
- Old R2.6: the submitted-corpus residual (0.7%; 17,667 refused substitutions, 546 drift).
  Already Table S8.
- Old R2.7: the Shariant breakdown (1,654 of 3,731 transcripts cited at more than one version;
  cdot-only 98.6% old RefSeq versions, 1.3% Ensembl; no ground truth). The headline stays in
  R3 as one sentence.
- Old R3.2: the hot/cold cache comparison. **Delete** rather than move.
- Old R4.2: the fix-type walkthrough. Table 2 stays in the main text; only the dropped-prefix
  sentence stays in R5.
- Old R4.3: the cleaning ceiling (887 queries, 2.7%: 63% incomplete, 28% fixable, 9% grammar
  gap; 81 non-HGVS excluded). Already Table S6, so just a pointer.
- Old R4.4: "database submissions fail by version age, typed input by formatting" (1,468
  repaired on the submitted corpus, 0 regressions). Moves to **Discussion D1**, not the
  supplement.
- Old R4.5: the injection control as a no-regression check. Already Table S5. R6 introduces
  the injected corpus.
- Old R5.2: whole-CDS vs partial drift, and the 5'→3' decline. Already Figure S1.
- Old R5.3: drift is concentrated (5.5% of RefSeq accessions ever drift; the top 5% hold 92.7%
  of drifting pairs).

### R2. cdot reproduces ClinVar's coordinates

- **One-line result:** across every c.HGVS in ClinVar, cdot reproduces ClinVar's own genomic
  coordinate for 99.95% of the variants it resolves.
- **Opens with:** every c.HGVS in ClinVar's variant_summary ({{ clinvar_vcf.n_pairs }} =
  4,583,798 pairs, essentially all RefSeq), scored against ClinVar's VCF coordinate.
- **Lands on:** the residual is small and listed in Table S4 (one clause), or the Ensembl
  clause.
- **Beats:**
  1. 99.46% resolve (`clinvar_vcf.resolved_pct`), 99.954% of those match
     (`clinvar_vcf.matched_of_resolved_pct`). Residual 0.05% (`clinvar_vcf.incorrect_pct`),
     Table S4.
  2. Ensembl, which UTA can't do at all: of 818 Ensembl pairs (`clinvar.n_ensembl`), cdot
     resolves 99.8% (`clinvar.cdot_ensembl_pct`), UTA 0% (`clinvar.uta_ensembl_pct`).
     [@Landrum2025].
- **Must keep:** the facts above; `[@Landrum2025]`. *Checked:* "effectively all RefSeq" (an
  [[ok]]), "UTA stores no Ensembl alignments" (`clinvar.uta_ensembl_pct == 0`).
- **Your version:**

### R3. Labs cite old transcript versions, which cdot has and UTA doesn't

- **One-line result:** most transcript versions labs submitted to ClinVar are superseded, so
  cdot resolves 99.3% of submitted strings against 81.9% for UTA. The same holds on Shariant's
  clinical data.
- **Opens with:** the corpus in R2 uses ClinVar's recomputed names, always at the current
  version. The strings labs actually submitted (from the VCV XML, Methods) are the historical
  record: 3,245,277 pairs (`clinvar_submitted.full_n`).
- **Lands on:** Shariant, one sentence: 36,010 strings (`historical.n_lines`) from clinical
  labs [@Tudini2022], cdot 99.6% vs UTA 93.7% (`historical.cdot_resolved_pct`,
  `historical.uta_resolved_pct`). Details in the supplement.
- **Beats:**
  1. 74.7% cite a version that is no longer current (`clinvar_submitted.version_not_current_pct`),
     and cdot holds 99.3% of those (`clinvar_submitted.not_current_in_cdot_pct`).
  2. cdot 99.3% vs UTA 81.9% against a local UTA release
     (`clinvar_submitted._prov.tools.uta`), because UTA has no alignment for 17.9% of the cited
     versions (`clinvar_submitted.uta_no_data_pct`).
  3. Shariant sentence.
- **Must keep:** the facts above; `[@Tudini2022]`. *Checked:* cdot "versus" UTA on submitted
  pairs, and on Shariant.
- **Watch:** 99.3% appears twice (versions held, and strings resolved). Say them so they
  aren't confused.
- **Your version:**

### R4. cdot is faster than UTA

- **One-line result:** cdot is about 4× faster than UTA locally and about 830× faster
  remotely; all of ClinVar takes 2 hours, where the public UTA would take over a year.
- **Opens with:** Table 1, same 500 pairs, engine and sequence layer held constant.
- **Lands on:** all of ClinVar: 4,583,798 pairs in about 121 minutes
  (`benchmark_fullscale.wall_min`); public UTA extrapolates to over a year
  (`benchmark_fullscale.uta_remote_days` = 408 days).
- **Beats:**
  1. Local: cdot JSON 741 HGVS/s vs UTA 187 (`benchmark.uta_local_tps`, about a quarter).
  2. Remote: cdot REST 108 vs public UTA 0.13, about 830× (`benchmark.cdot_rest_vs_uta_remote_ratio`).
  3. After one `prefetch()`, REST is within 0.7% of local JSON
     (`benchmark.rest_prefetch_vs_local_diff_pct`).
  4. GRCh38 file loads in about 5 s (`benchmark.grch38_load_time_s`), then all of ClinVar.
- **Must keep:** the facts above. *Checked:* "about a quarter of local-JSON throughput".
- **Overlap:** I2b already gives UTA remote at 0.13 HGVS/s. Don't repeat it, use the ratio.
- **Your version:**

### R5. Cleaning makes most failing real-world queries usable [private data]

- **One-line result:** on real search-box queries, cleaning makes 67% of the failing strings
  parseable (91.7% → 97.3%) with no regressions.
- **Opens with:** 32,671 real queries (`cleaning_corpus.corpus_n`) from VariantGrid-based
  platforms [@VariantGrid], whatever people pasted or typed.
- **Lands on:** the rarest fix, restoring a dropped accession prefix, needs the transcript
  data, which a string-only checker doesn't have (Table 2). Or the Table S6 pointer for what's
  left.
- **Beats:**
  1. The corpus, and briefly what the damage looks like (Word/PDF whitespace, casing,
     punctuation).
  2. 91.7% → 97.3% (`cleaning_corpus.as_submitted_pct`, `after_pct`), 1,828 rescued
     (`rescued`), 67% of failures (`rescued_share_pct`), 0 regressions (`regressions`).
  3. Table 2: punctuation and whitespace dominate. The dropped-prefix point.
  4. Optional clause: what still fails is mostly incomplete input (Table S6).
- **Must keep:** the facts above; `[@VariantGrid]`. *Checked:* "rescued, about … of those that
  failed".
- **Your version:**

### R6. Other tools repair only part of the same errors

- **One-line result:** on the same injected errors, LOVD repairs 69%, VariantValidator 54%
  and Mutalyzer 34% (weighted by the production error mix), and none of them changes a valid
  input.
- **Opens with:** the injected corpus: 998 clean ClinVar strings with each cleaning error
  injected, 3,419 cases (`lovd_comparison.n_cases`), run through LOVD
  (`lovd_comparison.lovd_version`), VariantValidator and Mutalyzer
  (`vv_mutalyzer_comparison.service_date`).
- **Lands on:** no tool falsely corrected a valid input.
- **Beats:**
  1. Caveat early: these are `clean_hgvs()`'s own targets, so it gets 100% by construction.
     The question is how much of that ground each tool also covers.
  2. LOVD 69% weighted (34% unweighted), VariantValidator 54%, Mutalyzer 34% (Table S5).
  3. No false corrections (`lovd_false_corrections`, `vv_false_corrections`,
     `mut_false_corrections` all 0).
- **Must keep:** `lovd_comparison.n_cases`, `lovd_comparison.lovd_version`,
  `lovd_comparison.lovd_top1_weighted_pct`, `lovd_comparison.lovd_top1_pct`,
  `vv_mutalyzer_comparison.vv_weighted_pct`, `vv_mutalyzer_comparison.mut_weighted_pct`;
  `[@LovdHgvsChecker]`, `[@Freeman2018]`, `[@Lefter2021]`. *Checked:* "no tool falsely
  corrected a valid input".
- **Moves from:** the Discussion (D3/D4 keep only what the numbers mean).
- **Your version:**

### R7. Version substitution: risky without a check, safe with cdot's

- **One-line result:** swapping to a neighbouring version moves coordinates for about 4% of
  RefSeq version pairs, and cdot's structure check accepted 8.4 million ClinVar substitutions
  with no coordinate changes.
- **Opens with:** the risk: of 19,393 consecutive RefSeq version pairs
  (`version_stability.refseq_pairs`), 96.0% keep every coding coordinate
  (`refseq_preserving_pct`); Ensembl 98.0% of 12,377.
- **Lands on:** 8,411,749 accepted substitutions (`version_safety_validation.n_safe_substitutions`),
  0 coordinate changes (`n_coordinate_changes`).
- **Beats:**
  1. The risk (one measure, per pair; the per-base figures go to the supplement).
  2. The check works across builds: structure identical on GRCh37/GRCh38/T2T for 99.4% of
     RefSeq and 100% of Ensembl versions (`refseq_struct_portable_pct`,
     `ensembl_struct_portable_pct`).
  3. On ClinVar: 8,411,749 accepted, 8,411,000 resolve under the substitute
     (`n_resolved_under_w`), 0 coordinate changes. Figure S1 pointer for how drift happens.
- **Must keep:** the facts above. *Checked:* "every one that resolves under the substitute",
  "was never aligned".
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
  1. Historical strings: the gain is in old transcript versions (R3).
  2. Speed: fast enough locally and over REST that no database is needed (R4).
  3. Cleaning: rescues typed input without breaking valid input (R5). Moved in from old R4.4:
     the two corpora fail differently. Database submissions are well-formed and fail on old
     transcript versions (cleaning repaired only 1,468 submitted strings, 0 regressions,
     `clinvar_submitted.rescued_by_fix`, `clinvar_submitted.regressions`); typed input fails
     on formatting.
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
