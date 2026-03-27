# Figures

*Application Note target: 1–2 main figures. Additional material goes to Supplementary.*
*Each legend must be self-contained.*

---

## Figure 1 — Transcript coverage and ClinVar resolution

**Type**: Two-panel figure

**Panel A** — Grouped bar chart: transcript counts by source (cdot RefSeq, cdot Ensembl, UTA) × genome build (GRCh37, GRCh38, T2T). Log scale on Y-axis. Annotate with "9×" arrow from UTA to cdot combined.

**Panel B** — Stacked bar chart or paired bars: ClinVar HGVS resolution rate for cdot vs UTA. Break the cdot bar into: resolved by RefSeq GRCh38 / resolved by Ensembl GRCh38 / resolved by GRCh37 / unresolved. Show UTA as comparison.

**Legend**: `Transcript coverage and ClinVar HGVS resolution. (A) Versioned transcript alignments in cdot (RefSeq and Ensembl) and UTA across genome builds. cdot provides {{ coverage.total_count | commas }} alignments, {{ coverage.improvement_fold | fmt('.0f') }}× more than UTA uta_20210129. T2T-CHM13v2.0 is absent from UTA. (B) Resolution of {{ clinvar.n_variants | commas }} ClinVar HGVS variant descriptions. cdot resolves {{ clinvar.cdot_resolution_pct | pct(1) }} versus {{ clinvar.uta_resolution_pct | pct(1) }} for UTA, with additional coverage from Ensembl and historical RefSeq annotation releases.`

*This is the single most important figure — it encapsulates both the coverage and practical utility arguments.*

---

## Figure 2 — Architecture (optional; use if word count allows)

**Type**: Schematic

**Content**: Compact diagram showing:
- Left: RefSeq GFF3 + Ensembl GTF → Snakemake pipeline → JSON.gz
- Centre: JSON.gz ↔ cdotlib.org REST API
- Right: biocommons/hgvs client and PyHGVS client, each reading from local JSON or REST

**Legend**: `cdot architecture. RefSeq and Ensembl annotation files are processed into versioned JSON.gz files accessible locally or via the cdotlib.org REST API. Python clients integrate with biocommons/hgvs and PyHGVS.`

*Include only if Figure 1 alone leaves the paper feeling abstract. If word count is tight, drop this and describe the architecture in text.*

---

## Supplementary Figures

**Figure S1** — Speed benchmark: throughput (transcripts/second) for cdot local, cdot REST, UTA local, UTA remote. Log-scale bar chart.

**Figure S2** — Historical RefSeq coverage: cumulative unique NM_ accession versions by annotation release year. Motivates multi-release ingestion.

**Figure S3** — T2T unique transcripts: genes present in T2T-CHM13v2.0 but absent from GRCh37/GRCh38.

---

*Notes:*
- *Figure 1 requires ClinVar benchmark data (#5) and Snakemake summary stats*
- *Speed numbers for Figure S1 need a scripted benchmark (commit the script)*
- *Consider combining panels A and B into one clean figure with a shared colour palette*
