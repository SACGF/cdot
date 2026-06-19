# cdot: resolving as many real-world HGVS strings as possible

# Abstract

*Application Note abstract format: Summary / Availability and Implementation / Contact /
Supplementary Information*
*Target: ~150 words. One paragraph per heading, no sub-bullets.*

---

**Summary:**

HGVS nomenclature is the international standard for describing sequence variants in
clinical reports and the databases that aggregate them. cdot is built around a single
goal: resolve as many real-world HGVS strings as possible. cdot does not itself perform
HGVS resolution; it supplies the versioned transcript-to-genome alignment data that the
established HGVS libraries already depend on, extending the existing ecosystem rather
than replacing it. The dominant such resource, UTA, requires PostgreSQL infrastructure,
lacks Ensembl support, and covers only ~{{ literature.uta_count | commas }} transcript
alignments. cdot instead provides {{ coverage.total_count | commas }} versioned
transcript/genome alignments from RefSeq and Ensembl across GRCh37, GRCh38, and
T2T-CHM13v2.0, distributed as a single gzipped JSON file or served on demand over HTTP
via a REST API (cdotlib.org), with no database required. To recover the malformed
strings that reach clinical search boxes and importers, cdot adds a parser-independent
cleaning step
(`clean_hgvs()`) that repairs common formatting errors before resolution, and an opt-in
adjacent-version fallback that substitutes a retired transcript version for the nearest
available one, with a build-independent check that the substitution does not move the
variant. Python
integrations are provided for both major HGVS libraries (biocommons/hgvs and PyHGVS).
cdot stores MANE Select and Ensembl canonical tags enabling gene-symbol HGVS lookup via
programmatic canonical transcript selection, and is the first transcript data source to
bring the T2T-CHM13v2.0 assembly to the Python HGVS libraries.

**Availability and Implementation:**

https://github.com/SACGF/cdot; `pip install cdot`; MIT licence. Data files (JSON.gz) at
cdotlib.org and [Zenodo DOI].

**Contact:** [email]

**Supplementary information:** Supplementary data are available at *Bioinformatics*
online.

---

*Fill in: Zenodo DOI, contact email. Run analysis scripts to populate coverage.csv and
benchmark.csv.*
