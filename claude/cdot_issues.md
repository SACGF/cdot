# cdot — issues, quirks, and improvement ideas

Notes from using [cdot](https://github.com/SACGF/cdot) (SACGF/cdot) extensively
for transcript annotation analysis. Collected to guide future contributions back
to the project.

---

## Tag name inconsistency between Ensembl and RefSeq

The same annotation concept uses different string representations depending on
the source consortium:

| Tag | Ensembl | RefSeq |
|---|---|---|
| MANE Select | `MANE_Select` | `MANE Select` (space) |
| GENCODE Basic | `gencode_basic` or `basic` | not applicable |
| GENCODE Primary | `gencode_primary` or `GENCODE Primary` | not applicable |

The MANE difference is the most problematic — downstream code must check both
spellings or normalise on load. See `cdot_utils.py:get_tags()` and
`cdot_data_fields.md` lines 88–95 for the workaround.

**Suggestion:** cdot could normalise tag names during GFF/GTF parsing so
consumers don't need to handle both forms. At minimum, `MANE Select` →
`MANE_Select` on the RefSeq side would unify the most clinically important tag.

Not yet filed as an issue.

---

## Biotype list length differs between Ensembl and RefSeq

Ensembl biotype is a 2-element list: `["mRNA", "protein_coding"]` (SO term +
internal biotype). RefSeq biotype is a 1-element list: `["mRNA"]`.

This means `bt[1]` is an `IndexError` on RefSeq data. Our accessor
(`cdot_utils.py:get_biotype()`) uses only `bt[0]`, which works for both, but the
inconsistency is a footgun for anyone accessing the raw JSON.

**Suggestion:** Either pad RefSeq to 2 elements (e.g. `["mRNA", "mRNA"]`) or
document the difference prominently. Related to
[#77](https://github.com/SACGF/cdot/issues/77) (document JSON format) and
[#78](https://github.com/SACGF/cdot/issues/78) (schema changes).

---

## Exon list order not guaranteed

The exons array in cdot JSON is not in any guaranteed order. Code must sort by
`exon[3]` (cdna_start) before iterating for transcript-order operations (NMD
prediction, cDNA extraction, etc.).

See `cdot_data_fields.md` lines 124–135 and `analyse_brca1.py:extract_cdna()`
for the sorting pattern.

**Suggestion:** cdot could sort exons by cdna_start during serialisation. This
would eliminate a class of silent bugs in downstream code. Related to
[#102](https://github.com/SACGF/cdot/issues/102) (document how exons work).

---

## Missing data that could be extracted from GTF/GFF

### APPRIS tier

Ensembl GTFs contain APPRIS annotations (`tag "appris_principal_1"`, etc.) that
indicate the principal protein isoform. These are not currently stored in cdot.
APPRIS is used by VEP for canonical transcript selection and would be valuable
for filtering analyses.

### NMD biotype

Ensembl classifies some transcripts as `nonsense_mediated_decay` in the biotype
field. This is captured in cdot via the biotype list, but a dedicated boolean
field would make NMD filtering more ergonomic.

### CDS phase / reading frame offset

[#76](https://github.com/SACGF/cdot/issues/76) tracks this — some transcripts
have non-zero CDS phase (e.g. ribosomal slippage). Currently not stored.

---

## BestRefSeq: all 368 BRCA1 transcripts are BestRefSeq

In RefSeq RS_2025_08, all 368 BRCA1 NM_ accessions have `source: ["BestRefSeq"]`.
This means the BestRefSeq designation alone cannot distinguish the 6 original
well-characterised transcripts from the 362 long-read-derived isoforms added in
RS_2023_03. CCDS cross-referencing (5 of 368 have CCDS) is the only way to
identify the multi-consortium-validated core set.

This is a RefSeq annotation decision rather than a cdot bug, but it's worth
noting because BestRefSeq is often treated as a quality filter and it does not
filter in this case.

---

## Documentation gaps

Several existing issues track documentation needs:

- [#77](https://github.com/SACGF/cdot/issues/77) — Document JSON.gz files and
  transcript format
- [#102](https://github.com/SACGF/cdot/issues/102) — Document how exons work
- [#78](https://github.com/SACGF/cdot/issues/78) — JSON breaking schema changes

From our experience, the most useful documentation additions would be:
1. Field type guarantees (e.g. biotype is always a list, source is always a list)
2. Which fields are Ensembl-only vs RefSeq-only vs shared
3. Which fields were added in which cdot version (breaking for anyone pinned to
   older data files)
4. Exon tuple layout with a worked example

---

## Summary of potential contributions

| Priority | Item | GitHub issue |
|---|---|---|
| High | Upgrade this project to 0.2.33 | — (internal) |
| Medium | Normalise MANE tag names across Ensembl/RefSeq | new issue |
| Medium | Guarantee exon ordering in JSON | extend #102 |
| Medium | Store APPRIS tier from Ensembl GTF | new issue |
| Low | Pad RefSeq biotype to 2 elements or document | extend #77/#78 |
| Low | Document field availability by cdot version | extend #77 |


Here's some code I used in another project:


```
def get_ensembl_cdot_path(cdot_data_dir: str, cdot_version: str, ensembl_version: str) -> str:
    """Return the path to the CDOT JSON for a given Ensembl release."""
    return os.path.join(cdot_data_dir, "ensembl", "GRCh38",
                        f"cdot-{cdot_version}.Homo_sapiens_GRCh38_Ensembl_{ensembl_version}.gtf.json.gz")


def get_refseq_cdot_path(cdot_data_dir: str, cdot_version: str, refseq_version: str) -> str:
    """Return the path to the CDOT JSON for a given RefSeq release."""
    return os.path.join(cdot_data_dir, "refseq", "GRCh38",
                        f"cdot-{cdot_version}.Homo_sapiens_GRCh38_RefSeq_{refseq_version}.json.gz")


def load_cdot_json(path: str) -> dict:
    """Load a gzipped CDOT JSON file and return the parsed dict."""
    log.info("Loading %s", path)
    return json.load(gzip.open(path))


def get_biotype(td: dict) -> str:
    """Return the primary biotype for a CDOT transcript dict.

    CDOT stores biotype as a list; the first element is used.
    Falls back to 'unknown' if absent or empty.
    """
    bt = td.get("biotype", [])
    if isinstance(bt, list) and bt:
        return bt[0]
    return bt or "unknown"


def get_tsl(td: dict) -> str:
    """Return the GRCh38 transcript support level as a single-character string.

    CDOT stores TSL under genome_builds.GRCh38.transcript_support_level as a
    string like '1', '2 (assigned to previous version 9)', 'NA', etc.
    Returns the first character when it is a digit, 'not assigned' otherwise.
    """
    gb = td.get("genome_builds", {}).get("GRCh38", {})
    tsl_raw = gb.get("transcript_support_level")
    if tsl_raw is None:
        return "not assigned"
    tsl_str = str(tsl_raw)
    if not tsl_str:
        return "not assigned"
    first_char = tsl_str[0]
    return first_char if first_char.isdigit() else "not assigned"


def get_tags(td: dict) -> list:
    """Return the list of GRCh38 GENCODE tags for a CDOT transcript dict.

    CDOT stores tags as a comma-separated string in genome_builds.GRCh38.tag.
    Returns an empty list if absent.
    """
    gb = td.get("genome_builds", {}).get("GRCh38", {})
    tag_str = gb.get("tag", "")
    return [t.strip() for t in tag_str.split(",") if t.strip()] if tag_str else []

```
