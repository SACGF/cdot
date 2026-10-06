# Checking for new annotation sources (maintainers)

Ensembl, RefSeq and UTA publish new annotation releases on their own schedules. Every source cdot
ingests is listed in [`cdot_transcripts.yaml`](../generate_transcript_data/cdot_transcripts.yaml),
so checking for new data means comparing that file against the upstream directory listings. Do this
every couple of months, and always before a [data release](data_release_workflow.md).

## Doing this with an AI agent

Point the agent at this file, eg "follow docs/checking_for_new_sources.md". It should:

1. Run the check script below. If a check errors, look at that listing by hand.
2. If everything is up to date, show the script's output and stop.
3. For a new Ensembl GRCh38 or RefSeq `RS_` release, make the changes listed under
   [Adding a new source](#adding-a-new-source), run the [parse check](#parse-check) and compare with
   the expected numbers, then run `python -m pytest tests/`. Work on a branch, not main, and ask
   before pushing or starting the data-release workflow.
4. For new historical alignments, a new UTA schema, a new Ensembl T2T geneset or an unknown RefSeq
   assembly accession, report it and ask before changing anything.
5. If anything here turned out to be wrong (a listing moved, the parse check numbers changed), fix
   this file and the script.

## Run the check

```bash
python3 generate_transcript_data/check_new_sources.py
```

It reads the yaml, fetches the upstream listings (no annotation files are downloaded, it takes a few
seconds) and prints one line per source. The exit code is 0 when everything is up to date and 1
when there is something new or a listing couldn't be read. Example output:

```
Ensembl GRCh38               yaml has 116            upstream latest 116            up to date
Ensembl GRCh37 (geneset)     yaml has 87             upstream latest 87             up to date
Ensembl T2T (rapid release)  yaml has 2022_07        upstream latest 2022_07        up to date
RefSeq GRCh37                yaml has RS_2024_09     upstream latest RS_2024_09     up to date
RefSeq GRCh38                yaml has RS_2025_08     upstream latest RS_2025_08     up to date
RefSeq T2T-CHM13v2.0         yaml has RS_2025_08     upstream latest RS_2025_08     up to date
RefSeq GRCh38 historical     yaml has RS_2024_08     upstream latest RS_2024_08     up to date
UTA                          yaml has uta_20241220   upstream latest uta_20241220   up to date
GENCODE (for HGNC IDs)       yaml has 50             upstream latest 50             up to date
```

## Where it looks

If the script errors (a site changed its layout), check these by hand:

| Source | Listing | What a new release looks like |
|---|---|---|
| Ensembl GRCh38 | <https://ftp.ensembl.org/pub/> | a new `release-N/`, with `gtf/homo_sapiens/Homo_sapiens.GRCh38.N.gtf.gz` in it |
| Ensembl GRCh37 | <https://ftp.ensembl.org/pub/grch37/> | a new `release-N/` appears with every Ensembl release, but the GTF inside is still `GRCh37.87`. GRCh37 is frozen, so a new geneset number would be a surprise |
| Ensembl T2T | <https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/> | a new `YYYY_MM/` directory |
| RefSeq (all builds) | <https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/> | a new `<assembly>-RS_YYYY_MM/`. `GCF_000001405.25` is GRCh37, `GCF_000001405.40` is GRCh38, `GCF_009914755.1` is T2T-CHM13v2.0 |
| RefSeq GRCh38 historical alignments (#51) | <https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/historical/GRCh38/> | a new `GCF_000001405.40-RS_YYYY_MM_historical/` |
| UTA | <https://hub.docker.com/r/biocommons/uta/tags> | a new `uta_YYYYMMDD` tag |
| GENCODE | <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/> | a new `release_N/` (only used to look up HGNC IDs for an Ensembl release) |

The script also flags a RefSeq `RS_` directory under an assembly accession it doesn't know (eg
`GCF_000001405.41` for a new GRCh38 patch, or a new T2T version). That needs a decision about how to
handle the new assembly, not just a yaml line, so raise an issue.

## Adding a new source

Sources are merged in yaml order and later entries win, so a new release goes at the end of its
build's block. Look at the previous release's entry and follow its naming. The commits that added
earlier releases are a good template: `git log --oneline -- generate_transcript_data/cdot_transcripts.yaml`.

### New Ensembl GRCh38 release

1. `cdot_transcripts.yaml`: add `Homo_sapiens_GRCh38_Ensembl_<N>.gtf` with its URL at the end of
   `ensembl: GRCh38`.
2. `ensembl_gencode.csv`: add the `<N>,<GENCODE release>` row at the top. Recent GENCODE releases
   are Ensembl minus 66 (116 is 50), but confirm it in the GENCODE release notes. Without this row the
   pipeline only prints a warning and that release's transcripts lose their HGNC IDs. The check
   script reports any Ensembl release in the yaml with no row.
3. Release assets: add the per-source file to the `files` list in
   [`.github/workflows/data_release.yml`](../.github/workflows/data_release.yml) (step "Assemble
   release assets and notes") and in
   [`github_release_upload.sh`](../generate_transcript_data/github_release_upload.sh). These two lists
   must stay the same. They currently ship the last five Ensembl GRCh38 releases individually, so
   drop the oldest when adding one.
4. Update the example release number in [create_data_from_scratch.md](create_data_from_scratch.md).
5. Parse check (below).
6. `CHANGELOG-data.md`: an `### Added` entry under `[unreleased]`, eg "Ensembl 117 (GRCh38)".

### New RefSeq annotation release (GRCh38, T2T or GRCh37)

1. `cdot_transcripts.yaml`: add `Homo_sapiens_<build>_RefSeq_RS_YYYY_MM.gff` with its URL at the end
   of `refseq: <build>`. The GRCh38 and T2T releases usually come out together, so check for both.
2. Release assets: for GRCh38 the latest RS release is shipped individually. In both lists from the
   Ensembl step 3, replace the old `RS_` file with the new one. GRCh37 and T2T per-source files
   aren't shipped, only the merged files.
3. Update the `RS_` example in [create_data_from_scratch.md](create_data_from_scratch.md).
4. Parse check (below).
5. `CHANGELOG-data.md`: an `### Added` entry.

### New historical alignments, UTA schema, or Ensembl T2T geneset

These are rarer and less mechanical, so talk them over before changing anything:

- **Historical alignments**: these sit low in the merge order on purpose (after UTA, before the
  official releases). Decide whether the new set replaces the old one or is added alongside it.
- **UTA**: the schema name appears in both the GRCh37 and GRCh38 `uta` entries (and in their keys).
  UTA sources are the lowest priority, so a new schema mostly adds old transcript versions. Compare
  per-source counts with the previous data release afterwards (step 4 of the
  [release workflow](data_release_workflow.md)).
- **Ensembl T2T**: add the new geneset to `ensembl: T2T-CHM13v2.0`.

## Parse check

Never run the full file locally. Instead, convert just the mitochondrial genome from the new file.
It's small, has coding transcripts on both strands, and uses a non-standard genetic code. If the
format changed (as the Ensembl GFF3 protein version did in 114, #135), this is where it shows up.

```bash
export PYTHONPATH=$(pwd)
# Gene info only adds summaries, so an empty one is fine here
echo '{"gene_info": {}, "api_retrieval_date": null}' | gzip > /tmp/empty_gene_info.json.gz

# Ensembl GTF: header lines plus contig MT
curl -s <GTF URL> | zcat | awk '/^#/ || $1 == "MT"' | gzip > /tmp/new_ensembl_mt.gtf.gz
generate_transcript_data/cdot_json.py gtf_to_json /tmp/new_ensembl_mt.gtf.gz --annotation-consortium ensembl \
    --url <GTF URL> --genome-build GRCh38 --gene-info-json /tmp/empty_gene_info.json.gz \
    --output /tmp/new_ensembl_mt.json.gz

# RefSeq GFF3: header lines plus NC_012920.1
curl -s <GFF URL> | zcat | awk '/^#/ || $1 == "NC_012920.1"' | gzip > /tmp/new_refseq_mt.gff.gz
generate_transcript_data/cdot_json.py gff3_to_json /tmp/new_refseq_mt.gff.gz --annotation-consortium refseq \
    --url <GFF URL> --genome-build GRCh38 --gene-info-json /tmp/empty_gene_info.json.gz \
    --output /tmp/new_refseq_mt.json.gz
```

Pass `--annotation-consortium` because the MT rows alone are too few for it to be detected, and the
generic conventions fail on RefSeq. Each download streams the whole compressed file but only keeps a
few hundred lines, and takes a few seconds.

What to expect (Ensembl 116 and RefSeq RS_2025_08 gave these): Ensembl has 37 transcripts, 13 with a
protein (eg `ENST00000361390.2` with `ENSP00000354687.2`), and RefSeq has 13 transcripts, all with a
protein (`fake-rna-ATP6` etc, the RefSeq mito fakes). Every coding transcript should have
`translation` `{"transl_table": 2}`. Fewer transcripts, missing proteins or a traceback means the
format changed. Then run the tests (`python -m pytest tests/`). The real build happens in the
[data-release workflow](data_release_workflow.md).
