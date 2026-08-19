# Building a data release (maintainers)

How official cdot transcript data releases are built and published. This is maintainer
documentation: users get data from the [GitHub releases](https://github.com/SACGF/cdot/releases)
or the REST API and never need any of this. To build data for yourself locally, see
[create data from scratch](create_data_from_scratch.md).

Data releases are built by the [`data-release` GitHub Actions workflow](../.github/workflows/data_release.yml),
which runs the [Snakemake pipeline](../generate_transcript_data/Snakefile) on hosted runners and
uploads the results to a **draft** release for manual validation. Nothing is public until the draft
is published: draft releases are invisible to the releases API, so `cdot.data_release` clients
cannot see the files until you press Publish.

## Release cycle

1. Update [`generate_transcript_data/cdot_transcripts.yaml`](../generate_transcript_data/cdot_transcripts.yaml)
   with any new sources (new Ensembl release, new RefSeq annotation release). Sources are merged in
   YAML order, later entries override earlier ones for the same transcript version, so position
   matters (see the comments in the file).
2. Bump `JSON_SCHEMA_VERSION` in
   [`generate_transcript_data/json_schema_version.py`](../generate_transcript_data/json_schema_version.py)
   and move the `[unreleased]` entries in [`CHANGELOG-data.md`](../CHANGELOG-data.md) under a heading
   for the new version. The version is embedded in every output filename and becomes the release tag
   (`data_v<version>`).
3. Push to main, then start the workflow: GitHub → Actions → "data-release" → "Run workflow"
   (or `gh workflow run data-release`). A full build parses every source and takes a few hours;
   the matrix runs one job per source file.
4. When the run finishes, read its summary page. It shows:
   - whether the merge printed the "different coding coordinates" report (transcripts whose codon
     coordinates differ between builds; the stale build coordinates are dropped, see
     [#54](https://github.com/SACGF/cdot/issues/54) for why the remaining cases are expected), and
   - per-file `url_counts`, how many transcripts each source contributed to each released file.
     Compare against the previous release's notes: a source whose count collapses indicates an
     import problem (this is how the missing-UTA regression in 0.2.33 would have been caught,
     [#95](https://github.com/SACGF/cdot/issues/95)).
5. Spot-check the draft release files, then Publish. Publishing creates the `data_v<version>` tag
   and makes the release visible to `cdot.data_release` clients. There is no tag-first flow: the
   tag is a product of publishing, not a trigger.
6. The cdotlib.org REST server ([cdot_rest](https://github.com/SACGF/cdot_rest)) is deployed
   separately and needs its data updated through its own process.

## Retrying failures

Individual jobs fail occasionally (source download hiccups, UTA server outages). Two mechanisms:

- **Same commit**: while the fix needs no code change, `gh run rerun <run-id> --failed` reruns just
  the failed jobs and everything downstream, reusing successful artifacts. Note reruns always use
  the run's original commit.
- **New commit**: if the fix required a code change, start a new run and pass the previous run's ID
  as the `reuse_run_id` input. Each parse job first tries to download its output from that run's
  artifacts; if present, Snakemake sees the output already exists and does nothing, so only the
  sources that actually failed (plus merge and the draft) are rebuilt. Artifacts are kept for
  7 days. This only works at the same data version, a version bump changes every filename and
  forces a full reparse.

## Pipeline notes

- One matrix job per source in `cdot_transcripts.yaml` (~77 jobs), capped at 10 in parallel so
  NCBI/Ensembl do not drop us for opening too many connections.
- The NCBI gene info file (Entrez summaries) is built once and shared via an artifact. Its filename
  embeds a date; the workflow pins `CDOT_GENE_INFO_DATE` so jobs straddling midnight agree on it.
- The UTA export query aggregates server-side for ~17 minutes before returning any rows. The psql
  connection sends TCP keepalives so NAT (on hosted runners and elsewhere) does not drop it during
  the silent phase, runs `ON_ERROR_STOP` so a partial export fails the rule instead of leaving a
  truncated CSV, and retries twice for transient outages.
- The merge/combine steps peak around 4-6 GB RAM, well within the 16 GB on hosted runners.
- Assets attached to the release are the two all-builds files, the six per-build merged files, and
  a small set of individual per-source files, the same list as
  [`github_release_upload.sh`](../generate_transcript_data/github_release_upload.sh) (the pre-CI
  manual upload script). The merge log and generated release notes are kept in the run's
  `release-files` artifact.
