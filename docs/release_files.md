# GitHub release file details

We distribute pre-built `.json.gz` files of transcripts via
[GitHub releases](https://github.com/SACGF/cdot/releases).

The main files contain **every historical transcript** - useful for maximum compatibility when you
want to resolve as many HGVS strings as possible.

The gene symbol associated with a gene (and thus a transcript version) can change between GTFs (i.e.
it is not part of what is frozen in a transcript version). When generating these historical files, we
use the gene/transcript version from the **most recent** GTF.

| Annotation Consortium | Build | File |
|------------------------|-------|------|
| Ensembl | GRCh37 | `cdot-0.2.36.ensembl.GRCh37.json.gz` |
| Ensembl | GRCh38 | `cdot-0.2.36.ensembl.GRCh38.json.gz` |
| Ensembl | T2T-CHM13v2.0 | `cdot-0.2.36.ensembl.T2T-CHM13v2.0.json.gz` |
| Ensembl | All builds | `cdot-0.2.36.all-builds-ensembl-grch37_grch38_t2t-chm13v2.0.json.gz` |
| RefSeq | GRCh37 | `cdot-0.2.36.refseq.GRCh37.json.gz` |
| RefSeq | GRCh38 | `cdot-0.2.36.refseq.GRCh38.json.gz` |
| RefSeq | T2T-CHM13v2.0 | `cdot-0.2.36.refseq.T2T-CHM13v2.0.json.gz` |
| RefSeq | All builds | `cdot-0.2.36.all-builds-refseq-grch37_grch38_t2t-chm13v2.0.json.gz` |

The all-builds file holds GRCh37, GRCh38 and T2T-CHM13v2.0 coordinates in one file, so you only need
it or the per-build files, not both. To find these files for the latest release in code, see
[Finding the latest data programmatically](local_json_files.md#finding-the-latest-data-programmatically).

Because the files above contain the most recent transcript, if you want the transcript versions /
gene symbols to match exactly what was in a release, you need to use `.json.gz` files produced from a
**single GTF/GFF release**. We publish some of these, eg the versions used in
[VariantGrid](https://github.com/SACGF/variantgrid) to match the
[release version used in Ensembl VEP](https://asia.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache_content)
annotation:

| Annotation Consortium | Build | File |
|------------------------|-------|------|
| Ensembl | GRCh37 | `cdot-0.2.36.Homo_sapiens_GRCh37_Ensembl_87.gtf.json.gz` |
| Ensembl | GRCh38 | `cdot-0.2.36.Homo_sapiens_GRCh38_Ensembl_112.gtf.json.gz` |
| Ensembl | GRCh38 | `cdot-0.2.36.Homo_sapiens_GRCh38_Ensembl_113.gtf.json.gz` |
| Ensembl | GRCh38 | `cdot-0.2.36.Homo_sapiens_GRCh38_Ensembl_114.gtf.json.gz` |
| Ensembl | GRCh38 | `cdot-0.2.36.Homo_sapiens_GRCh38_Ensembl_115.gtf.json.gz` |
| Ensembl | GRCh38 | `cdot-0.2.36.Homo_sapiens_GRCh38_Ensembl_116.gtf.json.gz` |
| RefSeq | GRCh37 | `cdot-0.2.36.Homo_sapiens_GRCh37_RefSeq_105.20201022.gff.json.gz` |
| RefSeq | GRCh37 | `cdot-0.2.36.Homo_sapiens_GRCh37_RefSeq_105.20220307.gff.json.gz` |
| RefSeq | GRCh38 | `cdot-0.2.36.Homo_sapiens_GRCh38_RefSeq_110.gff.json.gz` |
| RefSeq | GRCh38 | `cdot-0.2.36.Homo_sapiens_GRCh38_RefSeq_RS_2025_08.gff.json.gz` |

The list of released files is set in
[`.github/workflows/data_release.yml`](../.github/workflows/data_release.yml) (step "Assemble release
assets and notes"). Older releases used different file names (eg data_v0.2.25 used lowercase builds such as
`cdot-0.2.25.refseq.grch38.json.gz`).

## See also

- [Using local downloaded JSON.gz files](local_json_files.md)
- [JSON data format reference](json_data_format.md)
