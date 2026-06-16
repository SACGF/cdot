# GitHub release file details

We distribute pre-built `.json.gz` files of transcripts via
[GitHub releases](https://github.com/SACGF/cdot/releases).

The main files contain **every historical transcript** — useful for maximum compatibility when you
want to resolve as many HGVS strings as possible.

The gene symbol associated with a gene (and thus a transcript version) can change between GTFs (i.e.
it is not part of what is frozen in a transcript version). When generating these historical files, we
use the gene/transcript version from the **most recent** GTF.

| Annotation Consortium | Build | Example File |
|------------------------|-------|--------------|
| Ensembl | GRCh37 | `cdot-0.2.14.ensembl.grch37.json.gz` |
| Ensembl | GRCh38 | `cdot-0.2.14.ensembl.grch38.json.gz` |
| Ensembl | GRCh37/GRCh38 | `cdot-0.2.14.ensembl.grch37_grch38.json.gz` |
| RefSeq | GRCh37 | `cdot-0.2.14.refseq.grch37.json.gz` |
| RefSeq | GRCh38 | `cdot-0.2.14.refseq.grch38.json.gz` |
| RefSeq | GRCh37/GRCh38 | `cdot-0.2.14.refseq.grch37_grch38.json.gz` |

Because the files above contain the most recent transcript, if you want the transcript versions /
gene symbols to match exactly what was in a release, you need to use `.json.gz` files produced from a
**single GTF release**.

We host some of the versions used in [VariantGrid](https://github.com/SACGF/variantgrid), where they
match the [release version used in Ensembl VEP](https://asia.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache_content)
annotation, e.g.:

| Annotation Consortium | Build | Example File | VEP release |
|------------------------|-------|--------------|-------------|
| RefSeq | GRCh37 | `cdot-0.2.14.GCF_000001405.25_GRCh37.p13_genomic.105.20201022.gff.json.gz` | 100-109 |
| RefSeq | GRCh38 | `cdot-0.2.14.GCF_000001405.39_GRCh38.p13_genomic.109.20211119.gff.json.gz` | -108 |

## See also

- [Using local downloaded JSON.gz files](local_json_files.md)
- [JSON data format reference](json_data_format.md)
