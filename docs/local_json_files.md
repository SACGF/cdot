# Using local downloaded JSON.gz files

Download the JSON.gz files from [GitHub releases](https://github.com/SACGF/cdot/releases). See
[GitHub release file details](release_files.md) for what each file contains.

You only need ONE of the all-builds file (GRCh37, GRCh38 and T2T-CHM13v2.0 in one file) or the
individual genome build files, not both.

## Finding the latest data programmatically

`cdot.data_release` finds the newest data release compatible with your installed cdot (via the
GitHub releases API) and returns download URLs:

```python
from cdot.data_release import get_latest_release_version, get_latest_combo_file_urls

get_latest_release_version()  # eg '0.2.36'
get_latest_combo_file_urls(["refseq"], ["GRCh37", "GRCh38"])  # per-build files
get_latest_combo_file_urls(["refseq", "ensembl"], ["all-builds"])  # all-builds files
```

`genome_builds` values are matched case insensitively against the build in the filename. `"all-builds"`
selects the all-builds file, as does an underscore-joined combination such as `"grch37_grch38"` (the
name of the combined file in old releases).

## Loading in HGVS libraries

You can pass a list of JSON.gz files to the HGVS loaders, e.g. to load both RefSeq and Ensembl:

```python
from cdot.hgvs.dataproviders import JSONDataProvider

local_json = [
    "./cdot-0.2.36.all-builds-refseq-grch37_grch38_t2t-chm13v2.0.json.gz",
    "./cdot-0.2.36.all-builds-ensembl-grch37_grch38_t2t-chm13v2.0.json.gz",
]
hdp = JSONDataProvider(local_json)
```

## See also

- [Biocommons HGVS examples](examples_biocommons.md)
- [Advanced usage](advanced_usage.md) - fixing messy HGVS input and bulk retrieval
