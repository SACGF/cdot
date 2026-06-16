# Using local downloaded JSON.gz files

Download the JSON.gz files from [GitHub releases](https://github.com/SACGF/cdot/releases). See
[GitHub release file details](release_files.md) for what each file contains.

You only need ONE of the combined 37/38 file or the individual genome builds, not both.

## Loading in HGVS libraries

You can pass a list of JSON.gz files to the HGVS loaders, e.g. to load both RefSeq and Ensembl:

```python
from cdot.hgvs.dataproviders import JSONDataProvider

local_json = [
    "./cdot-0.2.32.refseq.grch37_grch38.json.gz",
    "./cdot-0.2.32.ensembl.grch37_grch38.json.gz",
]
hdp = JSONDataProvider(local_json)
```

## See also

- [Biocommons HGVS examples](examples_biocommons.md)
- [Advanced usage](advanced_usage.md) - fixing messy HGVS input and bulk retrieval
