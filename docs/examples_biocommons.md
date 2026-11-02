# Biocommons HGVS example code

If you get `HGVSDataNotAvailableError` (sequences not available), look into using
[FastaSeqFetcher](fasta_seqfetcher.md) to provide transcripts. To fix messy HGVS input or speed up
bulk processing, see [Advanced usage](advanced_usage.md).

## Convert transcript (c.) to genomic (g.)

```python
import hgvs
from hgvs.assemblymapper import AssemblyMapper
from cdot.hgvs.dataproviders import JSONDataProvider, RESTDataProvider

hdp = RESTDataProvider()  # Uses API server at cdotlib.org
# hdp = JSONDataProvider(["./cdot-0.2.36.refseq.GRCh37.json.gz"])  # Uses local JSON file

am = AssemblyMapper(hdp,
                    assembly_name='GRCh37',
                    alt_aln_method='splign', replace_reference=True)

hp = hgvs.parser.Parser()
var_c = hp.parse_hgvs_variant('NM_001637.3:c.1582G>A')
am.c_to_g(var_c)
```

## Convert genomic (g.) to transcript (c.)

```python
from cdot.hgvs.dataproviders import JSONDataProvider, RESTDataProvider
import hgvs.parser
import hgvs.assemblymapper

hp = hgvs.parser.Parser()
hdp = RESTDataProvider()
am = hgvs.assemblymapper.AssemblyMapper(hdp,
                                        assembly_name='GRCh37', alt_aln_method='splign',
                                        replace_reference=True)

var_g = hp.parse_hgvs_variant('NC_000007.13:g.36561662C>T')
am.g_to_c(var_g, 'NM_001637.3')
```

## T2T-CHM13v2.0 example

Resolving against the [T2T-CHM13v2.0](https://github.com/SACGF/cdot/issues/41) assembly needs a local
FASTA (see [FastaSeqFetcher](fasta_seqfetcher.md)).

### Setup / get data

```bash
# Get FASTA
wget --quiet -O - https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz | gzip -d | bgzip > GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
samtools faidx GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz

# Get cdot data
wget https://github.com/SACGF/cdot/releases/download/data_v0.2.36/cdot-0.2.36.refseq.T2T-CHM13v2.0.json.gz
```

### Python code

```python
import hgvs
from hgvs.assemblymapper import AssemblyMapper
from cdot.hgvs.dataproviders import JSONDataProvider
from cdot.hgvs.dataproviders.fasta_seqfetcher import FastaSeqFetcher

seqfetcher = FastaSeqFetcher("./GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz")
hdp = JSONDataProvider(["./cdot-0.2.36.refseq.T2T-CHM13v2.0.json.gz"],
                       seqfetcher=seqfetcher)

am = AssemblyMapper(hdp,
                    assembly_name='CHM13v2.0',
                    alt_aln_method='splign', replace_reference=True)

hp = hgvs.parser.Parser()
var_c = hp.parse_hgvs_variant('NM_001637.4:c.1582G>A')
print(am.c_to_g(var_c))
```

Output:

```
NC_060931.1:g.36662409C>T
```

## Ensembl Tark (no cdot data)

`EnsemblTarkDataProvider` implements the same biocommons interface from Ensembl's
[Tark](https://tark.ensembl.org/) transcript archive REST API, so nothing has to be downloaded.
Each transcript is fetched from Tark the first time it is used. Tark holds Ensembl and RefSeq
transcripts, but no RefSeq-to-genome alignments, so RefSeq transcripts that align with gaps are
not supported. Sequences come from `EnsemblTarkSeqFetcher` by default (Tark for transcripts,
biocommons `SeqFetcher` for genomic contigs); `get_ensembl_tark_fasta_seqfetchers` in
`cdot.hgvs.dataproviders.fasta_seqfetcher` builds fetchers that use a local genome FASTA instead.

```python
import hgvs
from hgvs.assemblymapper import AssemblyMapper
from cdot.hgvs.dataproviders import EnsemblTarkDataProvider

hdp = EnsemblTarkDataProvider()  # GRCh37 and GRCh38

am = AssemblyMapper(hdp,
                    assembly_name='GRCh38',
                    alt_aln_method='splign', replace_reference=True)

hp = hgvs.parser.Parser()
var_c = hp.parse_hgvs_variant('ENST00000617537.5:c.1582G>A')
print(am.c_to_g(var_c))
```

Output:

```
NC_000007.14:g.36522056C>T
```
