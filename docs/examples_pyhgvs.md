# PyHGVS example code

> **Note:** PyHGVS is no longer maintained — prefer the [biocommons HGVS examples](examples_biocommons.md).
> The PyHGVS integration exists for legacy compatibility.

Needs the `fasta` extra for `pysam`: `pip install cdot[fasta]`.

## HGVS to VCF

```python
import pyhgvs
from pysam.libcfaidx import FastaFile
from cdot.pyhgvs.pyhgvs_transcript import JSONPyHGVSTranscriptFactory, RESTPyHGVSTranscriptFactory

genome = FastaFile("/data/annotation/fasta/GCF_000001405.25_GRCh37.p13_genomic.fna.gz")
factory = RESTPyHGVSTranscriptFactory()
# factory = JSONPyHGVSTranscriptFactory(["./cdot-0.2.32.refseq.grch37.json.gz"])  # Uses local JSON file
pyhgvs.parse_hgvs_name('NM_001637.3:c.1582G>A', genome, get_transcript=factory.get_transcript_grch37)
```

## VCF to HGVS

```python
import pyhgvs as hgvs
from pysam.libcfaidx import FastaFile
from cdot.pyhgvs.pyhgvs_transcript import JSONPyHGVSTranscriptFactory

factory = JSONPyHGVSTranscriptFactory(["./cdot-0.2.32.refseq.grch37_grch38.json.gz"])

transcript = factory.get_transcript_grch37("NM_000352.3")  # change to 38 if you want that build
chrom, offset, ref, alt = ('chr11', 17496508, 'T', 'C')
genome = FastaFile("/data/annotation/fasta/GCF_000001405.25_GRCh37.p13_genomic.fna.gz")
hgvs_name = hgvs.format_hgvs_name(
    chrom, offset, ref, alt, genome, transcript)
```
