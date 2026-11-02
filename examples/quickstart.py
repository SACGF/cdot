"""
cdot quickstart: find the latest data release, download one file, load it into biocommons HGVS,
then clean and resolve a few real-world HGVS strings to genomic coordinates.

    pip install cdot
    python examples/quickstart.py                             # SeqRepo or NCBI for sequences
    python examples/quickstart.py --fasta GRCh38_genomic.fna.gz   # local genome FASTA (needs cdot[fasta])

The FASTA must use NCBI contig names (NC_000001.11 etc, eg the GCF_000001405.39 GRCh38 assembly),
bgzipped and indexed. See docs/fasta_seqfetcher.md.
"""
import argparse
import os

import hgvs.parser
import requests
from hgvs.assemblymapper import AssemblyMapper
from hgvs.dataproviders.seqfetcher import SeqFetcher

from cdot.data_release import get_latest_combo_file_urls
from cdot.hgvs import fix_hgvs
from cdot.hgvs.dataproviders import ChainedSeqFetcher, JSONDataProvider

GENOME_BUILD = "GRCh38"
EXAMPLES = [
    "BRCA2(NM_000059.4):c.68_69DEL p.(Glu23fs)",  # gene/transcript swapped, uppercase DEL, trailing p.
    "NM_001754.5 c.557T>A",  # RUNX1, missing colon
    "BRCA2:c.316+5G>A",  # gene symbol only, resolved to the MANE Select transcript
]


def download_latest_data(data_dir, annotation_consortium="refseq", genome_build=GENOME_BUILD):
    """ Downloads the per-build file from the latest compatible data release (skipped if already present) """
    url = get_latest_combo_file_urls([annotation_consortium], [genome_build])[0]
    path = os.path.join(data_dir, url.rsplit("/", 1)[-1])
    if not os.path.exists(path):
        print(f"Downloading {url}")
        with requests.get(url, stream=True, timeout=30) as response, open(path, "wb") as f:
            response.raise_for_status()
            for chunk in response.iter_content(chunk_size=1 << 20):
                f.write(chunk)
    return path


def build_provider(json_path, fasta=None):
    """ FastaSeqFetcher splices transcript sequences out of the genome, so transcripts that aren't in
        SeqRepo still resolve. It is chained after SeqFetcher so SeqRepo sequences win where they exist """
    seqfetcher = None
    if fasta:
        from cdot.hgvs.dataproviders import FastaSeqFetcher
        seqfetcher = ChainedSeqFetcher(SeqFetcher(), FastaSeqFetcher(fasta))
    return JSONDataProvider([json_path], seqfetcher=seqfetcher)


def resolve(hdp, genome_build, hgvs_strings):
    """ Yields (original, fixed, fixes, var_g) for each HGVS string """
    am = AssemblyMapper(hdp, assembly_name=genome_build, alt_aln_method="splign", replace_reference=True)
    hp = hgvs.parser.Parser()
    for hgvs_string in hgvs_strings:
        fixed, fixes = fix_hgvs(hgvs_string, hdp, genome_build)
        var_g = am.c_to_g(hp.parse_hgvs_variant(fixed))
        yield hgvs_string, fixed, fixes, var_g


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", default=".", help="where to keep the downloaded JSON.gz (default: cwd)")
    parser.add_argument("--fasta", help="bgzipped, indexed genome FASTA with NCBI contig names")
    args = parser.parse_args()

    json_path = download_latest_data(args.data_dir)
    hdp = build_provider(json_path, fasta=args.fasta)
    for hgvs_string, fixed, fixes, var_g in resolve(hdp, GENOME_BUILD, EXAMPLES):
        print(f"{hgvs_string!r} -> {fixed} -> {var_g}")
        for fix in fixes:
            print(f"    {fix}")


if __name__ == "__main__":
    main()
