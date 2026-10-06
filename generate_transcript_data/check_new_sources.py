#!/usr/bin/env python3
"""
Check the upstream FTP/registry listings for annotation releases newer than the ones in
cdot_transcripts.yaml. Read only: it never downloads an annotation file.

See docs/checking_for_new_sources.md for what to do with the output.
"""

import csv
import json
import os
import re
import sys
import urllib.request
from argparse import ArgumentParser

import yaml

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ENSEMBL_FTP = "https://ftp.ensembl.org/pub"
NCBI_ANNOTATION_RELEASES = "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/"
NCBI_HISTORICAL_GRCH38 = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/historical/GRCh38/"
ENSEMBL_T2T_GENESETS = f"{ENSEMBL_FTP}/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/"
GENCODE_FTP = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
UTA_DOCKER_TAGS = "https://hub.docker.com/v2/repositories/biocommons/uta/tags?page_size=100"

# RefSeq annotation release directories are named <assembly accession>-RS_<yyyy>_<mm>
REFSEQ_ASSEMBLIES = {
    "GCF_000001405.25": "GRCh37",
    "GCF_000001405.40": "GRCh38",
    "GCF_009914755.1": "T2T-CHM13v2.0",
}


def fetch(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "cdot check_new_sources"})
    with urllib.request.urlopen(req, timeout=120) as response:
        return response.read().decode()


def list_dirs(url: str) -> list[str]:
    """ Sub-directory names from an Apache/nginx style FTP-over-HTTP listing """
    return sorted(set(re.findall(r'href="([^"/?][^"/]*)/"', fetch(url))))


def load_config():
    with open(os.path.join(BASE_DIR, "cdot_transcripts.yaml")) as f:
        return yaml.safe_load(f)["config"]


def source_urls(config, consortium, build) -> list[str]:
    urls = []
    for source in (config[consortium].get(build) or {}).values():
        if source:
            urls.extend(source.get("urls", [source["url"]] if "url" in source else []))
    return urls


def report(label, have, latest, newer, notes=None):
    status = "NEW: " + ", ".join(newer) if newer else "up to date"
    print(f"{label:<28} yaml has {have:<14} upstream latest {latest:<14} {status}")
    for note in notes or []:
        print(f"{'':<28} ! {note}")
    return bool(newer)


def check_ensembl_grch38(config):
    have = [int(m.group(1)) for u in source_urls(config, "ensembl", "GRCh38")
            if (m := re.search(r"/release-(\d+)/", u))]
    releases = sorted(int(d.split("-")[1]) for d in list_dirs(ENSEMBL_FTP + "/") if re.fullmatch(r"release-\d+", d))
    newer = []
    notes = []
    for r in releases:
        if r <= max(have):
            continue
        gtf = f"Homo_sapiens.GRCh38.{r}.gtf.gz"
        if gtf in fetch(f"{ENSEMBL_FTP}/release-{r}/gtf/homo_sapiens/"):
            newer.append(str(r))
        else:
            notes.append(f"release-{r} directory exists but has no {gtf} yet")
    return report("Ensembl GRCh38", str(max(have)), str(releases[-1]), newer, notes)


def check_ensembl_grch37(config):
    """ GRCh37 is frozen at geneset 87 and re-published in every grch37/release-N, so look at the
        geneset version in the GTF filename, not the directory number """
    have = [int(m.group(1)) for u in source_urls(config, "ensembl", "GRCh37")
            if (m := re.search(r"GRCh37\.(\d+)\.gtf", u))]
    latest_dir = max(int(d.split("-")[1]) for d in list_dirs(f"{ENSEMBL_FTP}/grch37/") if re.fullmatch(r"release-\d+", d))
    listing = fetch(f"{ENSEMBL_FTP}/grch37/release-{latest_dir}/gtf/homo_sapiens/")
    genesets = sorted({int(g) for g in re.findall(r"Homo_sapiens\.GRCh37\.(\d+)\.gtf\.gz", listing)})
    newer = [str(g) for g in genesets if g > max(have)]
    return report("Ensembl GRCh37 (geneset)", str(max(have)), str(genesets[-1]), newer)


def check_ensembl_t2t(config):
    have = [m.group(1) for u in source_urls(config, "ensembl", "T2T-CHM13v2.0")
            if (m := re.search(r"/geneset/(\d{4}_\d{2})/", u))]
    genesets = [d for d in list_dirs(ENSEMBL_T2T_GENESETS) if re.fullmatch(r"\d{4}_\d{2}", d)]
    newer = [g for g in genesets if g > max(have)]
    return report("Ensembl T2T (rapid release)", max(have), genesets[-1], newer)


def check_refseq(config):
    found_new = False
    dirs = list_dirs(NCBI_ANNOTATION_RELEASES)
    rs_dirs = [d for d in dirs if "-RS_" in d]
    for accession, build in REFSEQ_ASSEMBLIES.items():
        have = [m.group(1) for u in source_urls(config, "refseq", build)
                if (m := re.search(rf"/{re.escape(accession)}-(RS_\d{{4}}_\d{{2}})/", u))]
        upstream = sorted(d.split("-", 1)[1] for d in rs_dirs if d.startswith(accession + "-"))
        have_latest = max(have) if have else "-"
        newer = [r for r in upstream if not have or r > max(have)]
        found_new |= report(f"RefSeq {build}", have_latest, upstream[-1] if upstream else "-", newer)

    # A new assembly patch (eg GCF_000001405.41) or T2T version gets a new accession prefix
    unknown = sorted({d.split("-", 1)[0] for d in rs_dirs} - set(REFSEQ_ASSEMBLIES))
    if unknown:
        found_new = True
        print(f"{'RefSeq':<28} ! unknown assembly accession(s) with RS_ releases: {', '.join(unknown)}")
    return found_new


def check_refseq_historical(config):
    """ NCBI historical transcript alignments (#51), published separately from annotation releases """
    have = [m.group(1) for u in source_urls(config, "refseq", "GRCh38")
            if (m := re.search(r"/GCF_000001405\.40-(RS_\d{4}_\d{2})_historical/", u))]
    upstream = sorted(m.group(1) for d in list_dirs(NCBI_HISTORICAL_GRCH38)
                      if (m := re.fullmatch(r"GCF_000001405\.40-(RS_\d{4}_\d{2})_historical", d)))
    newer = [r for r in upstream if r > max(have)]
    return report("RefSeq GRCh38 historical", max(have), upstream[-1], newer)


def check_uta(config):
    have = sorted({s["uta"]["schema"] for c in config.values() for b in c.values() if b
                   for s in b.values() if s and "uta" in s})
    tags = [t["name"] for t in json.loads(fetch(UTA_DOCKER_TAGS))["results"]]
    schemas = sorted(t for t in tags if re.fullmatch(r"uta_\d{8}[a-z]?", t))
    newer = [s for s in schemas if s[:12] > max(have)[:12]]
    return report("UTA", max(have), schemas[-1], newer)


def check_gencode_mapping(config):
    """ Every Ensembl release used needs a row in ensembl_gencode.csv, or its HGNC IDs are silently
        skipped (the Snakefile only prints a warning). Also report the latest GENCODE release. """
    with open(os.path.join(BASE_DIR, "ensembl_gencode.csv")) as f:
        mapped = {r["Ensembl release"].strip(): r["GENCODE release"].strip() for r in csv.DictReader(f)}
    used = sorted({m.group(1) for b in ("GRCh37", "GRCh38") for u in source_urls(config, "ensembl", b)
                   if (m := re.search(r"/release-(\d+)/", u))}, key=int)
    missing = [r for r in used if r not in mapped]
    latest = max(int(d.split("_")[1]) for d in list_dirs(GENCODE_FTP) if re.fullmatch(r"release_\d+", d))
    notes = [f"no ensembl_gencode.csv row for Ensembl {', '.join(missing)}"] if missing else []
    have = max(int(v) for v in mapped.values())
    return report("GENCODE (for HGNC IDs)", str(have), str(latest), [str(latest)] if latest > have else [], notes)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.parse_args()
    config = load_config()
    checks = [check_ensembl_grch38, check_ensembl_grch37, check_ensembl_t2t, check_refseq, check_refseq_historical,
              check_uta,
              check_gencode_mapping]
    found_new = False
    for check in checks:
        try:
            found_new |= check(config)
        except Exception as e:
            found_new = True
            print(f"{check.__name__:<28} ERROR {e!r} (check the listing by hand)")
    sys.exit(1 if found_new else 0)


if __name__ == "__main__":
    main()
