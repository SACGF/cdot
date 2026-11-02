"""
Find the latest cdot transcript data release on GitHub.

Data is released separately from the client code, as GitHub releases tagged ``data_v<version>``
(eg ``data_v0.2.36``). These helpers query the GitHub releases API for the newest data release whose
schema version matches this client (see ``cdot.get_data_schema_int``), and return download URLs for
its assets.

Release asset names (current scheme)::

    cdot-0.2.36.refseq.GRCh38.json.gz                                    # one genome build
    cdot-0.2.36.all-builds-refseq-grch37_grch38_t2t-chm13v2.0.json.gz    # all genome builds in one file

Older releases (data_v0.2.25 and earlier) used lowercase builds (``cdot-0.2.25.refseq.grch38.json.gz``),
and early ones (eg v0.2.14) had a combined file named ``cdot-<version>.refseq.grch37_grch38.json.gz``.
Matching is case insensitive, and both schemes are supported.
"""
import re
from typing import Iterable, Iterator, Optional

import requests
import cdot

from cdot import get_data_schema_int

ALL_BUILDS = "all-builds"

_PER_BUILD_FILE_PATTERN = re.compile(r"cdot-(\d+\.\d+\.\d+)\.(refseq|ensembl)\.(.+)\.json\.gz", re.IGNORECASE)
_ALL_BUILDS_FILE_PATTERN = re.compile(r"cdot-(\d+\.\d+\.\d+)\.all-builds-(refseq|ensembl)-(.+)\.json\.gz",
                                      re.IGNORECASE)


def get_latest_data_release_tag_name():
    """ Tag name (eg 'data_v0.2.36') of the latest data release compatible with this client, or None """
    latest_data_release = get_latest_data_release()
    return latest_data_release.get('tag_name')

def _get_version_from_tag_name(tag_name, data_version=False):
    """ Returns None if doesn't match required prefix """
    release_prefix = "v"
    if data_version:
        release_prefix = "data_" + release_prefix

    if not tag_name.startswith(release_prefix):
        return None
    return tag_name.removeprefix(release_prefix)


def get_latest_data_version_and_release() -> tuple[Optional[str], dict]:
    """ Returns (data_version, release) for the newest data release with the same data schema as this
        client, where release is the GitHub API release dict. Returns (None, {}) if there is none """
    client_data_schema = get_data_schema_int(cdot.__version__)

    url = "https://api.github.com/repos/SACGF/cdot/releases?per_page=100"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    json_data = response.json()
    for release in json_data:
        tag_name = release['tag_name']  # Should look like 'v0.2.25' for code or 'data_v0.2.25' for data
        # We require a data version
        data_version = _get_version_from_tag_name(tag_name, data_version=True)
        if data_version is None:
            continue

        data_schema = get_data_schema_int(data_version)
        if data_schema != client_data_schema:
            continue
        return data_version, release
    return None, {}


def get_latest_data_release() -> dict:
    """ GitHub API release dict of the latest compatible data release ({} if none) """
    _data_version, release = get_latest_data_version_and_release()
    return release


def get_latest_release_version() -> Optional[str]:
    """ Data version (eg '0.2.36') of the latest compatible data release, or None """
    if latest_data_release := get_latest_data_release():
        if tag_name := latest_data_release.get("tag_name"):
            return _get_version_from_tag_name(tag_name, data_version=True)
    return None


def get_latest_browser_urls() -> Iterator[str]:
    """ Download URLs of every asset in the latest compatible data release """
    if latest_data_release := get_latest_data_release():
        for asset in latest_data_release["assets"]:
            yield asset["browser_download_url"]


def _filename_matches(filename: str, annotation_consortia: set[str], genome_builds: set[str]) -> bool:
    """ annotation_consortia and genome_builds must already be lower case """
    if m := _ALL_BUILDS_FILE_PATTERN.fullmatch(filename):
        _version, annotation_consortium, file_builds = m.groups()
        if annotation_consortium.lower() not in annotation_consortia:
            return False
        file_builds = file_builds.lower()
        file_build_set = set(file_builds.split("_"))
        for requested in genome_builds:
            if requested in (ALL_BUILDS, file_builds):
                return True
            # A combination of builds (eg 'grch37_grch38', the old combined file) is served by an
            # all-builds file that contains them all
            requested_set = set(requested.split("_"))
            if len(requested_set) > 1 and requested_set <= file_build_set:
                return True
        return False

    if m := _PER_BUILD_FILE_PATTERN.fullmatch(filename):
        _version, annotation_consortium, genome_build = m.groups()
        return annotation_consortium.lower() in annotation_consortia and genome_build.lower() in genome_builds
    return False


def get_latest_combo_file_urls(annotation_consortia: Iterable[str], genome_builds: Iterable[str]) -> list[str]:
    """ Download URLs from the latest compatible data release for the merged ("combo") transcript files,
        ie the ones containing every historical transcript, for the requested consortia and builds.

        annotation_consortia: eg ["refseq", "ensembl"]
        genome_builds: any of
            - a single build, eg "GRCh37", "GRCh38", "T2T-CHM13v2.0" - the per-build file
            - "all-builds" (cdot.data_release.ALL_BUILDS) - the file containing all genome builds
            - an underscore-joined combination of builds, eg "grch37_grch38" - the old combined file, or
              in current releases the all-builds file (which contains them all)

        Matching is case insensitive. Eg get_latest_combo_file_urls(["refseq"], ["GRCh38"]) or
        get_latest_combo_file_urls(["refseq", "ensembl"], ["all-builds"])
    """
    # lower case everything to be case insensitive
    annotation_consortia = {x.lower() for x in annotation_consortia}
    genome_builds = {x.lower() for x in genome_builds}

    file_urls = []
    for browser_download_url in get_latest_browser_urls():
        filename = browser_download_url.rsplit("/")[-1]
        if _filename_matches(filename, annotation_consortia, genome_builds):
            file_urls.append(browser_download_url)
    return file_urls
