import re
import requests
import cdot

from cdot import get_data_schema_int


def get_latest_data_release_tag_name():
    latest_data_release = get_latest_data_release()
    return latest_data_release.get('tag_name')

def _get_version_from_tag_name(tag_name, data_version=False):
    """ Returns None if doesn't match required prefix """
    release_prefix = "v"
    if data_version:
        release_prefix = "data_" + release_prefix

    if not tag_name.startswith(release_prefix):
        return None
    return tag_name.lstrip(release_prefix)


def get_latest_data_release():
    client_data_schema = get_data_schema_int(cdot.__version__)

    url = "https://api.github.com/repos/SACGF/cdot/releases"
    response = requests.get(url)
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
        return release
    return {}


def get_latest_browser_urls():
    if latest_data_release := get_latest_data_release():
        for asset in latest_data_release["assets"]:
            yield asset["browser_download_url"]


def get_latest_combo_file_urls(annotation_consortia, genome_builds):
    # lower case everything to be case insensitive
    annotation_consortia = {x.lower() for x in annotation_consortia}
    genome_builds = {x.lower() for x in genome_builds}

    file_urls = []
    for browser_download_url in get_latest_browser_urls():
        filename = browser_download_url.rsplit("/")[-1]
        if m := re.match(r"cdot-(\d+\.\d+\.\d+)\.(refseq|ensembl)\.(.+)\.json\.gz", filename):
            _version, annotation_consortium, genome_build = m.groups()
            if annotation_consortium.lower() in annotation_consortia and genome_build.lower() in genome_builds:
                file_urls.append(browser_download_url)
    return file_urls