import unittest
from unittest.mock import patch, MagicMock

import cdot
from cdot import data_release

# Real asset names from the data_v0.2.36 release (current naming scheme)
CURRENT_ASSETS = [
    "cdot-0.2.36.all-builds-ensembl-grch37_grch38_t2t-chm13v2.0.json.gz",
    "cdot-0.2.36.all-builds-refseq-grch37_grch38_t2t-chm13v2.0.json.gz",
    "cdot-0.2.36.ensembl.GRCh37.json.gz",
    "cdot-0.2.36.ensembl.GRCh38.json.gz",
    "cdot-0.2.36.ensembl.T2T-CHM13v2.0.json.gz",
    "cdot-0.2.36.Homo_sapiens_GRCh37_Ensembl_87.gtf.json.gz",
    "cdot-0.2.36.Homo_sapiens_GRCh37_RefSeq_105.20220307.gff.json.gz",
    "cdot-0.2.36.Homo_sapiens_GRCh38_Ensembl_116.gtf.json.gz",
    "cdot-0.2.36.Homo_sapiens_GRCh38_RefSeq_RS_2025_08.gff.json.gz",
    "cdot-0.2.36.refseq.GRCh37.json.gz",
    "cdot-0.2.36.refseq.GRCh38.json.gz",
    "cdot-0.2.36.refseq.T2T-CHM13v2.0.json.gz",
]

# Older naming scheme: lowercase builds, single GTF files named after the source, combined 37/38 file
OLD_ASSETS = [
    "cdot-0.2.25.ensembl.grch37.json.gz",
    "cdot-0.2.25.ensembl.grch38.json.gz",
    "cdot-0.2.25.ensembl.grch37_grch38.json.gz",
    "cdot-0.2.25.ensembl.Homo_sapiens.GRCh38.111.gff3.json.gz",
    "cdot-0.2.25.ensembl.T2T-CHM13v2.0.json.gz",
    "cdot-0.2.25.GCF_000001405.40_GRCh38.p14_genomic.RS_2023_10.gff.json.gz",
    "cdot-0.2.25.refseq.grch37.json.gz",
    "cdot-0.2.25.refseq.grch38.json.gz",
    "cdot-0.2.25.refseq.grch37_grch38.json.gz",
    "cdot-0.2.25.refseq.T2T-CHM13v2.0.json.gz",
]

DOWNLOAD_BASE = "https://github.com/SACGF/cdot/releases/download"


def _data_version(patch_version: int) -> str:
    """ A data version with the same schema (major.minor) as this client """
    major, minor, _ = cdot.__version__.split(".")
    return f"{major}.{minor}.{patch_version}"


def _mock_releases_response(asset_names):
    tag = f"data_v{_data_version(99)}"
    releases = [
        {"tag_name": "v9.9.9", "assets": []},  # code release - skipped
        {"tag_name": "data_v999.0.0", "assets": []},  # incompatible schema - skipped
        {
            "tag_name": tag,
            "assets": [{"browser_download_url": f"{DOWNLOAD_BASE}/{tag}/{name}"} for name in asset_names],
        },
    ]
    response = MagicMock()
    response.json.return_value = releases
    return response


def _filenames(urls):
    return sorted(url.rsplit("/", 1)[-1] for url in urls)


class TestDataRelease(unittest.TestCase):
    def _combo(self, asset_names, annotation_consortia, genome_builds):
        with patch("cdot.data_release.requests.get", return_value=_mock_releases_response(asset_names)):
            return _filenames(data_release.get_latest_combo_file_urls(annotation_consortia, genome_builds))

    def test_latest_release_version(self):
        with patch("cdot.data_release.requests.get", return_value=_mock_releases_response(CURRENT_ASSETS)):
            self.assertEqual(data_release.get_latest_release_version(), _data_version(99))

    def test_current_per_build(self):
        self.assertEqual(self._combo(CURRENT_ASSETS, ["refseq"], ["grch38"]),
                         ["cdot-0.2.36.refseq.GRCh38.json.gz"])
        self.assertEqual(self._combo(CURRENT_ASSETS, ["RefSeq", "Ensembl"], ["GRCh37", "T2T-CHM13v2.0"]),
                         ["cdot-0.2.36.ensembl.GRCh37.json.gz",
                          "cdot-0.2.36.ensembl.T2T-CHM13v2.0.json.gz",
                          "cdot-0.2.36.refseq.GRCh37.json.gz",
                          "cdot-0.2.36.refseq.T2T-CHM13v2.0.json.gz"])

    def test_current_all_builds(self):
        expected = ["cdot-0.2.36.all-builds-refseq-grch37_grch38_t2t-chm13v2.0.json.gz"]
        self.assertEqual(self._combo(CURRENT_ASSETS, ["refseq"], [data_release.ALL_BUILDS]), expected)
        self.assertEqual(self._combo(CURRENT_ASSETS, ["refseq"], ["grch37_grch38_t2t-chm13v2.0"]), expected)
        # An old style combined request is served by the all-builds file
        self.assertEqual(self._combo(CURRENT_ASSETS, ["refseq"], ["grch37_grch38"]), expected)

    def test_current_all_builds_both_consortia(self):
        self.assertEqual(self._combo(CURRENT_ASSETS, ["refseq", "ensembl"], ["all-builds"]),
                         ["cdot-0.2.36.all-builds-ensembl-grch37_grch38_t2t-chm13v2.0.json.gz",
                          "cdot-0.2.36.all-builds-refseq-grch37_grch38_t2t-chm13v2.0.json.gz"])

    def test_old_naming(self):
        self.assertEqual(self._combo(OLD_ASSETS, ["refseq"], ["GRCh38"]),
                         ["cdot-0.2.25.refseq.grch38.json.gz"])
        self.assertEqual(self._combo(OLD_ASSETS, ["refseq"], ["grch37_grch38"]),
                         ["cdot-0.2.25.refseq.grch37_grch38.json.gz"])
        self.assertEqual(self._combo(OLD_ASSETS, ["refseq"], ["all-builds"]), [])

    def test_no_match(self):
        self.assertEqual(self._combo(CURRENT_ASSETS, ["refseq"], ["hg19"]), [])
        self.assertEqual(self._combo(CURRENT_ASSETS, ["gencode"], ["grch38"]), [])


if __name__ == "__main__":
    unittest.main()
