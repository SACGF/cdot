import json
import os.path
import re
from inspect import getsourcefile
from os.path import abspath

from cdot.hgvs.dataproviders.ensembl_tark_data_provider import EnsemblTarkDataProvider


class MockEnsemblTarkDataProvider(EnsemblTarkDataProvider):
    def __init__(self, assemblies: list[str] = None, mode=None, cache=None, seqfetcher=None):
        super().__init__(assemblies, mode, cache, seqfetcher)
        self._this_file_dir = os.path.dirname(abspath(getsourcefile(lambda: 0)))


    def _get_from_url(self, url):
        if not url.startswith(self.base_url):
            raise ValueError(f"{url} does not start with {self.base_url}")

        dirname = os.path.dirname(url)
        basename = os.path.basename(url)
        params = re.sub(r"^\?", "", basename)
        path = re.sub(f"^{self.base_url}/", "", dirname)
        filename = os.path.join(self._this_file_dir, "test_data", "ensembl_tark", path, f"{params}.json")
        if not os.path.exists(filename):
            raise FileNotFoundError(f"{filename} not found")

        with open(filename, "r") as f:
            return json.load(f)



