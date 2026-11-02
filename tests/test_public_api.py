""" The curated public API: __all__ in each package/module names things that exist, the providers
    are importable from the package, and a missing pysam gives a helpful error. """
import importlib
import sys
import unittest

import cdot
import cdot.hgvs
import cdot.hgvs.dataproviders
from cdot.hgvs.dataproviders import fasta_seqfetcher, json_data_provider, seqfetcher


class TestAllNamesResolve(unittest.TestCase):
    MODULES = [cdot.hgvs, cdot.hgvs.dataproviders, json_data_provider, seqfetcher, fasta_seqfetcher]

    def test_all_names_resolve(self):
        for module in self.MODULES:
            with self.subTest(module=module.__name__):
                self.assertTrue(module.__all__)
                for name in module.__all__:
                    self.assertTrue(hasattr(module, name), f"{module.__name__}.__all__ names missing {name!r}")

    def test_no_duplicates(self):
        for module in self.MODULES:
            with self.subTest(module=module.__name__):
                self.assertEqual(len(module.__all__), len(set(module.__all__)))

    def test_package_all_covers_submodule_all(self):
        """ Everything a submodule declares public is reachable from the package """
        package_all = set(cdot.hgvs.dataproviders.__all__)
        for module in (json_data_provider, seqfetcher, fasta_seqfetcher):
            with self.subTest(module=module.__name__):
                self.assertTrue(set(module.__all__) <= package_all)

    def test_star_import_is_only_public_names(self):
        namespace = {}
        exec("from cdot.hgvs.dataproviders import *", namespace)
        imported = {k for k in namespace if not k.startswith("__")}
        self.assertEqual(imported, set(cdot.hgvs.dataproviders.__all__))
        for incidental in ("abc", "requests", "re", "IntervalTree", "OrderedDict", "msgspec", "FastaFile"):
            self.assertNotIn(incidental, imported)

    def test_module_docstrings(self):
        for module in (cdot, cdot.hgvs, cdot.hgvs.dataproviders):
            with self.subTest(module=module.__name__):
                self.assertTrue(module.__doc__ and module.__doc__.strip())


class TestProvidersImportable(unittest.TestCase):
    def test_ensembl_tark_data_provider(self):
        from cdot.hgvs.dataproviders import EnsemblTarkDataProvider, EnsemblTarkSeqFetcher
        from cdot.hgvs.dataproviders.ensembl_tark_data_provider import EnsemblTarkDataProvider as Direct
        self.assertIs(EnsemblTarkDataProvider, Direct)
        self.assertIn("EnsemblTarkDataProvider", cdot.hgvs.dataproviders.__all__)
        self.assertTrue(EnsemblTarkSeqFetcher)

    def test_consumers_imports_still_work(self):
        """ The names README, docs, tests and cdot_rest import from the package """
        from cdot.hgvs.dataproviders import (  # noqa: F401
            AbstractJSONDataProvider, ChainedSeqFetcher, FastaSeqFetcher, JSONDataProvider,
            LocalDataProvider, RESTDataProvider, get_ac_name_map, get_name_ac_map)

    def test_provider_docstrings(self):
        for cls in (json_data_provider.JSONDataProvider, json_data_provider.RESTDataProvider):
            with self.subTest(cls=cls.__name__):
                for word in ("mode", "cache", "seqfetcher"):
                    self.assertIn(word, cls.__doc__)


class TestFastaSeqFetcherWithoutPysam(unittest.TestCase):
    """ Re-import cdot.hgvs.dataproviders with pysam made unimportable, then put everything back """

    PACKAGE = "cdot.hgvs.dataproviders"
    REPLACED = [PACKAGE, PACKAGE + ".fasta_seqfetcher"]
    BLOCKED = ["pysam", "pysam.libcfaidx"]

    def setUp(self):
        self.saved = {name: sys.modules.get(name) for name in self.REPLACED + self.BLOCKED}
        self.saved_attr = cdot.hgvs.dataproviders
        for name in self.REPLACED:
            sys.modules.pop(name, None)
        for name in self.BLOCKED:
            sys.modules[name] = None  # makes "import pysam" raise ModuleNotFoundError

    def tearDown(self):
        for name, module in self.saved.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module
        cdot.hgvs.dataproviders = self.saved_attr

    def test_stub_raises_helpful_import_error(self):
        package = importlib.import_module(self.PACKAGE)
        self.assertIsNot(package, self.saved_attr)
        self.assertIn("FastaSeqFetcher", package.__all__)
        self.assertNotIn("GenomeFastaSeqFetcher", package.__all__)
        for name in package.__all__:
            self.assertTrue(hasattr(package, name))
        # Importing is fine, the error comes on construction
        with self.assertRaisesRegex(ImportError, r"FastaSeqFetcher needs pysam: pip install cdot\[fasta\]"):
            package.FastaSeqFetcher("genome.fa")
        # The providers are unaffected
        self.assertTrue(package.JSONDataProvider)
        self.assertTrue(package.RESTDataProvider)
        self.assertTrue(package.EnsemblTarkDataProvider)


class TestFastaSeqFetcherOtherImportErrorNotMasked(unittest.TestCase):
    """ Only a missing pysam is turned into the stub; any other ImportError from fasta_seqfetcher
        is a real problem and must surface as itself """

    PACKAGE = "cdot.hgvs.dataproviders"
    REPLACED = [PACKAGE, PACKAGE + ".fasta_seqfetcher"]

    def setUp(self):
        self.saved = {name: sys.modules.get(name) for name in self.REPLACED + ["msgspec"]}
        self.saved_attr = cdot.hgvs.dataproviders
        for name in self.REPLACED:
            sys.modules.pop(name, None)
        sys.modules["msgspec"] = None

    def tearDown(self):
        for name, module in self.saved.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module
        cdot.hgvs.dataproviders = self.saved_attr

    def test_other_import_error_propagates(self):
        with self.assertRaises(ImportError) as cm:
            importlib.import_module(self.PACKAGE)
        self.assertEqual(cm.exception.name, "msgspec")


if __name__ == "__main__":
    unittest.main()
