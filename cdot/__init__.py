""" cdot: transcript data for the biocommons hgvs and PyHGVS variant annotation libraries.

    cdot converts RefSeq and Ensembl GTF/GFF3 annotation into JSON files (released on GitHub, see
    :mod:`cdot.data_release`) and serves them to the HGVS libraries as data providers:

    - :mod:`cdot.hgvs.dataproviders`: biocommons hgvs providers (``JSONDataProvider`` for local
      files, ``RESTDataProvider`` for cdotlib.org, ``EnsemblTarkDataProvider``) and sequence fetchers
    - :mod:`cdot.hgvs`: fixing real-world HGVS input (``clean_hgvs``, ``fix_hgvs``, ``fix_c_to_p``)
    - :mod:`cdot.models`: the typed JSON data format, and ``load()`` for reading a file
    - :mod:`cdot.pyhgvs`: the (legacy) PyHGVS transcript factories

    ``__version__`` is the client version released to PyPI. The data files carry their own version;
    a client can read data whose schema (major.minor) is not newer than its own, see
    :func:`get_data_schema_int`.
"""
__version__ = "0.2.37"
# Data version is kept in generate_transcript_data/json_schema_version.py

def get_data_schema_int(version: str) -> int:
    """ Return an int which increments upon breaking changes - ie anything other than patch """
    major, minor, patch = version.split(".")
    return 1000 * int(major) + int(minor)
