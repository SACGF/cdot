"""Write and read the paper's facts CSVs.

A facts file is vertical (``field,value,type``), so vibepaper checks each value against
its type when it loads it, and ``pct`` knows a percent from a fraction. The type comes
from the field name, a convention every script here follows:

    *_pct                       percent (already multiplied by 100)
    *_ratio, *_fold             ratio
    *_safety                    fraction (0-1)
    n_*, *_n, *_count, *_pairs  count
    build, *_version, *_date    text

Any other field is written untyped.
"""

import csv

from vibepaper.facts import Count, Fraction, Percent, Ratio, Text, write_facts


def fact_type(field):
    """The vibepaper type for ``field``, or None for an untyped fact."""
    if field.endswith("_pct"):
        return Percent
    if field.endswith(("_ratio", "_fold")):
        return Ratio
    if field.endswith("_safety"):
        return Fraction
    if field.startswith("n_") or field.endswith(("_n", "_count", "_pairs")):
        return Count
    if field == "build" or field.endswith(("_version", "_date")):
        return Text
    return None


def write_typed_facts(path, facts):
    """Write ``facts`` ({field: value}) to ``path`` as a vertical facts CSV."""
    typed = {}
    for field, value in facts.items():
        cls = fact_type(field)
        typed[field] = cls(value) if cls else value
    return write_facts(path, typed)


def read_facts(path):
    """A vertical facts CSV as {field: value}, with the values as written in the file."""
    with open(path, newline="") as fh:
        return {row["field"]: row["value"] for row in csv.DictReader(fh)}
