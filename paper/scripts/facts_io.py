"""The paper's facts convention: which type each field name implies.

Every script writes its facts with ``FACTS.write(path, facts)``, which writes a vertical
(``field,value,type``) CSV, so vibepaper checks each value against its type when it loads
it, and ``pct`` knows a percent from a fraction. The type comes from the field name:

    *_pct                       percent (already multiplied by 100)
    *_ratio, *_fold             ratio
    *_safety                    fraction (0-1)
    n_*, *_n, *_count, *_pairs  count
    build, *_version, *_date    text

Any other field holding a non-negative integer is a count too (``indel_mismatch``,
``refseq_grch38_releases``), unless its name marks a rate or a duration (``*_tps*``,
``*_s``, ``*_min``, ``*_max``, ``*_days``, ``*_pp``). Everything else is written untyped.
A None or NaN value is an error. Read another script's facts with
``vibepaper.facts.read_facts``.
"""

import numbers

from vibepaper.facts import Count, FactsWriter, Fraction, Percent, Ratio, Text

# Name parts that mark a measurement rather than a count, so an integer value stays untyped.
NOT_COUNT_SUFFIXES = ("_s", "_min", "_max", "_days", "_pp")


class PaperFactsWriter(FactsWriter):
    def type_for(self, field, value=None):
        cls = super().type_for(field, value)
        if cls is None and is_count_value(field, value):
            return Count
        return cls


def is_count_value(field, value):
    """Whether ``value`` is a non-negative integer that counts something."""
    if isinstance(value, bool) or not isinstance(value, numbers.Integral) or value < 0:
        return False
    return "_tps" not in field and not field.endswith(NOT_COUNT_SUFFIXES)


FACTS = PaperFactsWriter({
    "*_pct": Percent,
    "*_ratio": Ratio, "*_fold": Ratio,
    "*_safety": Fraction,
    "n_*": Count, "*_n": Count, "*_count": Count, "*_pairs": Count,
    "build": Text, "*_version": Text, "*_date": Text,
})
