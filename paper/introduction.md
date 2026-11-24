# Introduction

HGVS is the clinical standard for describing variants [@DenDunnen2016]. The c. form is
typically preferred in clinical genetics because medical scientists find it easier to reason
about. In genomics, we need to convert them to genomic coordinates. To do this, we require the
reference transcript to be at a specific version. These versions of the transcript are aligned
to the genome by RefSeq and Ensembl. The exon coordinates from the transcript genome
alignments can be used to convert them to genomic coordinates.

One of the surprising things we've seen is that HGVS strings are often broken in a way that
suggests that humans have either just typed them out by hand or copy-pasted them
incorrectly. Mutalyzer looked at 26 million descriptions submitted over five years and found
only half were correct, {{ literature.mutalyzer_error_pct | pct(0) }} had errors, and it
could only fix {{ literature.mutalyzer_autocorrect_pct | pct(0) }} [@Lefter2021].

The two main Python HGVS libraries are biocommons HGVS [@Wang2018] and PyHGVS. Biocommons
gets its transcript genome data from the data provider interface, which has become a bit of
a standard, having been adopted by a number of other tools, including VariantValidator
[@Freeman2018] and MaveDB [@Rubin2025]. We found
{{ adoption.n_interface_implementers | number_word }} by searching public GitHub
repositories, so this is a lower bound. PyHGVS provides its own simple data structures.

The default implementation of the biocommons HGVS data provider is UTA, the Universal
Transcript Archive [@Hart2020]. This works by performing its own transcript to genome
alignments and storing them in a Postgres database.

This is a great resource, but has a number of drawbacks. Firstly, doing its own alignment
means that it's a large amount of effort to add new transcripts, and for instance, Ensembl
is not supported. The total number of transcript alignments is
{{ coverage.uta_count | commas }}, RefSeq only, and it provides a limited amount of
historical support. The other drawback of a PostgreSQL database is that if you set it up
locally, it requires quite a bit of configuration work. Many workplaces, clinical labs and
hospitals firewall off all connections other than HTTPS, which means that they are unable
to connect to the public Postgres instance. The protocol for communication with UTA requires
several round trips to perform the conversion, which means the performance is
{{ benchmark.uta_remote_tps | dp(2) }} HGVS/s (Table 1).

cdot was split out of the Australian Genomics Shariant project [@Tudini2022]. The goal of
Shariant is to collect variant classifications from the different labs across Australia. We
found the most common form that these variants were stored in was HGVS, but they came from
dozens of different systems. There were multiple commercial systems, custom built software
solutions, and also spreadsheets. These often went back many, many years. So one of the
great problems in Shariant was how to resolve all of these HGVSs including transcript
versions that are many years old.

We also found people entering HGVSs in the search box, looking for existing
classifications. People were using all kinds of invalid HGVSs, and we tried to resolve these
to what the people actually intended.

We had to solve a number of problems. There was the cleaning. There was also collecting as
many historical transcript versions as we could. And the other problem was some labs were on
GRCh38 and some on GRCh37, and we had to convert between them. There were often no official
releases of a transcript version against a particular genome build. So what we did then was
try and find the nearest transcript version available for that build that still had the same
exonic structure and thus was safe to substitute (opt-in and flagged).

cdot reads the official GFF or GTF releases from RefSeq and Ensembl, and then converts them
into a widely supported JSON format that allows us to collect all of the
{{ coverage.total_count | commas }} versioned transcript alignments across GRCh37, GRCh38,
and T2T genome builds [@Nurk2022]. Almost all languages have libraries that have been highly
optimized so that you can load it very rapidly. The same data can also be transferred across
the internet via a REST server, which means that workplaces with web access can retrieve them
without local copies. We provide our cleaning library, which fixes a lot of common mistakes,
giving a record of what fixes were made.

We validate and benchmark this by comparing cdot versus UTA on the historical ClinVar
submissions, the VCV XML, which records the original transcript version that was submitted
at the time, giving a historical test record. We measure performance in HGVS resolved per
second, as well as the number resolved. To test the HGVS cleaning and fixes, we take public
ClinVar data and then synthetically modify it to match the problems we found in the invalid
HGVSs from our own private search data, and we compare our fixes to other tools, including
LOVD, VariantValidator and Mutalyzer.
