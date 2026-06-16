# Design notes & project direction

Background on *why* the cdot JSON format looks the way it does, and where the project is heading. For
the concrete format, see the [JSON data format reference](json_data_format.md) and
[Coordinates & exon alignments](coordinates_and_exons.md).

## Why JSON?

In bioinformatics, gene/transcript information is traditionally stored in GTF/GFF — good for genome
browsers and manipulable with unix tools, but very verbose and not easy to use.

JSON has become the standard for the web, and Python has very fast implementations — so storing
gene/transcript information in JSON seems like a good thing to do. Ideally we'd like to make it a
standard, and will pitch GA4GH about it.

The transcript format was originally built for the [PyReference](https://github.com/SACGF/pyreference)
project, then adapted to PyHGVS, and finally to biocommons HGVS and cdot (REST). Having multiple
consumers keeps us honest and portable.

There are example transcripts on the REST API at [cdotlib.org](https://cdotlib.org). An example JSON
file is in the repo:
[`tests/test_data/cdot.refseq.grch37.json`](../tests/test_data/cdot.refseq.grch37.json).

## Known issues / potential changes

* The `exons` array tracks transcript start/end; this could be rebuilt from the alignment gaps.
* The `exons` array uses `gap=None` to indicate a perfect alignment. biocommons HGVS always provides
  e.g. `M100` for a perfect 100-base match.
* The `exons` array contains `exon_id`; this could be derived at runtime via reverse/enumerate.
* Genome build patches may have base changes that alter splice sites — this may make historical GTFs
  subtly wrong.
* In `cdot_json` `merge_builds`, we throw away earlier coordinates when there's a conflict — worth
  investigating further (likely build-sequence changes altering exons).
* We store coordinates per genome build; storing by contig would remove redundancy for shared contigs
  (e.g. chrM), at the cost of being slightly harder for humans to read.

## Other examples in the wild

Ensembl provide a [lookup API](https://rest.ensembl.org/documentation/info/lookup), e.g.
[ENSG00000179348](https://rest.ensembl.org/lookup/id/ENSG00000179348?expand=1;content-type=application/json),
which looks to have enough to create HGVS records — however Ensembl do not provide all historical
transcript versions.

## Project goals / directions

This project was spun out of the Australian Genomics [Shariant](https://shariant.org.au) project,
which involves resolving legacy HGVS from many different labs. We feel there is a need to obtain a
JSON representation of historical transcript versions for use in resolving HGVS, and this project can
fill that gap.

Long term, we'd like to make it **obsolete**. That would require:

* Perhaps moving cdot into biocommons (if they'll have us), or the client code into the HGVS projects.
* Proposing a GA4GH transcript JSON format.
* Having RefSeq/Ensembl provide APIs for historical transcript versions (ideally in the GA4GH standard
  above).

We'll try to do the above, and as things change, we can be the glue that ties it together.
