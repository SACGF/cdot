# Figure planning scratchpad

Working notes for the paper figures. Iterate on the text graphs here first;
only regenerate the real figure (`figure1.svg`) once the structure is settled.
A draft rendering of Option A already exists at `paper/figures/figure1.svg`
(open in a browser to view).

## Constraints

- Bioinformatics: 4 pages, ~2,000 words + 1 figure, or ~2,600 words with none.
  README targets 2 main figures, so space is tight; anything that doesn't earn
  its place goes to supplementary.
- Numbers drawn into a figure (release counts, 1.62 M, var/s) are frozen at
  draw time; they will NOT update with the facts pipeline. Either keep numbers
  out of the figure, or accept re-drawing on each data refresh.

## Options considered

**A. One two-panel architecture figure (current draft)**
Panel A: build-time data generation (sources -> pipeline -> JSON.gz -> release/REST).
Panel B: client-side resolution (HGVS string -> clean -> biocommons/hgvs -> g.),
with providers + sequence layer feeding the library, joined to panel A by two
vertical "download"/"HTTPS" arrows.
Pro: one figure covers the whole system; the vertical arrows show the key claim
(same JSON drives local and REST). Con: dense; the novel string-handling story
(clean/fallback) is reduced to two small boxes.

**B. Pure architecture diagram (components + interfaces, no flow)**
Class-diagram flavour: Interface at centre, providers around it.
Con: restates the module list; a cold reader learns less than from prose. Weak.

**C. Flow chart of resolving one HGVS string**
Lifecycle with decision points: parse fail -> clean_hgvs -> reparse;
version missing -> get_best_transcript_version -> safety check -> retry;
gene symbol only -> MANE lookup. This is the paper's novel contribution
(R4/R5) and makes a good story figure, but hides where the data comes from.

**D. Two figures: A-style architecture as Fig 1 + C-style resolution flow as Fig 2**
Fig 1 stays high level (shrink panel B: drop helper boxes, keep providers),
Fig 2 gets the decision logic with the fix/fallback loop. Matches the
"target: 2 main figures" note. Risk: overlap between Fig 1 panel B and Fig 2;
keep Fig 1 panel B strictly structural (who talks to whom) and Fig 2 strictly
behavioural (what happens to one string) to avoid it.

## Decisions so far (2026-08-13)

- **Two figures** (Option D): Fig 1 architecture/dataflow, Fig 2 resolution flow.
- **Node encoding**: distinguish *processes/code* from *data/files*.
  Proposal: shape carries type (rounded rect = process/code, document/cylinder
  shape = data/file), colour carries ownership (blue = cdot, grey = external).
  In the text graphs below: `[square brackets] = process`, `(parens) = data`.
- **Drop the "data release" node**: JSON.gz is simply hosted on GitHub; say so
  inside the JSON.gz node. Drop the Zenodo mention (note: abstract still has a
  "[Zenodo DOI]" placeholder; revisit for consistency at submission).
- **JSON.gz node mentions subset files**: combined file plus subsets per
  source (RefSeq/Ensembl), per build, per release.
- **Shape vocabulary (v3)**: rounded rect = process/code; dog-ear (folded
  corner) = file; cylinder = database (UTA, SeqRepo); parallelogram = plain
  string input/output (HGVS in, converted HGVS out). No wavy-bottom document
  shape.
- **Output node is bidirectional**: "converted HGVS, c. -> g. or g. -> c.",
  not just genomic coordinates.
- **cdot_rest arrow** labelled "loaded into Redis".
- **Panel B flipped**: the groups that feed biocommons/hgvs (providers,
  sequence layer, helpers) sit at the TOP of panel B, directly under panel A,
  so the download/HTTPS arrows are short. The string-resolution flow runs
  left-to-right along the BOTTOM. Whole figure then reads top-to-bottom:
  sources -> data -> providers -> resolution.

## Text graph: Fig 1 v2 (architecture / dataflow)

```
PANEL A - data generation (build time)

 (RefSeq GFF3: 35 hist. releases)--\
 (Ensembl GTF: 40 rel. 81 to 115)--+-> [Snakemake pipeline]--> (cdot JSON.gz)      ---loads--->  [cdot_rest]
 (UTA: algns missing from GTF)-----/    parse (HTSeq)           1.62M alignments                  REST API
                                        normalise contigs       GRCh37/38/T2T                     cdotlib.org
                                        CDS, gene info          hosted on GitHub
                                        merge, newest wins      + subset files
                                                                (source/build/release)
========================================================================|=========================|=========
PANEL B - variant resolution (client)                                   | download                | HTTPS
                                                                        v                         v
 [opt-in helpers]          [sequence layer]          [cdot data providers (biocommons/hgvs Interface)]
  resolve_gene_hgvs()       ChainedSeqFetcher:        JSONDataProvider | EnsemblTark    | RESTDataProvider
   gene symbol -> MANE       (SeqRepo)                 in-memory,      | DataProvider   |  per-transcript
  get_best_transcript_       [FastaSeqFetcher]         interval trees  | (TARK service) |  + prefetch()
   version()                  <- (genome FASTA)              \               |               /
   fallback + safety check         |                          \    exons + alignments      /
       |                           | transcript                v             v            v
       | rewrites input            |   sequence          +--------------------------------+
       v                           +-------------------> |                                |
 (HGVS string) ---> [clean_hgvs()] --------------------> [biocommons/hgvs] ---> (genomic (g.) coordinates)
  possibly           repair +                             parse, validate,
  malformed          HGVSFix audit                        map c. <-> g.

 footnote row: also JSONPyHGVSTranscriptFactory -> PyHGVS (legacy)
```

## Text graph: Fig 2 candidate (resolution flow, Option C/D)

```
              gene-symbol HGVS? --yes--> resolve_gene_hgvs() -> MANE transcript, rewrite
                     |no
  input string -> parse OK? --no--> clean_hgvs() -> HGVSFix records -> reparse
                     |yes                              (ERROR-level fix -> report, stop)
                     v
        transcript version in data? --no--> get_best_transcript_version()
                     |yes                    -> candidate version -> CDS safety check
                     v                          safe? --no--> report, stop
              map c. -> g.                       |yes -> rewrite + retry
              (provider: exons/alignments;
               seqfetcher: transcript seq)
                     v
              genomic coordinates + audit trail of every change
```

## Open questions

1. ~~One figure or two?~~ DECIDED: two.
2. Numbers in the figure (35/40 releases, 1.62 M, 540-665 var/s): keep, round
   harder ("~1.6 M"), or move all numbers to the caption?
3. How prominent should PyHGVS be? (currently a one-line footnote)
4. EnsemblTarkDataProvider: keep in Fig 1, or footnote it like PyHGVS? It has
   no arrow from panel A (talks to TARK directly), which clutters the story.
5. Does panel A need the pipeline-step detail (5 bullets), or just
   "parse + normalise + merge, newest wins"?
6. Sequence layer: essential to the figure's claim, or caption material?
7. Shape encoding for data vs process: document shape, cylinder (DB-ish), or
   just square vs rounded corners? (Cylinders may read as "database" which is
   wrong for flat files; document shape is the conventional flowchart glyph.)
8. cdot_rest arrow: with the release node gone, is "loads" from JSON.gz on
   GitHub accurate, or does cdot_rest import into its own store? (affects
   the arrow label only)
