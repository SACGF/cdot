Paper feedback

ABSTRACT

I don’t like “underpinning >3,000,000 ClinVar submissions and ACMG/AMP classification guidelines” - I like the hgvs opening in intro focusing on clinical reporting - if you can reuse that without sounding like obvious duplication

We should insert that cdot does not do HGVS resolution, but instead extends existing tools and fits into the existing ecosystem (I think after “...as possible” - but wherever you think it fits)

“files loadable at 540–560 transcripts/second — no database required” - I think rewrite this to say can distribute as a gzipped file, or obtained from REST server over the internet (abstract does not currently mention REST server)

INTRODUCTION

“The clinical scale of this requirement is large: ClinVar holds >3,000,000 variants described in HGVS notation” - ClinVar is also distributed as VCF so I don’t think it is a hard requirement or HGVS underpins ClinVar - basically enough to say it’s used in reporting

“ACMG/AMP guidelines assume reliable conversion” - This is kind of strange thing to claim

PyHGVS does not use a data provider. It gives example on page of how to load transcripts (from exported UCSC data)

“Remote server frequently blocked by clinical firewalls” - should mention that it’s postgres protocol and many hospitals / workplaces block eveything but http/https

The part about “resolving as many real world transcripts as possible” bold is a bit too much, and we already say this in the title. - tone this section down reads “too excited” and kinda cringe

“Those naming transcripts” - not sure about “naming” here

“Cdot addresses these gaps” - we just said it was built for these gaps

DATA SOURCES AND GENERATION 

We also use HGNC and RefSeq gene info. Methods kind of mention HGNC at the end

JSON FORMAT
If we are discussing JSON format - should justify the “genome_builds” design with a key for genome build is used even in single build JSONs, so that the file format remains consistent whether it is 1 genome build (used as replacement for single GTF) - can combine easily from multiple builds (eg when we return GRCh37, GRCh38 and T2T from cdot rest)
—-------

Design decisions
I feel we should separate out some of the cleaning/transcript bumping etc to be always unambiguous and possibly making a wrong decision (eg transcript bumping, looking up MANE transcript)
The potentially wrong ones should probably not be on by default, and perhaps we should return some kind of easy way to check “is unambiguous” etc
In the paper, these should be grouped together in a section talking about risk etc and not being on by default
First mention of Get_best_transcript_version does not mention ambiguity really.
“Reconstruction and gene/transcript repair” - BRCA(NM_0059.4)c: - I think these may be officially supported HGVS nomenclature standard - I guess just not by the tool? Should distinguish this

Local JSON - “Throughput of 540–560 transcripts/second is 540–560× faster than UTA remote access (~1 transcript/second), comparable to SeqRepo’s 1300× speedup over remote sequence retrieval (Hart and Prlić 2020).” - we should only really compare local cdot vs UTA and REST vs public postgres - unfair otherwise - people can see the total times in the table

FastaSeqFetcher - we should discuss how “pasting together exons from the genome” works for Ensembl, but is not guaranteed to match transcript RefSeq - the HGVS c. is transcript reference (though variants are typically discovered in NGS data mapped to genomes so the “bad mapping” can go the other way typically without question) 
Canonical transcript selection - should give an example eg BRCA1 variant and the transcript conversion

Transcript coverage and ClinVar resolution - ClinVar uses the latest so it’s not a demonstration so much of providing lots of historical references, but rather a sanity check that our GTF -> JSON -> Biocommons HGVS provision works at ClinVar scale
“They are reported as frozen constants” - implementation issue readers will not care

For the bad real world HGVS - should mention that this was from a search box on a webapp, that takes people to variants/classifications - people use as shortcut to get to information quickly (background on what exactly this data is)
Fixes - “whitespace / non printable” - big issue probably from word documents, PDFS etc and copy/pasting around in clinical environment
Table 1 - fixes should all have simple examples - makes it clearer
R3 - ceiling of cleaning - CLASS is way too long, and takes up too much of table (may need code refactor) table should have counts as well as percentage. Examples would go a long way again
I feel the NON_HGVS should just be taken out of the table, as it’s probably just a bad regex on my search thinking it perhaps needs to be rescued HGVS but was just random crap/other
GENOMIC_REF_IN_PARENS - we could probably fix this?
“An LLM applying” - use the full name of the model (you!)
R4: Transcript version fallback
This doesn’t mention any numbers. As this is risky - I’d like to know how many times version were bumped out of total clinvar and how many times it mattered
Really, we should do some kind of analysis on how often a transcript version bump is going to work - could see by checking how often exon coords change - how many c. bases end up at different g. Ones by walking along transcript - not sure if anyone has ever looked at this?)

—---------
“JSON file loads in ~11s” somewhere else mentions “typically 10.1 seconds for GRCh38 RefSeq” - these should be consistent - and use same amount of decimal places I am fine with ~11s


“Orders of magnitude faster than remote server” - don’t do direct local vs remote comparisons - it is unfair, and if people want to look - table tells that story.

DISCUSSION
Motivation for Ensembl TARK is that it allows people to use Ensembl’s official REST API source

VariantValidator comparison is strange. Maybe just mention it and Mutalyzer are used to test HGVS correctness in discussion don’t suggest they should adopt us

Discussion seems short - Another use for cdot is that JSON is convenient to parse (much easier/faster than GTF) and trivial to load over the internet (REST) so we are in effect hosting an easy gene/transcript reference. We distribute 1 GTF/GFF release versions on the github release page - you can use this as a replacement for a GTF/GFF that loads much more rapidly.

The REST serves JSON so fast (and in batches) that you can distribute software without requiring large downloads. And use AI to call it - much gzipped GTF/GFF. Ensembl provide a REST service, but only Ensembl, and only latest version for each transcript. 
—---
Bibliography doesn’t seem sorted right. I think Biocommons HGVS should be #2 after HGVS nomeclature

