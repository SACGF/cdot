The paper has lost it's "North star" - the goal to "resolve as many HGVS as possible" - which is the central part of the story - and now reads like a mismash of cleaning, transcript bumping and transcript provision. If you think that shouldn't be the title, then it should definitely be in the opening "motivation" section - which also doesn't mention Shariant. So it doesn't actually contain the motivation

This perhaps needs a whole paper scan to see that we have proper framing so the story hangs together

ABSTRACT

Motivation - UTA doesn't require a database - only for local, there is public one. so tradeoff is local install complexity vs slowness using internet resources (and it's firewalled off for many users who are restricted to http(s))

results - "nearly 4 orders of magnitude" - we should only ever compare like for like ie remote/remote vs local/local it is unfair to UTA to do our local vs their remote

Availability - the JSON.gz files are on github page (releases)

INTRODUCTION

1st paragraph too long. "Transcript choice also matters" - not sure this is relevant.

Do we need to cite both biocommons HGVS papers, or just latest?

METHODS

"because a NM version retired years ago is still cited..." - this is rehash of motivation doesn't belong in methods

Json format is very long, feel could be a sentence or two and I think we mention it a few times later.

"Version fallback" May be best referred to as version substitition? We also call it version bumping a few times - I feel substitution is best

One of the questions is - WHY do we need version fallback - we have gotten every file off the FTP sites - so what explains it? My suspicion is that perhaps those labs are using tools that align transcripts to genomes themselves - it may be worth seeing if transcripts we don't have come from only a small amount of labs, or what to try and work out what's going on. I think this is a new analysis

benchmarking section seems too long

RESULTS

Instead of Tier 1 / Tier 2 maybe just call them public data/private data?

"Ensembl VEP can already generate HGVS" - does this belong in methods? Also - Ensembl VEP don't handle RefSeq alignment gaps - known issue documented on their repo.

R2 - using Ensembl in ClinVar test is not fair on UTA, no point and it conflates not having Ensembl vs missing historical versions. 500 refseq and 818 Ensembl? WTF? That does not match the reality of what is in clinvar. I think this whole analysis is using the OLD ClinVar VCF not the VCV XLM - and needs to be updated to use a sample of historical records (fair eg taken from time buckets, then maybe second run that just uses whole file random sampling - which will have recency bias due to growth over time)

We don't mention why we didn't do UTA on whole of ClinVar (would have been far too slow) - I think we could possibly do it on the local UTA though? Only take a few hours?

"What laboratories actually submit" is a wacky way to talk about historical submissions - it's not about what they do now but what they did in the past

You say "Not one submitted string cites an Ensembl transcript" is wrong, and contradicts that you used 818 Ensembl transcripts a few paragraphs before

R3 Throughput - cdot REST still faster

The results discusses the UTA vs cdot differences - not sure if should be in this section - but if you are discussing it - the main reason I think is that cdot only needs 1 json object to get all the info it needs to perform transcript coordinate projection - while UTA requires multiple queries. Maybe this belongs in JSON section or discussion.

You say "sequence cache conditions differed" - maybe we should just do 1 run then discard it to have hot cache (report that time and how it is only a tiny bit slower)

R4

"typed into the HGVS search box" - to be pedantic it's a search box and we are taking the subset of strings that matched a broad HGVS regex (ie were loosely hgvs - so could be repaired)

Residual errors - "non HGVS input - pasted urls etc" - I feel this is a data collection error, and we should remove these from the test corpus before counting ith, as it doesn't reflect tool correctness - analysis change I think

"version bump can move a variant" - might be too slangy? - perhaps change what coordinate projected between trascript/genome. 

DISCUSSION

Initial paragraph should be about general source of transcripts I think - all of them, easy to parse, general use.  

Ensembl TARK section too big. Just offhand mention required






