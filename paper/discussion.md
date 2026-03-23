# Discussion

*Target: ~300 words. Positioning, limitations, future. Do not repeat implementation detail.*

---

cdot fills a specific gap in the biocommons HGVS stack: comprehensive, versioned, multi-source transcript coordinate data accessible without database infrastructure. The analogous role for biological sequences is played by SeqRepo [Hart 2020]; together, cdot + SeqRepo + biocommons/hgvs constitute a complete offline HGVS processing stack.

Compared to Ensembl TARK — the other REST-based transcript archive — cdot provides RefSeq coverage, offline capability, alignment gap storage, and support for T2T-CHM13v2.0. VariantValidator [Freeman 2018], which provides web-based HGVS validation built on biocommons/hgvs, is complementary rather than competing: cdot could serve as VariantValidator's data backend, combining cdot's coverage with VariantValidator's validation logic.

The current version has known limitations. CDS phase at the exon level is not modelled, affecting protein HGVS for a small number of transcripts with ribosomal frameshifts. Ensembl GFF3 files contain unreliable protein version information; cdot uses GTF files for Ensembl data as a workaround. cdot provides coordinate data only — HGVS syntax validation should be performed by a dedicated tool (Mutalyzer [Lefter 2021], VariantValidator [Freeman 2018]).

Looking ahead, the Human Pangenome Reference [Liao 2023] and increasing use of haplotype-resolved assemblies motivate cdot's multi-build architecture: additional reference assemblies can be added without format changes. VEP-compatible data releases would allow cross-referencing between VEP annotations and HGVS-specific tools using matched transcript versions, directly addressing the tool-disagreement problem documented by McCarthy et al. [2014].

---

*If over word count, cut the Pangenome/VEP sentence — move to supplementary or drop entirely.*
