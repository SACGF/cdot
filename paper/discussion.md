# Discussion

*Target: ~300 words. Positioning, limitations, future. Do not repeat implementation
detail.*

---

cdot fills a specific gap in the biocommons HGVS stack: comprehensive, versioned,
multi-source transcript coordinate data accessible without database infrastructure.
Together with SeqRepo [@Hart2020] for biological sequences, cdot and biocommons/hgvs
form a complete offline HGVS processing stack.

cdot also integrates Ensembl TARK, the Ensembl transcript archive: its
`EnsemblTarkDataProvider` is, to our knowledge, the only client that exposes TARK through
the biocommons/hgvs data-provider interface, letting existing HGVS pipelines query TARK
without bespoke code. Beyond what TARK's REST service offers, cdot adds RefSeq coverage,
fully offline operation, and support for T2T-CHM13v2.0. VariantValidator [@Freeman2018],
which provides web-based HGVS validation built on biocommons/hgvs, is complementary: cdot
could serve as VariantValidator's data backend, combining cdot's coverage with
VariantValidator's validation logic.

