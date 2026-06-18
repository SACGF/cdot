# Discussion

*Target: ~300 words. Positioning, limitations, future. Do not repeat implementation
detail.*

---

cdot fills a specific gap in the biocommons HGVS stack: comprehensive, versioned,
multi-source transcript coordinate data accessible without database infrastructure.
Together with SeqRepo [@Hart2020] for biological sequences, cdot and biocommons/hgvs
form a complete offline HGVS processing stack.

Compared to Ensembl TARK, the other REST-based transcript archive, cdot provides
RefSeq coverage, offline capability, alignment gap storage, and support for
T2T-CHM13v2.0. VariantValidator [@Freeman2018], which provides web-based HGVS validation
built on biocommons/hgvs, is complementary: cdot could serve as
VariantValidator's data backend, combining cdot's coverage with VariantValidator's
validation logic.

