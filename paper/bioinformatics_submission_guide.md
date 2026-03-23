# Bioinformatics (Oxford) — Submission Guide Summary

Source: https://academic.oup.com/bioinformatics/pages/instructions_for_authors

---

## Article Types

### Original Paper
- **Page limit**: 7 pages (~5,000 words), excluding figures and references
- **Abstract structure**: Structured with the following headings:
  - **Motivation**: Why the problem matters
  - **Results**: What the paper demonstrates
  - **Availability and Implementation**: URL, language, licence
  - **Contact**: Corresponding author email
  - **Supplementary information**: If applicable
- **Body sections**: Introduction, System and methods (or Algorithm / Implementation as appropriate), Discussion, References
- Best fit for: novel algorithms, significant new methods, substantial benchmarks

### Application Note
- **Page limit**: 4 pages (~2,600 words), or ~2,000 words + 1 figure
- **Abstract structure**:
  - **Summary**: One paragraph, no subheadings
  - **Availability and Implementation**: URL, language, licence
  - **Contact**: Corresponding author email
  - **Supplementary Information**: If applicable
- Best fit for: new software tools, databases, web servers with demonstrated utility
- cdot is likely an **Application Note** (tool + data resource), but the benchmarks may push it toward an Original Paper

---

## Key Rules

- Manuscripts **>20% over the page limit** are returned without review — do not exceed
- A **cover letter is required**; omission may lead to editorial rejection without review
- Formatting is **not required at submission** — only needed at revision after acceptance
- Supplementary material should be **referenced in the abstract** if included
- All code must be **openly available** (GitHub link required in abstract)
- Data must be **publicly accessible** (Zenodo, FigShare, or equivalent for files too large for GitHub)

---

## Recommended Article Type for cdot

**Application Note** is the most natural fit: cdot is a software tool + data resource with a REST API and Python integrations. However, given the breadth of benchmarks planned (ClinVar coverage, speed comparison, T2T novelty, MANE canonical selection), an **Original Paper** is defensible and gives more room for the methods and results sections.

**Recommendation**: Draft as Original Paper (7 pages). If the content fits comfortably in 4 pages after trimming, resubmit as Application Note. Do not exceed the page limit.

---

## Abstract Template (Original Paper)

```
Motivation: [Why HGVS coordinate conversion is difficult; the UTA/infrastructure problem]

Results: [What cdot provides: 1.3M+ transcripts, JSON, REST, speed, T2T, MANE]

Availability and Implementation: https://github.com/SACGF/cdot; pip install cdot; MIT licence; data at [zenodo/cdotlib.org]

Contact: [email]

Supplementary information: Supplementary data are available at Bioinformatics online.
```

---

## Cover Letter Checklist

- [ ] Summarise the contribution in 2–3 sentences
- [ ] State why Bioinformatics is the right venue (tool paper, bioinformatics community)
- [ ] State that the manuscript has not been submitted elsewhere
- [ ] Confirm all authors agree to submission
- [ ] Suggest 3–4 potential reviewers (optional but helpful)

---

## References

- Use numbered references, cited in order of appearance: `[1]`, `[2]`, etc.
- Reference format: Author(s). Title. *Journal* Year;Volume:Pages. doi:...
- Software citations: cite the primary paper, not just the URL
- All PMIDs collected in `literature/literature_review.md` are ready for formatted reference list

---

## Figures

- Submit as separate high-resolution files at revision (not required at submission)
- Figures embedded in manuscript PDF for review is acceptable
- Colour figures: no extra charge (online-only colour is free)
- Each figure needs a descriptive legend that is self-contained (reader should understand the figure without reading the main text)
- Aim for ≤5 main figures; additional figures go to Supplementary

---

## Data Availability

Include a **Data Availability** statement:
> The cdot JSON.gz data files are available at [cdotlib.org / Zenodo DOI]. All data are versioned and publicly accessible without registration.
