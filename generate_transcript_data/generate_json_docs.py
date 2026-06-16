#!/usr/bin/env python3
"""
Generate documentation for the cdot JSON data format from the typed models in
``cdot/models.py`` (issue #37).

Produces two artifacts in ``docs/``:

* ``cdot-json-schema.json`` - a JSON Schema (language-agnostic, for validation /
  codegen in any language), enriched with per-field descriptions.
* ``json_data_format.md`` - a human-readable reference table per structure.

Both are derived from the single source of truth (the structs), so they stay in
sync with the format automatically. Regenerate with::

    python generate_transcript_data/generate_json_docs.py
"""
import ast
import json
import os

import msgspec
import msgspec.inspect as mi

from cdot import models, __version__

# Document the structures top-down
STRUCTS = [models.CdotData, models.Transcript, models.GenomeBuild, models.Exon, models.Gene]

# Curated notes for "strangeness" that can't be read off the struct types (issue #77):
# how the object maps are keyed, and source-specific quirks.
FIELD_NOTES = {
    ("CdotData", "transcripts"): "Keyed by transcript accession (e.g. `NM_001637.3`).",
    ("CdotData", "genes"): "Keyed by **gene ID** (e.g. RefSeq `80167`), *not* the symbol - the "
                           "symbol is in each gene's `gene_symbol` field. For UTA-sourced data the "
                           "gene ID is unknown, so a placeholder symbol is used as the key instead.",
}

QUIRKS = """\
## Source-specific notes

* **Gene map keys.** The top-level `genes` map is keyed by gene ID (e.g. RefSeq `80167`),
  not by gene symbol. The human-readable symbol is stored in each record's `gene_symbol` field.
* **UTA-derived data.** When data is built from UTA (a `url` like
  `postgresql://uta.biocommons.org/uta_20210129`), the gene ID is not available - only the symbol.
  In that case a placeholder symbol is used as the `genes` map key.
* **Coordinate systems.** Genomic coordinates (`alt_start`/`alt_end`, build `cds_start`/`cds_end`,
  `start`/`stop`) are 0-based. Transcript (cDNA) coordinates inside each exon
  (`cds_start`/`cds_end`) are 1-based.
"""
MODELS_PY = os.path.join(os.path.dirname(models.__file__), "models.py")
DOCS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs")


def extract_field_docs():
    """{ClassName: {field_name: description}} from the ``\"\"\"...\"\"\"`` after each field."""
    tree = ast.parse(open(MODELS_PY).read())
    out = {}
    for node in tree.body:
        if not isinstance(node, ast.ClassDef):
            continue
        field_docs = {}
        body = node.body
        for i, stmt in enumerate(body):
            if isinstance(stmt, ast.AnnAssign) and isinstance(stmt.target, ast.Name):
                nxt = body[i + 1] if i + 1 < len(body) else None
                if (isinstance(nxt, ast.Expr) and isinstance(nxt.value, ast.Constant)
                        and isinstance(nxt.value.value, str)):
                    desc = " ".join(nxt.value.value.split()).replace("``", "`")
                    field_docs[stmt.target.id] = desc
        out[node.name] = field_docs
    return out


def type_str(t, *, link=True):
    """Render a msgspec.inspect type as a readable string (with Markdown links to structs)."""
    if isinstance(t, mi.Metadata):
        t = t.type
    if isinstance(t, mi.StrType):
        return "string"
    if isinstance(t, mi.IntType):
        return "integer"
    if isinstance(t, mi.FloatType):
        return "number"
    if isinstance(t, mi.BoolType):
        return "boolean"
    if isinstance(t, mi.NoneType):
        return "null"
    if isinstance(t, mi.AnyType):
        return "any"
    if isinstance(t, mi.ListType):
        return f"array of {type_str(t.item_type, link=link)}"
    if isinstance(t, mi.DictType):
        return f"object ({type_str(t.value_type, link=link)} values)"
    if isinstance(t, mi.StructType):
        name = t.cls.__name__
        return f"[{name}](#{name.lower()})" if link else name
    if isinstance(t, mi.UnionType):
        return " or ".join(type_str(x, link=link) for x in t.types)
    return getattr(t, "cls", type(t)).__name__


def field_description(cls_name, field, docs):
    """Field docstring + any curated FIELD_NOTES, combined."""
    parts = [p for p in (docs.get(field, ""), FIELD_NOTES.get((cls_name, field), "")) if p]
    return " ".join(parts)


def first_paragraph(doc):
    doc = (doc or "").strip()
    para = doc.split("\n\n", 1)[0]
    return " ".join(para.split())


def clean_doc(doc):
    """Light RST -> Markdown cleanup for class docstrings."""
    doc = (doc or "").strip().replace("::", ":").replace("``", "`")
    # collapse the indented example block into an inline code fence
    return doc


def build_schema(field_docs):
    schema = msgspec.json.schema(models.CdotData)
    defs = schema.get("$defs", {})
    for cls in STRUCTS:
        name = cls.__name__
        docs = field_docs.get(name, {})
        d = defs.get(name, {})
        props = d.get("properties", {})
        for field in props:
            if desc := field_description(name, field, docs):
                props[field]["description"] = desc
        # Exon is a positional array (prefixItems) - annotate each position by name
        if "prefixItems" in d:
            for idx, field in enumerate(cls.__struct_fields__):
                if idx < len(d["prefixItems"]):
                    item = d["prefixItems"][idx]
                    item["title"] = field
                    if field in docs:
                        item["description"] = docs[field]
    return schema


def render_markdown(field_docs):
    lines = [
        "# cdot JSON data format",
        "",
        "> Auto-generated from the typed models in [`cdot/models.py`](../cdot/models.py) "
        "by `generate_transcript_data/generate_json_docs.py`. Do not edit by hand.",
        "",
        f"Generated from cdot **{__version__}**. "
        "A machine-readable [JSON Schema](cdot-json-schema.json) is generated alongside this file.",
        "",
        "```python",
        "from cdot import models",
        "",
        'data = models.load("cdot-0.2.32.refseq.GRCh38.json.gz")',
        'tx = data.transcripts["NM_001637.3"]',
        "print(tx.gene_name, tx.protein)",
        "for build_name, build in tx.genome_builds.items():",
        "    for exon in build.exons:",
        "        print(exon.alt_start, exon.alt_end, exon.cds_start, exon.gap)",
        "```",
        "",
    ]

    for cls in STRUCTS:
        name = cls.__name__
        info = mi.type_info(cls)
        docs = field_docs.get(name, {})
        lines.append(f"## {name}")
        lines.append("")
        lines.append(clean_doc(cls.__doc__))
        lines.append("")

        if isinstance(info, mi.StructType) and getattr(info, "array_like", False):
            # Positional array (e.g. Exon)
            lines.append("Encoded as a JSON array (positional):")
            lines.append("")
            lines.append("| # | Field | Type | Description |")
            lines.append("|---|-------|------|-------------|")
            for idx, f in enumerate(info.fields):
                lines.append(f"| {idx} | `{f.name}` | {type_str(f.type)} | {docs.get(f.name, '')} |")
        else:
            lines.append("| Field | Type | Required | Description |")
            lines.append("|-------|------|----------|-------------|")
            for f in info.fields:
                req = "yes" if f.required else "no"
                desc = field_description(name, f.name, docs)
                lines.append(f"| `{f.name}` | {type_str(f.type)} | {req} | {desc} |")
        lines.append("")
    lines.append(QUIRKS)
    return "\n".join(lines)


def main():
    os.makedirs(DOCS_DIR, exist_ok=True)
    field_docs = extract_field_docs()

    schema = build_schema(field_docs)
    schema_path = os.path.join(DOCS_DIR, "cdot-json-schema.json")
    with open(schema_path, "w") as f:
        json.dump(schema, f, indent=2)
        f.write("\n")

    md = render_markdown(field_docs)
    md_path = os.path.join(DOCS_DIR, "json_data_format.md")
    with open(md_path, "w") as f:
        f.write(md)

    print(f"wrote {schema_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()
