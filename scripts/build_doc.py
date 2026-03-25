#!/usr/bin/env python3
"""Build a Word document from the paper markdown files.

Run from the project root:
    .venv/bin/python3 scripts/build_doc.py --output-dir output/2026-03-21_8b1471f
"""

import argparse
import os
import shutil
import subprocess
import zipfile
from datetime import date
from pathlib import Path
from lxml import etree

SECTIONS = [
    "paper/abstract.md",
    "paper/introduction.md",
    "paper/implementation.md",
    "paper/discussion.md",
    "paper/references.md",
    "paper/figures.md",
]

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"


def build_docx(docx: Path):
    docx.parent.mkdir(parents=True, exist_ok=True)
    print(f"Building {docx} ...")
    subprocess.run(
        [
            "pandoc",
            *SECTIONS,
            "--from", "markdown",
            "--to", "docx",
            "--output", str(docx),
            "--resource-path", ".",
            "-V", "geometry:margin=2.5cm",
        ],
        check=True,
    )


def strip_bookmarks(docx: Path):
    """Remove Word bookmark elements that pandoc inserts for heading anchors.
    These cause warnings when opening in Google Docs / older Word versions.
    """
    print(f"Stripping bookmarks from {docx} ...")
    tmp = docx.with_suffix(".tmp.docx")
    shutil.copy(docx, tmp)

    with zipfile.ZipFile(tmp, "r") as zin, \
         zipfile.ZipFile(docx, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == "word/document.xml":
                tree = etree.fromstring(data)
                for tag in (f"{{{W}}}bookmarkStart", f"{{{W}}}bookmarkEnd"):
                    for el in tree.iter(tag):
                        el.getparent().remove(el)
                data = etree.tostring(
                    tree, xml_declaration=True, encoding="UTF-8", standalone=True
                )
            zout.writestr(item, data)

    tmp.unlink()


def main():
    parser = argparse.ArgumentParser(description="Build Word doc from paper markdown files.")
    parser.add_argument(
        "--output-dir", required=True, metavar="DIR",
        help="Analysis output directory (used to resolve figure paths)",
    )
    args = parser.parse_args()

    os.chdir(Path(__file__).parent.parent)  # run from project root

    output_dir = Path(args.output_dir)
    docx = output_dir / f"cdot_{date.today()}.docx"

    build_docx(docx)
    strip_bookmarks(docx)
    print(f"Done: {docx}")


if __name__ == "__main__":
    main()
