"""Generate figures for the cdot talk. Palette = dataviz reference (light mode)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import os

HERE = os.path.dirname(os.path.abspath(__file__))
# Generated figures are artefacts, so they go to the repo's ignored output/ directory.
OUT = os.path.normpath(os.path.join(HERE, "..", "..", "output", "talk_figs"))
os.makedirs(OUT, exist_ok=True)

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8880"
BLUE, ORANGE, AQUA, YELLOW, MAGENTA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "axes.edgecolor": "#d9d7d0",
    "axes.labelcolor": INK2,
    "text.color": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": 13,
})


def finish(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=200, facecolor=SURFACE,
                bbox_inches="tight", pad_inches=0.18)
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- 1. coverage
def fig_coverage():
    fig, ax = plt.subplots(figsize=(11.5, 3.5))
    segs = [
        ("RefSeq GRCh37", 197_847, BLUE),
        ("RefSeq GRCh38", 518_810, ORANGE),
        ("Ensembl GRCh37", 196_501, AQUA),
        ("Ensembl GRCh38", 854_940, YELLOW),
        ("T2T-CHM13v2.0", 186_644, MAGENTA),
    ]
    left = 0
    for label, val, colour in segs:
        ax.barh(1, val, left=left, height=0.52, color=colour, label=label,
                edgecolor=SURFACE, linewidth=2)
        left += val
    ax.barh(0, 141_000, height=0.52, color=MUTED, edgecolor=SURFACE, linewidth=2)

    ax.text(1_954_742 + 30_000, 1, "1,954,742", va="center", ha="left",
            fontsize=15, fontweight="bold", color=INK)
    ax.text(141_000 + 30_000, 0, "~141,000  (RefSeq only)", va="center", ha="left",
            fontsize=15, color=INK2)

    ax.set_yticks([0, 1])
    ax.set_yticklabels(["UTA", "cdot"], fontsize=15, color=INK)
    ax.set_ylim(-0.55, 1.55)
    ax.set_xlim(0, 2_400_000)
    ax.set_xlabel("versioned transcript/genome alignments")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{v/1e6:g}M" if v else "0"))
    ax.grid(axis="x", color="#eceae4", linewidth=1)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=3,
              frameon=False, fontsize=11.5, handlelength=1.1, handleheight=1.1)
    finish(fig, "coverage.png")


# ------------------------------------------------- 2. ClinVar resolution split
def fig_clinvar():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 3.7),
                                   gridspec_kw={"width_ratios": [1, 1.5], "wspace": 0.32})

    # left: overall sample
    for i, (lbl, val, colour) in enumerate([("cdot", 99.0, BLUE), ("UTA", 70.3, ORANGE)]):
        ax1.bar(i, val, width=0.55, color=colour, edgecolor=SURFACE, linewidth=2)
        ax1.text(i, val + 2, f"{val}%", ha="center", fontsize=15, fontweight="bold", color=INK)
    ax1.set_xticks([0, 1]); ax1.set_xticklabels(["cdot", "UTA"], fontsize=13, color=INK)
    ax1.set_ylim(0, 112); ax1.set_yticks([0, 25, 50, 75, 100])
    ax1.set_ylabel("resolved (%)")
    ax1.set_title("Seeded ClinVar sample\n(n = 2,818)", fontsize=13, color=INK2, pad=12)
    ax1.grid(axis="y", color="#eceae4", linewidth=1); ax1.set_axisbelow(True)
    ax1.tick_params(axis="x", length=0)

    # right: per source
    groups = ["RefSeq (NM_)\nn = 500", "Ensembl (ENST)\nn = 818"]
    cdot_v, uta_v = [100.0, 99.8], [100.0, 0.0]
    xs = range(len(groups))
    ax2.bar([x - 0.19 for x in xs], cdot_v, width=0.34, color=BLUE, label="cdot",
            edgecolor=SURFACE, linewidth=2)
    ax2.bar([x + 0.19 for x in xs], uta_v, width=0.34, color=ORANGE, label="UTA",
            edgecolor=SURFACE, linewidth=2)
    for x, v in zip(xs, cdot_v):
        ax2.text(x - 0.19, v + 2, f"{v}%", ha="center", fontsize=13, fontweight="bold", color=INK)
    for x, v in zip(xs, uta_v):
        ax2.text(x + 0.19, v + 2, f"{v}%", ha="center", fontsize=13, fontweight="bold", color=INK)
    ax2.set_xticks(list(xs)); ax2.set_xticklabels(groups, fontsize=12, color=INK)
    ax2.set_ylim(0, 112); ax2.set_yticks([0, 25, 50, 75, 100])
    ax2.set_ylabel("resolved (%)")
    ax2.set_title("Committed test pairs, by transcript source", fontsize=13, color=INK2, pad=12)
    ax2.grid(axis="y", color="#eceae4", linewidth=1); ax2.set_axisbelow(True)
    ax2.tick_params(axis="x", length=0)
    ax2.legend(frameon=False, fontsize=12.5, loc="upper center",
               bbox_to_anchor=(0.5, -0.14), ncol=2)
    finish(fig, "clinvar.png")


# -------------------------------------------------------------- 3. throughput
def fig_throughput():
    fig, ax = plt.subplots(figsize=(11.5, 3.5))
    rows = [
        ("UTA: public remote DB", 0.1, ORANGE, "~0.1"),
        ("UTA: local PostgreSQL", 24, ORANGE, "~24"),
        ("cdot REST (per transcript)", 39, BLUE, "~39"),
        ("cdot local JSON", 665, BLUE, "540–665"),
        ("cdot REST (after prefetch)", 731, BLUE, "~731"),
    ]
    ys = range(len(rows))
    for y, (label, val, colour, txt) in zip(ys, rows):
        ax.barh(y, val, height=0.55, color=colour, edgecolor=SURFACE, linewidth=2)
        ax.text(val * 1.25, y, txt, va="center", ha="left", fontsize=14,
                fontweight="bold", color=INK)
    ax.set_yticks(list(ys))
    ax.set_yticklabels([r[0] for r in rows], fontsize=13, color=INK)
    ax.set_xscale("log")
    ax.set_xlim(0.05, 4000)
    ax.set_xticks([0.1, 1, 10, 100, 1000])
    ax.set_xticklabels(["0.1", "1", "10", "100", "1,000"])
    ax.set_xlabel("HGVS resolved per second  (log scale)")
    ax.grid(axis="x", color="#eceae4", linewidth=1); ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    handles = [plt.Rectangle((0, 0), 1, 1, color=BLUE), plt.Rectangle((0, 0), 1, 1, color=ORANGE)]
    ax.legend(handles, ["cdot", "UTA"], frameon=False, fontsize=12,
              loc="lower right", bbox_to_anchor=(1.0, -0.02))
    finish(fig, "throughput.png")


# ---------------------------------------------------------------- 4. cleaning
def fig_cleaning():
    fig, ax = plt.subplots(figsize=(11.5, 3.9))
    cats = [
        ("Prefix / kind restoration", 11),
        ("Other (del/dup count, case)", 21),
        ("Genomic ref in parentheses", 57),
        ("Protein-suffix stripping", 130),
        ("Gene/transcript wrapper repair", 209),
        ("Structural-punctuation repair", 310),
        ("Base re-casing", 553),
        ("Whitespace / non-printable", 940),
    ]
    ys = range(len(cats))
    ax.barh(list(ys), [c[1] for c in cats], height=0.62, color=BLUE,
            edgecolor=SURFACE, linewidth=2)
    for y, (lbl, v) in zip(ys, cats):
        ax.text(v + 18, y, f"{v}   ({v/1678*100:.0f}%)", va="center", fontsize=12.5, color=INK)
    ax.set_yticks(list(ys)); ax.set_yticklabels([c[0] for c in cats], fontsize=12.5, color=INK)
    ax.set_xlim(0, 1250)
    ax.set_xlabel("rescued queries of 1,678  (a string may need several fixes)")
    ax.grid(axis="x", color="#eceae4", linewidth=1); ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    finish(fig, "cleaning.png")


# ---------------------------------------------------------------- 5. residual
def fig_residual():
    fig, ax = plt.subplots(figsize=(11.5, 3.9))
    rows = [
        ("Insertion (length only)", 30, AQUA),
        ("Grammar gap (library limit)", 81, MAGENTA),
        ("Non-HGVS input", 81, MAGENTA),
        ("Trailing / concatenated", 85, ORANGE),
        ("Edit syntax", 113, ORANGE),
        ("Bad accession", 167, ORANGE),
        ("No reference sequence", 277, AQUA),
        ("Truncated", 284, AQUA),
    ]
    ys = range(len(rows))
    for y, (lbl, v, colour) in zip(ys, rows):
        ax.barh(y, v, height=0.62, color=colour, edgecolor=SURFACE, linewidth=2)
        ax.text(v + 5, y, f"{v}  ({v/1118*100:.1f}%)", va="center", fontsize=12, color=INK)
    ax.set_yticks(list(ys)); ax.set_yticklabels([r[0] for r in rows], fontsize=12, color=INK)
    ax.set_xlim(0, 400)
    ax.set_xlabel("residual queries after cleaning (n = 1,118)")
    ax.grid(axis="x", color="#eceae4", linewidth=1); ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (AQUA, ORANGE, MAGENTA)]
    ax.legend(handles,
              ["Information never supplied (~53%)",
               "Fixable: frontier for new rules (~33%)",
               "Out of scope for cleaning (~14%)"],
              frameon=False, fontsize=11.5, loc="lower right")
    finish(fig, "residual.png")


# ---------------------------------------------------------- 6. version safety
def fig_version_safety():
    fig, ax = plt.subplots(figsize=(10.4, 3.7))
    groups = ["RefSeq\n19,393 version bumps", "Ensembl\n12,377 version bumps"]
    per_pair = [96.0, 98.0]
    per_base = [96.8, 98.7]
    xs = range(len(groups))
    ax.bar([x - 0.19 for x in xs], per_pair, width=0.34, color=BLUE,
           label="% of version bumps preserving every coding coordinate",
           edgecolor=SURFACE, linewidth=2)
    ax.bar([x + 0.19 for x in xs], per_base, width=0.34, color=ORANGE,
           label="% chance a random coding variant is unmoved",
           edgecolor=SURFACE, linewidth=2)
    for x, v in zip(xs, per_pair):
        ax.text(x - 0.19, v + 1.5, f"{v}%", ha="center", fontsize=13.5, fontweight="bold", color=INK)
    for x, v in zip(xs, per_base):
        ax.text(x + 0.19, v + 1.5, f"{v}%", ha="center", fontsize=13.5, fontweight="bold", color=INK)
    ax.set_xticks(list(xs)); ax.set_xticklabels(groups, fontsize=12.5, color=INK)
    ax.set_ylim(0, 118); ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("coordinate-safe (%)")
    ax.grid(axis="y", color="#eceae4", linewidth=1); ax.set_axisbelow(True)
    ax.tick_params(axis="x", length=0)
    ax.legend(frameon=False, fontsize=11.5, loc="upper center",
              bbox_to_anchor=(0.5, -0.22), ncol=1)
    finish(fig, "version_safety.png")


for f in (fig_coverage, fig_clinvar, fig_throughput, fig_cleaning,
          fig_residual, fig_version_safety):
    f()
