#!/usr/bin/env python3
"""Render Figure S1 (positional drift along the CDS) from positional_drift.csv.

Draws the per-decile preserved-coding-base curves produced by
compute_version_stability.py as a two-panel SVG line chart: panel A conditioned
on the partial-drift version bumps (where any positional effect must live, since
whole-CDS relocation is position independent), panel B the unconditioned curve
over every compared bump. Shared 0-100% y scale so the contrast between the two
panels is visible at a glance.

Dependency-free (writes SVG directly, like the hand-built figure1.svg; the repo
deliberately has no plotting library). Colors are categorical slots 1 and 2 of
the validated reference dataviz palette (blue #2a78d6 / green #008300, a
documented CVD-safe adjacent pair on a light surface); series identity is also
carried by marker shape (circle vs square), direct end labels, and the legend,
never by color alone. Static print figure: light surface only.

Usage:
    python paper/scripts/make_positional_figure.py \
        [--facts paper/empirical_results/positional_drift.csv] \
        [--out paper/figures/figure_s1_positional_drift.svg]
"""
import argparse
import csv
from pathlib import Path

# Chart tokens (reference dataviz palette, light mode, print surface)
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#ffffff"
SERIES = {"refseq": ("#2a78d6", "RefSeq", "circle"),
          "ensembl": ("#008300", "Ensembl", "square")}
FONT = "font-family=\"system-ui, -apple-system, 'Segoe UI', sans-serif\""

N_BINS = 10
# Panel geometry
PLOT_W, PLOT_H = 280, 200
MARGIN_L, MARGIN_T = 64, 64
PANEL_GAP = 96
WIDTH = MARGIN_L + PLOT_W + PANEL_GAP + PLOT_W + 48
HEIGHT = MARGIN_T + PLOT_H + 58


def read_curves(path):
    with open(path, newline="") as fh:
        row = next(csv.DictReader(fh))
    curves = {}
    for label in SERIES:
        for cond in ("partial", "all"):
            curves[(label, cond)] = [
                float(row[f"{label}_{cond}_decile{i}_pct"]) for i in range(1, N_BINS + 1)]
    n_pairs = {label: int(row[f"{label}_partial_pairs"]) for label in SERIES}
    return curves, n_pairs


def xpos(x0, i):
    return x0 + (i + 0.5) * PLOT_W / N_BINS


def ypos(pct):
    return MARGIN_T + PLOT_H * (1 - pct / 100.0)


def marker(shape, cx, cy, color):
    ring = f'stroke="{SURFACE}" stroke-width="2"'
    if shape == "square":
        return (f'<rect x="{cx - 4:.1f}" y="{cy - 4:.1f}" width="8" height="8" '
                f'fill="{color}" {ring}/>')
    return f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="4" fill="{color}" {ring}/>'


def panel(x0, title, subtitle, curves, cond, end_labels):
    parts = [f'<text x="{x0}" y="{MARGIN_T - 38}" fill="{INK}" font-size="13" '
             f'font-weight="600">{title}</text>',
             f'<text x="{x0}" y="{MARGIN_T - 22}" fill="{MUTED}" '
             f'font-size="11">{subtitle}</text>']
    # y gridlines + ticks
    for pct in (0, 25, 50, 75, 100):
        y = ypos(pct)
        parts.append(f'<line x1="{x0}" y1="{y:.1f}" x2="{x0 + PLOT_W}" y2="{y:.1f}" '
                     f'stroke="{GRID}" stroke-width="1"/>')
        parts.append(f'<text x="{x0 - 8}" y="{y + 3.5:.1f}" fill="{MUTED}" '
                     f'font-size="11" text-anchor="end">{pct}</text>')
    # x axis (baseline) + tick labels at each decile
    y_base = ypos(0)
    parts.append(f'<line x1="{x0}" y1="{y_base:.1f}" x2="{x0 + PLOT_W}" '
                 f'y2="{y_base:.1f}" stroke="{AXIS}" stroke-width="1"/>')
    for i in range(N_BINS):
        parts.append(f'<text x="{xpos(x0, i):.1f}" y="{y_base + 16:.1f}" '
                     f'fill="{MUTED}" font-size="11" text-anchor="middle">{i + 1}</text>')
    parts.append(f'<text x="{x0 + PLOT_W / 2:.1f}" y="{y_base + 34:.1f}" '
                 f"fill=\"{INK2}\" font-size=\"11.5\" text-anchor=\"middle\">"
                 "CDS position decile (5&#8242; &#8594; 3&#8242;)</text>")
    # series lines, markers, optional end labels
    for label, (color, name, shape) in SERIES.items():
        vals = curves[(label, cond)]
        pts = [(xpos(x0, i), ypos(v)) for i, v in enumerate(vals)]
        d = " ".join(f"{'M' if i == 0 else 'L'}{px:.1f},{py:.1f}"
                     for i, (px, py) in enumerate(pts))
        parts.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="2" '
                     'stroke-linejoin="round" stroke-linecap="round"/>')
        for px, py in pts:
            parts.append(marker(shape, px, py, color))
        if end_labels:
            px, py = pts[-1]
            parts.append(f'<text x="{px + 9:.1f}" y="{py + 4:.1f}" fill="{INK2}" '
                         f'font-size="11.5" font-weight="600">{name} {vals[-1]:.0f}%</text>')
    return "\n".join(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--facts", default="paper/empirical_results/positional_drift.csv")
    ap.add_argument("--out", default="paper/figures/figure_s1_positional_drift.svg")
    args = ap.parse_args()

    curves, n_pairs = read_curves(args.facts)
    x_a = MARGIN_L
    x_b = MARGIN_L + PLOT_W + PANEL_GAP

    # Legend inside panel A's empty lower-left corner (stacked, clear of the curves)
    legend = []
    lx, ly = x_a + 14, ypos(18)
    for label, (color, name, shape) in SERIES.items():
        legend.append(f'<line x1="{lx}" y1="{ly:.1f}" x2="{lx + 18}" y2="{ly:.1f}" '
                      f'stroke="{color}" stroke-width="2"/>')
        legend.append(marker(shape, lx + 9, ly, color))
        legend.append(f'<text x="{lx + 26}" y="{ly + 4:.1f}" fill="{INK2}" '
                      f'font-size="12">{name}</text>')
        ly += 20

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {WIDTH} {HEIGHT}"
     {FONT} font-size="12">
<rect width="{WIDTH}" height="{HEIGHT}" fill="{SURFACE}"/>
{''.join(legend)}
<text x="14" y="{MARGIN_T + PLOT_H / 2:.0f}" fill="{INK2}" font-size="11.5"
      text-anchor="middle" transform="rotate(-90 14 {MARGIN_T + PLOT_H / 2:.0f})">Coding bases preserved (%)</text>
{panel(x_a, "A&#160;&#160;Partial-drift bumps only",
       f"n = {n_pairs['refseq']} RefSeq / {n_pairs['ensembl']} Ensembl pairs",
       curves, "partial", end_labels=True)}
{panel(x_b, "B&#160;&#160;All version bumps",
       "unconditioned; whole-CDS drift dominates", curves, "all", end_labels=False)}
</svg>
"""
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(svg)
    print(f"Written: {out}")


if __name__ == "__main__":
    main()
