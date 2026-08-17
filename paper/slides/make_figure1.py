"""Render paper/figures/figure1.svg to slide-ready PNGs (whole figure + each band)."""
import os
import cairosvg
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.normpath(os.path.join(HERE, "..", "..", "output", "talk_figs"))
SRC = os.path.normpath(os.path.join(HERE, "..", "figures", "figure1.svg"))
SURFACE = "#fcfcfb"
S = 4  # supersample factor: svg user units -> px

os.makedirs(FIGS, exist_ok=True)

master = os.path.join(FIGS, "figure1_full.png")
cairosvg.svg2png(url=SRC, write_to=master, output_width=1190 * S,
                 output_height=734 * S, background_color=SURFACE)

im = Image.open(master)
# svg content bounds: bands span x 20..1170, y 36..710
crops = {
    # whole figure, outer padding trimmed
    "figure1.png": (14, 28, 1176, 718),
    # band A: data generation (rect y 36..276)
    "figure1_band_a.png": (14, 28, 1176, 279),
    # band B: variant resolution (rect y 306..710)
    "figure1_band_b.png": (14, 298, 1176, 718),
}
for name, (x0, y0, x1, y1) in crops.items():
    out = os.path.join(FIGS, name)
    im.crop((x0 * S, y0 * S, x1 * S, y1 * S)).save(out)
    w, h = Image.open(out).size
    print(f"{name}: {w}x{h}  aspect {w/h:.2f}")
os.remove(master)
