#!/usr/bin/env python3
"""Generate the animated SVGs used in the profile README.

    python3 tools/make_svgs.py        # writes assets/*.svg

Everything is plain SVG + CSS keyframes (no JavaScript), so GitHub renders
and animates it inside an <img> tag. Tweak the constants and re-run.
"""

from __future__ import annotations

import math
import random
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "assets"
MONO = "ui-monospace,SFMono-Regular,Menlo,Consolas,'Liberation Mono',monospace"

CYAN = "#36BCF7"
VIOLET = "#a371f7"
REDUCED_MOTION = (
    "@media (prefers-reduced-motion: reduce){*{animation:none!important}}"
)


def f(x: float) -> str:
    """Compact float formatting for SVG coordinates."""
    s = f"{x:.2f}".rstrip("0").rstrip(".")
    return "0" if s == "-0" else s


# --------------------------------------------------------------------------
# 1. Scrolling DNA double-helix divider
# --------------------------------------------------------------------------
def dna_divider() -> str:
    W, H = 800, 56
    cy, amp = H / 2, 17
    period = 64  # px per full helix turn
    rungs_per_period = 8
    # The base sequence repeats every 4 turns, so the loop is seamless when we
    # slide the helix left by exactly that distance.
    seq = "GATTACACGTTAGCCATGGATCCGTACTAGCA"  # 32 bp, GATTACA included ;)
    repeat = period * len(seq) // rungs_per_period  # 256 px
    span = W + repeat + period
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    base_color = {"A": "#3fb950", "T": "#f85149", "G": "#d29922", "C": "#58a6ff"}

    def y1(x: float) -> float:
        return cy + amp * math.sin(2 * math.pi * x / period)

    def y2(x: float) -> float:
        return cy - amp * math.sin(2 * math.pi * x / period)

    # Split each strand into half-turn segments so the strand "in front"
    # can be drawn brighter/thicker than the one behind it.
    half = period / 2
    front, back = [], []
    x = -period / 4
    seg_i = 0
    while x < span:
        xs = [x + i * half / 16 for i in range(17)]
        p1 = "M" + " L".join(f"{f(v)} {f(y1(v))}" for v in xs)
        p2 = "M" + " L".join(f"{f(v)} {f(y2(v))}" for v in xs)
        # cos(theta) > 0 -> strand 1 is in front for this half-turn
        if seg_i % 2 == 0:
            front.append((p1, CYAN))
            back.append((p2, VIOLET))
        else:
            front.append((p2, VIOLET))
            back.append((p1, CYAN))
        x += half
        seg_i += 1

    rungs = []
    step = period / rungs_per_period
    n = int(span // step) + 1
    for i in range(n):
        rx = i * step + step / 2
        b1 = seq[i % len(seq)]
        b2 = comp[b1]
        a, b = y1(rx), y2(rx)
        mid = (a + b) / 2
        depth = abs(math.sin(2 * math.pi * rx / period))
        op = 0.35 + 0.6 * depth
        rungs.append(
            f'<line x1="{f(rx)}" y1="{f(a)}" x2="{f(rx)}" y2="{f(mid)}" '
            f'stroke="{base_color[b1]}" stroke-opacity="{f(op)}"/>'
            f'<line x1="{f(rx)}" y1="{f(mid)}" x2="{f(rx)}" y2="{f(b)}" '
            f'stroke="{base_color[b2]}" stroke-opacity="{f(op)}"/>'
        )

    back_svg = "".join(
        f'<path d="{d}" stroke="{c}" stroke-opacity=".45" stroke-width="2.2"/>'
        for d, c in back
    )
    front_svg = "".join(
        f'<path d="{d}" stroke="{c}" stroke-width="3.4"/>' for d, c in front
    )
    secs = repeat / 16  # ~16 px/s drift

    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-labelledby="t">
<title id="t">Animated DNA double helix</title>
<style>
.helix{{animation:slide {f(secs)}s linear infinite}}
@keyframes slide{{to{{transform:translateX(-{repeat}px)}}}}
{REDUCED_MOTION}
</style>
<defs>
<linearGradient id="fade" x1="0" x2="1">
<stop offset="0" stop-color="#fff" stop-opacity="0"/>
<stop offset=".12" stop-color="#fff"/>
<stop offset=".88" stop-color="#fff"/>
<stop offset="1" stop-color="#fff" stop-opacity="0"/>
</linearGradient>
<mask id="m"><rect width="{W}" height="{H}" fill="url(#fade)"/></mask>
</defs>
<g mask="url(#m)">
<g class="helix" fill="none" stroke-linecap="round">
<g>{back_svg}</g>
<g stroke-width="2.4">{''.join(rungs)}</g>
<g>{front_svg}</g>
</g>
</g>
</svg>
"""


# --------------------------------------------------------------------------
# 2. "Imputing a genotype matrix" terminal animation
# --------------------------------------------------------------------------
GENO = {
    0: ("0/0", "#79c0ff", "#58a6ff"),
    1: ("0/1", "#d2a8ff", "#bc8cff"),
    2: ("1/1", "#ffa657", "#f0883e"),
}


def make_matrix(seed: int = 721):
    rng = random.Random(seed)
    n_pop, per_pop, n_loci = 2, 4, 18
    geno = []
    freqs = []
    for _ in range(n_loci):
        base = rng.uniform(0.1, 0.9)
        diff = rng.uniform(-0.45, 0.45)  # population structure
        freqs.append((min(0.95, max(0.05, base - diff)), min(0.95, max(0.05, base + diff))))
    for p in range(n_pop):
        for _ in range(per_pop):
            row = []
            for l in range(n_loci):
                q = freqs[l][p]
                row.append(int(rng.random() < q) + int(rng.random() < q))
            geno.append(row)
    cells = [(r, c) for r in range(n_pop * per_pop) for c in range(n_loci)]
    missing = set()
    while len(missing) < 15:
        r, c = rng.choice(cells)
        # avoid more than 2 missing calls per locus so the sweep looks even
        if sum(1 for _, cc in missing if cc == c) < 2:
            missing.add((r, c))
    return geno, sorted(missing), per_pop, n_loci


def imputation_card() -> str:
    geno, missing, per_pop, n_loci = make_matrix()
    n_rows = len(geno)
    W = 860
    x0, pitch_x, cw = 128, 39, 33
    y0, pitch_y, ch = 104, 28, 22
    grid_w = pitch_x * n_loci - (pitch_x - cw)
    grid_bottom = y0 + pitch_y * (n_rows - 1) + ch
    bar_y = grid_bottom + 26
    status_y = bar_y + 26
    legend_y = status_y + 34
    H = legend_y + 26

    T = 10.0  # seconds per loop
    t0, t1 = 0.8, 7.0  # scanner sweep window
    reset = 94.0  # % of loop where imputed calls fade back out

    def pct(t: float) -> float:
        return t / T * 100

    def reveal_time(c: int) -> float:
        cx = c * pitch_x + cw / 2
        return t0 + cx / grid_w * (t1 - t0)

    css = [
        f"text{{font-family:{MONO}}}",
        ".imp,.ring{transform-box:fill-box;transform-origin:center}",
        f".scan{{animation:scan {T}s linear infinite}}",
        f"@keyframes scan{{0%{{transform:translateX({x0}px);opacity:0}}"
        f"{f(pct(t0) - 2)}%{{opacity:0}}{f(pct(t0))}%{{transform:translateX({x0}px);opacity:1}}"
        f"{f(pct(t1))}%{{transform:translateX({x0 + grid_w}px);opacity:1}}"
        f"{f(pct(t1) + 3)}%,100%{{transform:translateX({x0 + grid_w}px);opacity:0}}}}",
        f".bar{{transform-origin:{x0}px 0;animation:bar {T}s linear infinite}}",
        f"@keyframes bar{{0%,{f(pct(t0))}%{{transform:scaleX(0);opacity:1}}"
        f"{f(pct(t1))}%{{transform:scaleX(1)}}{reset}%{{transform:scaleX(1);opacity:1}}"
        f"98%,100%{{transform:scaleX(1);opacity:0}}}}",
    ]

    # one keyframe set per locus column (all cells in a column pop together)
    cols = sorted({c for _, c in missing})
    for c in cols:
        p = pct(reveal_time(c))
        css.append(
            f"@keyframes c{c}{{0%,{f(p - 0.2)}%{{opacity:0;transform:scale(.3)}}"
            f"{f(p + 0.8)}%{{opacity:1;transform:scale(1.3)}}"
            f"{f(p + 2.6)}%,{reset}%{{opacity:1;transform:scale(1)}}"
            f"{reset + 3}%,100%{{opacity:0;transform:scale(1)}}}}"
            f".c{c}{{animation:c{c} {T}s ease-out infinite}}"
            f"@keyframes r{c}{{0%,{f(p)}%{{opacity:0;transform:scale(1)}}"
            f"{f(p + 0.5)}%{{opacity:.9;transform:scale(1)}}"
            f"{f(p + 6)}%,100%{{opacity:0;transform:scale(1.9)}}}}"
            f".r{c}{{animation:r{c} {T}s ease-out infinite}}"
        )

    # "missing calls: N" counter that ticks down as the scanner passes
    events = sorted(reveal_time(c) for _, c in missing)
    breakpoints = [0.0] + [pct(t) for t in events]
    counter = []
    total = len(missing)
    remaining_at = {}
    for i, bp in enumerate(breakpoints):
        remaining_at[bp] = total - i
    uniq = sorted(remaining_at)
    for i, start in enumerate(uniq):
        n_left = remaining_at[start]
        end = uniq[i + 1] if i + 1 < len(uniq) else reset + 3
        kf = f"k{i}"
        if i == 0:
            frames = (
                f"0%,{f(end - 0.01)}%{{opacity:1}}{f(end)}%{{opacity:0}}"
                f"{f(reset + 3)}%{{opacity:0}}{f(reset + 3.01)}%,100%{{opacity:1}}"
            )
        else:
            frames = (
                f"0%,{f(start - 0.01)}%{{opacity:0}}{f(start)}%,{f(end - 0.01)}%{{opacity:1}}"
                f"{f(end)}%,100%{{opacity:0}}"
            )
        css.append(f"@keyframes {kf}{{{frames}}}.{kf}{{animation:{kf} {T}s linear infinite}}")
        counter.append(
            f'<text class="{kf}" x="{x0 + grid_w}" y="{status_y}" text-anchor="end" '
            f'font-size="12.5" fill="{"#3fb950" if n_left == 0 else "#8b949e"}">'
            f"missing calls: {n_left}</text>"
        )

    done = pct(t1) + 1
    css.append(
        f"@keyframes s0{{0%,{f(done - 0.01)}%{{opacity:1}}{f(done)}%,{f(reset + 2.99)}%{{opacity:0}}"
        f"{f(reset + 3)}%,100%{{opacity:1}}}}.s0{{animation:s0 {T}s linear infinite}}"
        f"@keyframes s1{{0%,{f(done - 0.01)}%{{opacity:0}}{f(done)}%,{f(reset + 2.99)}%{{opacity:1}}"
        f"{f(reset + 3)}%,100%{{opacity:0}}}}.s1{{animation:s1 {T}s linear infinite}}"
    )
    css.append(REDUCED_MOTION)

    out = []
    # window chrome
    out.append(
        f'<rect x=".5" y=".5" width="{W - 1}" height="{H - 1}" rx="12" fill="#0d1117"/>'
        f'<path d="M12.5 .5h{W - 25}a12 12 0 0 1 12 12V36H.5V12.5a12 12 0 0 1 12-12z" fill="#161b22"/>'
        f'<path d="M.5 36.5h{W - 1}" stroke="#30363d"/>'
        '<circle cx="22" cy="18.5" r="6" fill="#ff5f57"/>'
        '<circle cx="42" cy="18.5" r="6" fill="#febc2e"/>'
        '<circle cx="62" cy="18.5" r="6" fill="#28c840"/>'
        f'<text x="{W / 2}" y="23" text-anchor="middle" font-size="12.5" fill="#8b949e">'
        "genotypes.vcf — deep-learning imputation</text>"
    )
    out.append(
        f'<text x="24" y="66" font-size="13" fill="#8b949e"><tspan fill="#3fb950">#</tspan> '
        f"{n_rows} samples × {n_loci} SNP loci · 2 populations · "
        f'<tspan fill="{CYAN}">filling in the gaps…</tspan></text>'
    )

    # column indices
    for c in range(n_loci):
        out.append(
            f'<text x="{f(x0 + c * pitch_x + cw / 2)}" y="{y0 - 10}" text-anchor="middle" '
            f'font-size="10" fill="#484f58">{c + 1}</text>'
        )
    out.append(f'<text x="24" y="{y0 - 10}" font-size="10" fill="#484f58">sample \\ locus</text>')

    # sample labels + population bars
    for r in range(n_rows):
        pop = "A" if r < per_pop else "B"
        idx = r % per_pop + 1
        out.append(
            f'<text x="30" y="{y0 + r * pitch_y + ch / 2}" dy=".35em" font-size="12" '
            f'fill="#c9d1d9">pop{pop}_{idx:02d}</text>'
        )
    for p, color in ((0, CYAN), (1, VIOLET)):
        top = y0 + p * per_pop * pitch_y
        h = per_pop * pitch_y - (pitch_y - ch)
        out.append(f'<rect x="18" y="{top}" width="4" height="{h}" rx="2" fill="{color}"/>')

    # the scanner sits under the cells' text but over their backgrounds
    scanner = (
        f'<g class="scan"><rect x="-26" y="{y0 - 4}" width="52" '
        f'height="{grid_bottom - y0 + 8}" fill="url(#beam)"/>'
        f'<rect x="-1" y="{y0 - 6}" width="2" height="{grid_bottom - y0 + 12}" '
        f'fill="{CYAN}" opacity=".9"/></g>'
    )

    miss = set(missing)
    cells_bg, cells_fg, imputed = [], [], []
    for r, row in enumerate(geno):
        for c, g in enumerate(row):
            x = x0 + c * pitch_x
            y = y0 + r * pitch_y
            label, txt, edge = GENO[g]
            if (r, c) in miss:
                cells_bg.append(
                    f'<rect x="{x}" y="{y}" width="{cw}" height="{ch}" rx="4" fill="#161b22" '
                    f'stroke="#30363d" stroke-dasharray="3 2"/>'
                )
                cells_fg.append(
                    f'<text x="{f(x + cw / 2)}" y="{f(y + ch / 2)}" dy=".35em" text-anchor="middle" '
                    f'font-size="11" fill="#484f58">./.</text>'
                )
                imputed.append(
                    f'<rect class="ring r{c}" x="{x}" y="{y}" width="{cw}" height="{ch}" rx="4" '
                    f'fill="none" stroke="{CYAN}" stroke-width="1.5"/>'
                    f'<g class="imp c{c}"><rect x="{x}" y="{y}" width="{cw}" height="{ch}" rx="4" '
                    f'fill="#0d1117"/><rect x="{x}" y="{y}" width="{cw}" height="{ch}" rx="4" '
                    f'fill="{edge}" fill-opacity=".22" stroke="{CYAN}" stroke-width="1.6" '
                    f'filter="url(#glow)"/><text x="{f(x + cw / 2)}" y="{f(y + ch / 2)}" dy=".35em" '
                    f'text-anchor="middle" font-size="11" font-weight="700" fill="{txt}">{label}</text></g>'
                )
            else:
                cells_bg.append(
                    f'<rect x="{x}" y="{y}" width="{cw}" height="{ch}" rx="4" fill="{edge}" '
                    f'fill-opacity=".1" stroke="{edge}" stroke-opacity=".3"/>'
                )
                cells_fg.append(
                    f'<text x="{f(x + cw / 2)}" y="{f(y + ch / 2)}" dy=".35em" text-anchor="middle" '
                    f'font-size="11" fill="{txt}">{label}</text>'
                )
    out += cells_bg + [scanner] + cells_fg + imputed

    # progress bar + status line
    out.append(
        f'<rect x="{x0}" y="{bar_y}" width="{grid_w}" height="6" rx="3" fill="#21262d"/>'
        f'<rect class="bar" x="{x0}" y="{bar_y}" width="{grid_w}" height="6" rx="3" fill="url(#prog)"/>'
        f'<text x="24" y="{bar_y + 6}" font-size="11.5" fill="#8b949e">progress</text>'
        f'<text class="s0" x="{x0}" y="{status_y}" font-size="12.5" fill="#8b949e">'
        f'<tspan fill="{CYAN}">●</tspan> predicting missing genotypes…</text>'
        f'<text class="s1" x="{x0}" y="{status_y}" font-size="12.5" fill="#3fb950">'
        "✓ all missing genotypes imputed</text>"
    )
    out += counter

    # legend
    lx = x0
    legend = [
        ("0/0", "#79c0ff", "#58a6ff", "hom. ref", False),
        ("0/1", "#d2a8ff", "#bc8cff", "het", False),
        ("1/1", "#ffa657", "#f0883e", "hom. alt", False),
        ("./.", "#484f58", None, "missing", False),
        ("0/1", "#d2a8ff", "#bc8cff", "imputed", True),
    ]
    for label, txt, edge, name, is_imp in legend:
        if edge is None:
            out.append(
                f'<rect x="{lx}" y="{legend_y - 14}" width="{cw}" height="20" rx="4" '
                f'fill="#161b22" stroke="#30363d" stroke-dasharray="3 2"/>'
            )
        elif is_imp:
            out.append(
                f'<rect x="{lx}" y="{legend_y - 14}" width="{cw}" height="20" rx="4" fill="{edge}" '
                f'fill-opacity=".22" stroke="{CYAN}" stroke-width="1.6" filter="url(#glow)"/>'
            )
        else:
            out.append(
                f'<rect x="{lx}" y="{legend_y - 14}" width="{cw}" height="20" rx="4" fill="{edge}" '
                f'fill-opacity=".1" stroke="{edge}" stroke-opacity=".3"/>'
            )
        out.append(
            f'<text x="{f(lx + cw / 2)}" y="{legend_y - 4}" dy=".35em" text-anchor="middle" '
            f'font-size="10.5" fill="{txt}">{label}</text>'
            f'<text x="{lx + cw + 8}" y="{legend_y - 4}" dy=".35em" font-size="11.5" '
            f'fill="#8b949e">{name}</text>'
        )
        lx += cw + 8 + len(name) * 7.2 + 26

    out.append(
        f'<rect x=".5" y=".5" width="{W - 1}" height="{H - 1}" rx="12" fill="none" stroke="#30363d"/>'
    )

    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-labelledby="t d">
<title id="t">Imputing missing genotypes</title>
<desc id="d">A genotype matrix of {n_rows} samples by {n_loci} SNP loci. A scanner sweeps across the loci and fills each missing call (./.) with a predicted genotype.</desc>
<style>{''.join(css)}</style>
<defs>
<linearGradient id="beam" x1="0" x2="1">
<stop offset="0" stop-color="{CYAN}" stop-opacity="0"/>
<stop offset=".5" stop-color="{CYAN}" stop-opacity=".22"/>
<stop offset="1" stop-color="{CYAN}" stop-opacity="0"/>
</linearGradient>
<linearGradient id="prog" x1="0" x2="1">
<stop offset="0" stop-color="{CYAN}"/>
<stop offset="1" stop-color="{VIOLET}"/>
</linearGradient>
<filter id="glow" x="-50%" y="-50%" width="200%" height="200%">
<feDropShadow dx="0" dy="0" stdDeviation="2.5" flood-color="{CYAN}" flood-opacity=".7"/>
</filter>
</defs>
{chr(10).join(out)}
</svg>
"""


# --------------------------------------------------------------------------
# 3. Slurm / HPC icon styled to sit next to skillicons.dev icons
# --------------------------------------------------------------------------
def slurm_icon() -> str:
    units = []
    for i in range(2):
        y = 40 + i * 54
        units.append(
            f'<rect x="56" y="{y}" width="144" height="40" rx="9" fill="none" '
            f'stroke="{CYAN}" stroke-width="11"/>'
            f'<circle cx="83" cy="{y + 20}" r="8" fill="#3fb950"/>'
            f'<path d="M108 {y + 20}h68" stroke="#c9d1d9" stroke-width="10" stroke-linecap="round"/>'
        )
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="256" height="256" viewBox="0 0 256 256" role="img" aria-labelledby="t">
<title id="t">Slurm HPC</title>
<rect width="256" height="256" rx="60" fill="#242938"/>
{''.join(units)}
<text x="128" y="212" text-anchor="middle" font-family="-apple-system,'Segoe UI',Helvetica,Arial,sans-serif" font-weight="800" font-size="58" letter-spacing="1" fill="#fff">SLURM</text>
</svg>
"""


def main() -> None:
    OUT.mkdir(exist_ok=True)
    files = {
        "dna-divider.svg": dna_divider(),
        "genotype-imputation.svg": imputation_card(),
        "slurm.svg": slurm_icon(),
    }
    for name, svg in files.items():
        (OUT / name).write_text(svg, encoding="utf-8")
        print(f"wrote assets/{name} ({len(svg.encode()) / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
