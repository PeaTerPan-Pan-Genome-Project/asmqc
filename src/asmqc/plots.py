"""Report plots (SPEC §7.4) as PNG bytes for embedding.

Colours follow one fixed scheme: categorical slot 1 (blue) for a single
series, slots 1 and 2 (blue, orange) for A/T vs G/C; text in ink colours;
thin marks, a surface gap between stacked segments, one y-axis per chart.
"""

import io

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
SERIES = ("#2a78d6", "#eb6834")  # slot 1 blue, slot 2 orange


def _figure(width=6.4, height=3.2):
    fig, ax = plt.subplots(figsize=(width, height))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK2, labelcolor=INK2, labelsize=8)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    return fig, ax


def _png(fig) -> bytes:
    buf = io.BytesIO()
    fig.tight_layout()
    fig.savefig(buf, format="png", dpi=130, facecolor=SURFACE)
    plt.close(fig)
    return buf.getvalue()


def _bp(v: float) -> str:
    for unit, div in (("Gb", 1e9), ("Mb", 1e6), ("kb", 1e3)):
        if abs(v) >= div:
            return f"{v / div:.3g} {unit}"
    return f"{v:.0f} bp"


def cumulative_contigs(lengths: list[int], method: str) -> bytes:
    """Cumulative contig length against contig rank (longest first)."""
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    cum, acc = [], 0
    for n in lengths:
        acc += n
        cum.append(acc)
    fig, ax = _figure()
    ax.step(range(1, len(cum) + 1), cum, where="post", color=SERIES[0], linewidth=2)
    ax.axhline(total / 2, color=INK2, linewidth=0.8, linestyle=(0, (3, 3)))
    ax.text(len(cum), total / 2, " 50 %", va="bottom", ha="right", color=INK2, fontsize=8)
    ax.set_xlabel(f"contig rank ({method})", color=INK2, fontsize=9)
    ax.set_ylabel("cumulative length", color=INK2, fontsize=9)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: _bp(v)))
    ax.set_xscale("log" if len(cum) > 50 else "linear")
    return _png(fig)


def class_bars(classes: list[str], bp: list[int], title: str) -> bytes:
    """Horizontal bars, one series, value labels at the bar ends."""
    fig, ax = _figure(height=0.45 * len(classes) + 1.0)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", color=GRID, linewidth=0.6)
    y = list(range(len(classes)))[::-1]
    ax.barh(y, bp, height=0.6, color=SERIES[0])
    top = max(bp) if bp and max(bp) > 0 else 1
    for yi, v in zip(y, bp, strict=True):
        ax.text(v + top * 0.01, yi, _bp(v), va="center", color=INK, fontsize=8)
    ax.set_yticks(y, classes, color=INK, fontsize=8)
    ax.set_xlim(0, top * 1.18)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: _bp(v)))
    ax.set_title(title, color=INK, fontsize=9, loc="left")
    return _png(fig)


HP_BINS = ["4-8", "9-10", "11-12", "13-15", "16-20", ">20"]


def hp_errors(rows: list[dict], title: str) -> bytes:
    """HP errors by assembly run-length bin, stacked A/T vs G/C (SPEC §8.11)."""
    at = {b: 0 for b in HP_BINS}
    gc = {b: 0 for b in HP_BINS}
    for r in rows:
        if r["class"] == "hp":
            (at if r["base"] == "AT" else gc)[r["run_length_or_unit"]] += int(r["count"])
    fig, ax = _figure()
    x = range(len(HP_BINS))
    a = [at[b] for b in HP_BINS]
    g = [gc[b] for b in HP_BINS]
    ax.bar(x, a, width=0.7, color=SERIES[0], edgecolor=SURFACE, linewidth=2, label="A/T")
    ax.bar(x, g, width=0.7, bottom=a, color=SERIES[1], edgecolor=SURFACE, linewidth=2,
           label="G/C")
    ax.set_xticks(list(x), HP_BINS, color=INK2, fontsize=8)
    ax.set_xlabel("homopolymer length in the assembly (bp)", color=INK2, fontsize=9)
    ax.set_ylabel("hom-alt HP errors", color=INK2, fontsize=9)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK)
    ax.set_title(title, color=INK, fontsize=9, loc="left")
    return _png(fig)


def scatter(points: list[tuple[str, float, float]], xlabel: str, ylabel: str) -> bytes:
    """One series; every point labelled with its assembly."""
    fig, ax = _figure(height=3.6)
    ax.grid(axis="x", color=GRID, linewidth=0.6)
    xs = [p[1] for p in points]
    ys = [p[2] for p in points]
    ax.scatter(xs, ys, s=40, color=SERIES[0], edgecolor=SURFACE, linewidth=2, zorder=3)
    for label, x, y in points:
        ax.annotate(label, (x, y), xytext=(5, 4), textcoords="offset points",
                    fontsize=7, color=INK)
    ax.set_xlabel(xlabel, color=INK2, fontsize=9)
    ax.set_ylabel(ylabel, color=INK2, fontsize=9)
    return _png(fig)


CHROMS = [f"chr{i}" for i in range(1, 8)]


def dotplot(pts: list[dict], ref_lengths: dict[str, int], lengths: dict[str, int],
            title: str, small: bool = False) -> bytes:
    """BUSCO-anchored dotplot: reference chr1-chr7 (x) against the assembly
    chr1-chr7 (y), unplaced hits in a band on top; slot 1 same, slot 2
    opposite gene orientation."""
    def offsets(ls):
        off, pos = {}, 0
        for c in CHROMS:
            off[c] = pos
            pos += int(ls.get(c, 0) or 0)
        return off, pos

    ro, rtot = offsets(ref_lengths)
    so, stot = offsets(lengths)
    band = stot * 0.04
    fig, ax = _figure(width=2.6 if small else 6.0, height=2.7 if small else 6.0)
    ax.grid(False)
    groups = {"yes": ([], []), "no": ([], [])}
    for p in pts:
        x = ro[p["ref_chromosome"]] + int(p["ref_pos"])
        if p["seq_id"] in so:
            y = so[p["seq_id"]] + int(p["pos"])
        else:  # spread unplaced hits inside the band
            y = stot + band * (0.25 + 0.5 * (sum(map(ord, p["seq_id"])) % 97) / 97)
        groups[p["same_strand"]][0].append(x)
        groups[p["same_strand"]][1].append(y)
    size = 0.6 if small else 1.6
    for key, colour, label in (("yes", SERIES[0], "same orientation"),
                               ("no", SERIES[1], "opposite orientation")):
        ax.scatter(*groups[key], s=size, color=colour, linewidths=0, rasterized=True,
                   label=label)
    for c in CHROMS[1:]:
        ax.axvline(ro[c], color=GRID, lw=0.6)
        ax.axhline(so[c], color=GRID, lw=0.6)
    ax.axhline(stot, color=INK2, lw=0.6)
    ax.set_xlim(0, rtot)
    ax.set_ylim(0, stot + band)
    fs = 6 if small else 8
    lab = [c.removeprefix("chr") for c in CHROMS]
    ax.set_xticks([ro[c] + ref_lengths[c] / 2 for c in CHROMS], lab, fontsize=fs, color=INK2)
    ax.set_yticks([so[c] + int(lengths.get(c, 0) or 0) / 2 for c in CHROMS] + [stot + band / 2],
                  [*lab, "u"], fontsize=fs, color=INK2)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_color(GRID)
    ax.set_title(title, color=INK, fontsize=7 if small else 9, loc="left")
    if not small:
        ax.set_xlabel("Caméor v2 chromosomes", color=INK2, fontsize=9)
        ax.set_ylabel("assembly chromosomes (u = unplaced)", color=INK2, fontsize=9)
        ax.legend(loc="upper left", fontsize=7, markerscale=6, frameon=False, labelcolor=INK)
    return _png(fig)
