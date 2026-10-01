# Adapted from scripts/telomere_karyoplot.py in
# https://github.com/kavonrtep/ont_genome_assembly_pipeline
# at commit 2144dc8b467b7c606c52b6e8667b8ce24640b077 (GPL-3.0).
# Changes for asmqc (SPEC §8.4): chr1-chr7 only, in chromosome order; bands
# come from asmqc.m04_telomeres (one 50 kb terminal rule, no terminal_frac);
# every band is drawn (terminal and interstitial) and coloured by status, not
# by strand; absent arms are marked; no AGP reorientation.
"""Telomere karyoplot for module 4."""

from pathlib import Path

from asmqc.validate import CHROMOSOMES

# status palette (fixed roles); each colour is paired with a legend label
COLOURS = {
    "capped": "#0ca30c",
    "wrong_orientation": "#d03b3b",
    "interstitial": "#fab219",
    "other": "#8a8a85",  # within 50 kb of an end but not the outermost band
}
LABELS = {
    "capped": "capped arm",
    "wrong_orientation": "wrong orientation",
    "interstitial": "interstitial array",
    "other": "other band",
}
BAR, EDGE, INK, MUTED, SURFACE = "#e9e9e6", "#b4b4ae", "#2b2b28", "#6b6b66", "#fcfcfb"


def draw(path: Path, lengths: dict[str, int], bands_by_seq: dict, arms: list[tuple],
         interstitial: list[tuple], title: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    status = {(a[0], a[1]): a[2] for a in arms}
    inter = {(i[0], i[1]) for i in interstitial}
    max_len = max(lengths.get(c, 0) for c in CHROMOSOMES)
    min_w = max_len * 0.003  # keep short bands visible

    fig, ax = plt.subplots(figsize=(10, 0.45 * len(CHROMOSOMES) + 1.2))
    ax.set_facecolor(SURFACE)
    used = set()
    for row, c in enumerate(CHROMOSOMES):
        y = len(CHROMOSOMES) - row
        n = lengths.get(c, 0)
        ax.add_patch(plt.Rectangle((0, y - 0.28), n, 0.56, facecolor=BAR, edgecolor=EDGE,
                                   linewidth=0.6))
        bs = bands_by_seq.get(c, [])
        for i, b in enumerate(bs):
            if (c, b.start) in inter:
                kind = "interstitial"
            elif i == 0 and status[(c, "start")] != "absent":
                kind = status[(c, "start")]
            elif i == len(bs) - 1 and status[(c, "end")] != "absent":
                kind = status[(c, "end")]
            else:
                kind = "other"
            used.add(kind)
            ax.add_patch(plt.Rectangle((b.start - 1, y - 0.28), max(b.end - b.start + 1, min_w),
                                       0.56, facecolor=COLOURS[kind], edgecolor="none"))
        for arm, x in (("start", 0), ("end", n)):
            if status[(c, arm)] == "absent":
                used.add("absent")
                ax.plot([x], [y + 0.42], marker="x", color=MUTED, markersize=6,
                        markeredgewidth=1.4)

    ax.set_ylim(0.4, len(CHROMOSOMES) + 0.8)
    ax.set_xlim(-max_len * 0.01, max_len * 1.01)
    ax.set_yticks([len(CHROMOSOMES) - r for r in range(len(CHROMOSOMES))])
    ax.set_yticklabels(CHROMOSOMES, fontsize=9, color=INK)
    scale, unit = (1e6, "Mb") if max_len >= 5e6 else (1e3, "kb")
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
        lambda v, _: f"{v / scale:g}"))
    ax.set_xlabel(f"position ({unit})", color=INK)
    ax.tick_params(colors=MUTED, labelcolor=INK)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(EDGE)
    ax.set_title(title, color=INK, fontsize=11, loc="left")

    handles = [Patch(color=COLOURS[k], label=LABELS[k]) for k in COLOURS if k in used]
    if "absent" in used:
        handles.append(Line2D([], [], marker="x", color=MUTED, linestyle="none",
                              markeredgewidth=1.4, label="absent arm"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.18),
              ncol=len(handles), frameon=False, fontsize=8, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
