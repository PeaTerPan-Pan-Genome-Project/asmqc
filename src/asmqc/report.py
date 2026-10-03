"""Self-contained HTML reports (SPEC §7.4): one per assembly, and the
combined report written by `asmqc aggregate`.

Everything is read from a result directory (<outdir>/<label>/), so the
combined report works on directories sent back by other groups. The one
exception is M2's contig lengths for the cumulative plot, taken from the
workdir while it still exists.
"""

import base64
import bisect
import csv
import gzip
import json
import re
from pathlib import Path

import jinja2

from asmqc import docs, plots, schema, summary
from asmqc import flags as fl

MODULE_TITLES = {
    "m01": "Integrity, format and ENA rules",
    "m02": "Contiguity and anchoring",
    "m04": "Telomeres",
    "m05": "Organelle and rDNA inventory",
    "m06": "Gene-space completeness (BUSCO)",
    "m07": "Redundancy of unplaced scaffolds",
    "m08": "QV and k-mer completeness (Merqury)",
    "m09": "Read-back structural validation (CRAQ)",
    "m11": "Homopolymer and short-STR errors",
}



# --- loading -------------------------------------------------------------------------
def read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load(res: Path) -> dict:
    header, rows = summary.read_summary(res / "qc_summary.tsv")
    return {
        "dir": res,
        "header": header,
        "row": rows[0],
        "flags": fl.read(res / "flags.tsv"),
        "manifest": json.loads((res / "run_manifest.json").read_text()),
    }


def png_uri(data: bytes) -> str:
    return "data:image/png;base64," + base64.b64encode(data).decode()


def file_uri(path: Path) -> str | None:
    return png_uri(path.read_bytes()) if path.exists() else None


# --- cross-checks (report only) ----------------------------------------------------
def dup_buscos_on_duplicates(res: Path) -> int | None:
    """M6 duplicated BUSCOs with an unplaced copy on an M7 `duplicate` scaffold."""
    busco, red = res / "m06_busco" / "busco_derived.tsv", res / "m07_redundancy" / "redundancy.tsv"
    if not (busco.exists() and red.exists()):
        return None
    dups = {r["scaffold"] for r in read_tsv(red) if r["class"] == "duplicate"}
    return sum(1 for r in read_tsv(busco)
               if r["status"] == "Duplicated" and int(r["n_on_unplaced"]) > 0
               and dups & set(r["sequences"].split(";")))


def asm_only_near_errors(res: Path, window: int = 20) -> tuple[int, int] | None:
    """(asm-only k-mers within `window` bp of an M11 error, all asm-only k-mers)."""
    kmers = res / "m08_merqury" / "asm_only_kmers.bed.gz"
    errors = res / "m11_homopolymer" / "errors.bed.gz"
    if not (kmers.exists() and errors.exists()):
        return None
    pos: dict[str, list[int]] = {}
    with gzip.open(errors, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            pos.setdefault(f[0], []).append(int(f[1]))
    for v in pos.values():
        v.sort()
    near = total = 0
    with gzip.open(kmers, "rt") as fh:
        for line in fh:
            seq, s, e = line.split("\t")[:3]
            s, e = int(s), int(e)
            total += 1
            ps = pos.get(seq, [])
            i = bisect.bisect_left(ps, s - window)
            near += i < len(ps) and ps[i] <= e + window
    return near, total


def synteny_dotplot(res: Path, label: str, small: bool = False) -> bytes | None:
    from asmqc import synteny

    pts = read_tsv(res / "m06_busco" / "synteny_points.tsv")
    tab = read_tsv(res / "m06_busco" / "synteny.tsv")
    anchors = synteny.anchors_path()
    if not (pts and tab and anchors.exists()):
        return None
    ref_lengths, _ = synteny.read_anchors(anchors)
    lengths = {r["chromosome"]: int(r["length"]) for r in tab
               if r["chromosome"] != "unplaced" and r["length"] != "NA"}
    title = label if small else f"{label} against Caméor v2 ({len(pts)} BUSCOs)"
    return plots.dotplot(pts, ref_lengths, lengths, title, small=small)


# --- explicit values ------------------------------------------------------------------
GLANCE = ["m01_total_bp", "m02_anchored_pct", "m02_contig_n50", "m01_n_gaps_ge10",
          "m04_t2t_chromosomes", "m04_capped_arms", "m06_complete_pct", "m06_duplicated_pct",
          "m08_qv", "m08_completeness_pct", "m09_r_aqi", "m09_s_aqi", "m11_hp_errors_per_mb"]


def metric(col: str, raw: str) -> dict:
    """A value with its plain label and description (docs.COLUMNS)."""
    label, _, help_ = docs.COLUMNS.get(col, (col, "", ""))
    return {"col": col, "label": label, "value": docs.fmt(col, raw), "help": help_}


def at_a_glance(row: dict) -> list[dict]:
    out = []
    for c in GLANCE:
        if row.get(c, "NA") == "NA":
            continue
        m = metric(c, row[c])
        if docs.COLUMNS[c][1] == "bp":  # short form on the cards; full value in the tooltip
            m["value"] = docs.human_bp(int(row[c]))
            m["help"] = f"{int(row[c]):,} bp. {m['help']}"
        out.append(m)
    return out


def per_chromosome(res: Path) -> list[dict]:
    """One row per chromosome from M2, M4, M6 and M8 outputs (whatever exists)."""
    from asmqc.validate import CHROMOSOMES

    m2 = {r["chromosome"]: r for r in read_tsv(res / "m02_contiguity" / "per_chromosome.tsv")}
    tel: dict[str, dict] = {}
    for r in read_tsv(res / "m04_telomeres" / "telomeres.tsv"):
        tel.setdefault(r["chromosome"], {})[r["arm"]] = r
    syn = {r["chromosome"]: r for r in read_tsv(res / "m06_busco" / "synteny.tsv")}
    qv = {r["seq_id"]: r for r in read_tsv(res / "m08_merqury" / "per_chromosome_qv.tsv")}
    busco: dict[str, int] = {}
    table = res / "m06_busco" / "full_table.tsv"
    if table.exists():
        from asmqc.m06_busco import read_full_table

        for r in read_full_table(table):
            if r["status"] == "Complete":
                busco[r["sequence"]] = busco.get(r["sequence"], 0) + 1
    rows = []
    for c in CHROMOSOMES:
        if c not in m2 and c not in tel:
            continue
        r2, t, s = m2.get(c, {}), tel.get(c, {}), syn.get(c, {})
        rows.append({
            "chromosome": c,
            "length": docs.human_bp(int(r2["length"])) if r2 else "",
            "contigs": r2.get("contigs", ""),
            "gaps": r2.get("gaps", ""),
            "contig_n50": docs.human_bp(int(r2["contig_n50"])) if r2 else "",
            "start": t.get("start", {}).get("status", "").replace("_", " "),
            "end": t.get("end", {}).get("status", "").replace("_", " "),
            "t2t": t.get("start", {}).get("t2t") or (  # results before 0.1.0rc5
                ("yes" if t.get("start", {}).get("status") == t.get("end", {}).get("status")
                 == "capped" and r2.get("gaps") == "0" else "no") if t and r2 else ""),
            "busco": busco.get(c, ""),
            "cameor": (f"{s['best_ref_chromosome']} ({float(s['frac_on_best']):.0%}, "
                       f"{s['orientation']})" if s and s["best_ref_chromosome"] != "NA"
                       else ""),
            "qv": f"{float(qv[c]['qv']):.1f}" if c in qv else "",
        })
    return rows


def rdna_arrays(res: Path, n: int = 15) -> list[dict]:
    rows = [r for r in read_tsv(res / "m05_organelle_rdna" / "rdna_arrays.tsv")
            if r.get("class", "array") == "array"]
    rows.sort(key=lambda r: -int(r["copies"]))
    for r in rows:
        r["span"] = docs.human_bp(int(r["end"]) - int(r["start"]) + 1)
        r["where"] = f"{r['seq_id']}:{int(r['start']):,}-{int(r['end']):,}"
    return rows[:n]


def rdna_plot(res: Path) -> bytes | None:
    arrays = [r for r in read_tsv(res / "m05_organelle_rdna" / "rdna_arrays.tsv")
              if r.get("class", "array") == "array"]
    lengths = {r["seq_id"]: int(r["length"])
               for r in read_tsv(res / "m01_integrity" / "sequences.tsv")}
    if not arrays or not lengths:
        return None
    return plots.rdna_karyoplot(arrays, lengths)


def cse_rows(res: Path, n: int = 100) -> list[dict]:
    rows = read_tsv(res / "m09_craq" / "craq_derived.tsv")
    for r in rows:
        r["where"] = f"{r['seq_id']}:{int(r['start']):,}-{int(r['end']):,}"
        d = r["distance_to_agp_junction_bp"]
        r["distance"] = f"{int(d):,} bp" if d not in ("NA", "") else "no AGP"
    return rows[:n]


CHANGE_RE = re.compile(r"([+-])(\d+)")


def busco_cds_errors(res: Path, top: int = 20) -> dict | None:
    """M11 homopolymer and dinucleotide-repeat errors inside BUSCO coding exons (M6)."""
    cds_file = res / "m06_busco" / "busco_cds.bed.gz"
    err_file = res / "m11_homopolymer" / "errors.bed.gz"
    if not (cds_file.exists() and err_file.exists()):
        return None
    cds: dict[str, list[tuple[int, int, str]]] = {}
    genes = set()
    with gzip.open(cds_file, "rt") as fh:
        for line in fh:
            seq, s, e, bid = line.split("\t")[:4]
            cds.setdefault(seq, []).append((int(s), int(e), bid))
            genes.add(bid)
    cds_bp = 0
    starts = {}
    for seq, ivs in cds.items():
        ivs.sort()
        starts[seq] = [s for s, _, _ in ivs]
        end = -1
        for s, e, _ in ivs:  # merged length
            cds_bp += max(0, e - max(s, end))
            end = max(end, e)
    hits: dict[str, dict] = {}
    n = {"hp": 0, "str2": 0, "frameshift": 0}
    with gzip.open(err_file, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            seq, pos, cls = f[0], int(f[1]), f[3]
            ivs = cds.get(seq)
            if not ivs:
                continue
            k = bisect.bisect_right(starts[seq], pos) - 1
            if k < 0 or not (ivs[k][0] <= pos < ivs[k][1]):
                continue
            m = CHANGE_RE.match(f[6])
            size = int(m[2]) if m else 0
            shift = size % 3 != 0
            n[cls] = n.get(cls, 0) + 1
            n["frameshift"] += shift
            g = hits.setdefault(ivs[k][2], {"busco_id": ivs[k][2], "seq": seq, "errors": 0,
                                            "frameshift": 0, "examples": []})
            g["errors"] += 1
            g["frameshift"] += shift
            if len(g["examples"]) < 3:
                g["examples"].append(f"{f[6]} at {pos + 1:,} ({f[5]}-bp run)" if cls == "hp"
                                     else f"{f[6]} at {pos + 1:,}")
    desc = {}
    table = res / "m06_busco" / "full_table.tsv"
    if table.exists():
        from asmqc.m06_busco import read_full_table

        desc = {r["busco_id"]: r.get("Description", "") for r in read_full_table(table)}
    worst = sorted(hits.values(), key=lambda g: (-g["errors"], g["busco_id"]))[:top]
    for g in worst:
        g["description"] = desc.get(g["busco_id"], "")
    total = n["hp"] + n["str2"]
    return {"genes": len(genes), "cds_bp": cds_bp, "cds": docs.human_bp(cds_bp),
            "hp": n["hp"], "str2": n["str2"], "frameshift": n["frameshift"],
            "genes_affected": len(hits), "per_mb": total / (cds_bp / 1e6) if cds_bp else 0,
            "worst": worst}


# --- per-assembly report -------------------------------------------------------------
def module_sections(data: dict, contig_lengths: list[int] | None) -> list[dict]:
    row, res = data["row"], data["dir"]
    sections = []
    for m in schema.MODULES:
        status = row[f"{m}_status"]
        reason = data["manifest"].get("modules", {}).get(m, {}).get("reason")
        sec = {"id": m, "title": MODULE_TITLES[m], "status": status, "reason": reason,
               "text": docs.MODULE_TEXT[m], "metrics": [], "images": [], "notes": []}
        if status == "ok":
            sec["metrics"] = [metric(c, row[c]) for c in schema.module_columns(m)]
        sections.append(sec)
        if status != "ok":
            continue
        if m == "m02" and contig_lengths:
            sec["images"].append(("Cumulative contig length", png_uri(
                plots.cumulative_contigs(contig_lengths, row["m02_contig_method"]))))
        elif m == "m04":
            sec["images"].append(("Telomere karyoplot",
                                  file_uri(res / "m04_telomeres" / "karyoplot.png")))
            sec["interstitial"] = read_tsv(res / "m04_telomeres" / "interstitial.tsv")
        elif m == "m05":
            rp = rdna_plot(res)
            if rp:
                caption = ("rDNA arrays on chr1–chr7 (marker area ∝ copies; arrays on unplaced "
                           "sequences are in the table)")
                sec["images"].append((caption, png_uri(rp)))
            sec["rdna"] = rdna_arrays(res)
        elif m == "m06":
            dp = synteny_dotplot(res, row["label"])
            if dp:
                caption = "Synteny with Caméor v2 (Complete single-copy BUSCOs; report only)"
                sec["images"].append((caption, png_uri(dp)))
                sec["synteny"] = read_tsv(res / "m06_busco" / "synteny.tsv")
            n = dup_buscos_on_duplicates(res)
            if n is not None:
                sec["notes"].append(f"Duplicated BUSCOs with an unplaced copy on an M7 "
                                    f"duplicate scaffold: {n} of "
                                    f"{row['m06_dup_any_on_unplaced']}.")
        elif m == "m07":
            classes = ["duplicate", "partial_overlap", "repeat_like", "unique", "short"]
            bp = [int(row[f"m07_{c}_bp"]) for c in classes]
            sec["images"].append(("Unplaced bp per class", png_uri(
                plots.class_bars([c.replace("_", " ") for c in classes], bp,
                                 "unplaced sequence by class"))))
        elif m == "m08":
            if row["m08_reads_independent"] == "no":
                sec["notes"].append("The reads were used to build the assembly: the QV is a "
                                    "self-consistency QV, higher than the true accuracy.")
            if row["m08_low_coverage"] == "yes":
                sec["notes"].append(f"k-mer coverage is {row['m08_kmer_coverage']}×, below "
                                    "20×: completeness is underestimated and the QV less "
                                    "reliable.")
            sec["images"] += [("Spectra-cn", file_uri(res / "m08_merqury" / "spectra-cn.png")),
                              ("Spectra-asm",
                               file_uri(res / "m08_merqury" / "spectra-asm.png"))]
        elif m == "m09":
            sec["cse"] = cse_rows(res)
        elif m == "m11":
            rows = read_tsv(res / "m11_homopolymer" / "errors.tsv")
            sec["images"].append(("HP errors by run length", png_uri(
                plots.hp_errors(rows, "hom-alt homopolymer errors"))))
            sec["exons"] = busco_cds_errors(res)
            cc = asm_only_near_errors(res)
            if cc and cc[1]:
                sec["notes"].append(f"Merqury assembly-only k-mers within 20 bp of an error: "
                                    f"{cc[0]:,} of {cc[1]:,} ({100 * cc[0] / cc[1]:.1f} %). "
                                    "Errors found by both methods are the most certain.")
        sec["images"] = [(c, u) for c, u in sec["images"] if u]
    return sections


def environment() -> jinja2.Environment:
    from asmqc.runner import asmqc_home

    return jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(asmqc_home() / "templates")),
        autoescape=True, trim_blocks=True, lstrip_blocks=True)


def render_assembly(res: Path, work: Path | None = None) -> Path:
    data = load(res)
    lengths_file = work / "m02" / "contig_lengths.txt" if work else None
    contig_lengths = ([int(x) for x in lengths_file.read_text().split()]
                      if lengths_file and lengths_file.exists() else None)
    flags = data["flags"]
    html = environment().get_template("report.html.j2").render(
        row=data["row"], manifest=data["manifest"],
        blocking=[f for f in flags if f.severity == fl.ENA_BLOCKING],
        warnings=[f for f in flags if f.severity == fl.WARNING],
        info=[f for f in flags if f.severity == fl.INFO],
        sections=module_sections(data, contig_lengths),
        glance=at_a_glance(data["row"]), chromosomes=per_chromosome(data["dir"]),
        flag_help=docs.FLAGS, glossary=docs.GLOSSARY, files=docs.FILES,
        ena=metric("ena_rules", data["row"]["ena_rules"]),
        params=json.dumps(data["manifest"].get("parameters", {}), indent=1))
    out = res / "report.html"
    out.write_text(html)
    return out


# --- combined report -------------------------------------------------------------------
def render_combined(results: list[dict], out: Path) -> Path:
    groups = [("identity", [c for c in schema.HEADER if c[:3] not in schema.MODULES])]
    groups += [(m, [c for c in schema.HEADER if c.startswith(f"{m}_")]) for m in schema.MODULES]
    karyoplots, hp_plots, dotplots, points = [], [], [], []
    for r in results:
        label, row, res = r["row"]["label"], r["row"], r["dir"]
        k = file_uri(res / "m04_telomeres" / "karyoplot.png")
        if k:
            karyoplots.append((label, k))
        rows = read_tsv(res / "m11_homopolymer" / "errors.tsv")
        if rows:
            hp_plots.append((label, png_uri(plots.hp_errors(rows, label))))
        dp = synteny_dotplot(res, label, small=True)
        if dp:
            dotplots.append((label, png_uri(dp)))
        x, y = row["m06_internal_stop_pct"], row["m11_hp_errors_per_mb"]
        if x != schema.NA and y != schema.NA:
            points.append((label, float(x), float(y)))
    scatter = png_uri(plots.scatter(points, "M6 internal stop codons (% of BUSCOs)",
                                    "M11 HP errors per Mb")) if points else None
    html = environment().get_template("combined.html.j2").render(
        rows=[r["row"] for r in results], groups=groups, titles=MODULE_TITLES, docs=docs.COLUMNS,
        karyoplots=karyoplots, hp_plots=hp_plots, dotplots=dotplots, scatter=scatter,
        n_points=len(points))
    path = out / "report.html"
    path.write_text(html)
    return path
