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
from pathlib import Path

import jinja2

from asmqc import flags as fl
from asmqc import plots, schema, summary

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

HOW_TO_READ = {
    "m01": "ENA rules fail on any sequence shorter than 20 bp, any leading or trailing N, "
           "duplicate names and non-IUPAC characters; everything else is a warning or "
           "information. Round lengths and round AGP cut coordinates point to breaks placed on "
           "a fixed grid during manual Hi-C curation; the chance rate of a round cut is about "
           "0.1 %. They never affect the ENA result.",
    "m02": "In a chromosome-level assembly the scaffold N50 is about the length of one "
           "chromosome; the contig N50 is the discriminating number. Contigs are the AGP "
           "components when the AGP matches the FASTA, otherwise the pieces between runs of "
           "at least 10 N. A difference between the two contig N50 values means N-runs inside "
           "AGP components.",
    "m04": "An arm is capped when a telomeric band lies within 50 kb of the end and is "
           "dominated by the expected strand (CCCTAAA at the start, TTTAGGG at the end). A "
           "band of the other strand is reported as wrong orientation. Interstitial arrays "
           "can be genuine in pea. Unplaced sequences with a terminal band are chromosome "
           "ends that exist but were not anchored.",
    "m05": "Report only. Nuclear insertions of organelle DNA are real biology; a very large "
           "block on a chromosome can indicate a misjoin. The 45S NORs are expected on chr4 "
           "and chr7. Assembled rDNA copy number measures how much of each array was "
           "captured, not the copy number in the plant.",
    "m06": "In high-quality assemblies the complete percentage saturates near 100, so the "
           "duplicated percentage and the placement of complete and duplicated BUSCOs carry "
           "the signal. A duplicated BUSCO with a copy on an unplaced sequence is a likely "
           "false duplication. The internal stop-codon percentage tracks frameshifting "
           "indels in genes.",
    "m07": "A measurement, never a purge. A scaffold is a duplicate only when one collinear "
           "alignment covers at least 90 % of it at 99 % identity and MAPQ 20; summed "
           "coverage would call repetitive scaffolds duplicates in a genome that is about "
           "85 % repeats.",
    "m08": "QV and completeness compare assembly k-mers (k = 21) with read k-mers. QVs are "
           "comparable only between assemblies with the same reads declaration. Below 20x "
           "k-mer coverage the numbers are less reliable.",
    "m09": "CRAQ classifies clipped read alignments as regional (CRE) or structural (CSE) "
           "errors; R-AQI and S-AQI summarise them per Mb (100 = none). A structural "
           "breakpoint near an AGP gap or junction suggests a scaffolding error; one inside "
           "a contig suggests an assembler error.",
    "m11": "Hom-alt indels in homopolymer runs and dinucleotide repeats, called from accurate "
           "reads in the callable region. An insertion relative to the assembly means the "
           "assembled run is too short, the typical error of ONT-based assemblies. "
           "Heterozygous calls in an inbred line mostly reflect mismapping in repeats and are "
           "never counted as errors.",
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


# --- per-assembly report -------------------------------------------------------------
def module_sections(data: dict, contig_lengths: list[int] | None) -> list[dict]:
    row, res = data["row"], data["dir"]
    sections = []
    for m in schema.MODULES:
        status = row[f"{m}_status"]
        reason = data["manifest"].get("modules", {}).get(m, {}).get("reason")
        sec = {"id": m, "title": MODULE_TITLES[m], "status": status, "reason": reason,
               "text": HOW_TO_READ[m], "metrics": [], "images": [], "notes": []}
        if status == "ok":
            sec["metrics"] = [(c.removeprefix(f"{m}_").replace("_", " "), row[c])
                              for c in schema.module_columns(m)]
        sections.append(sec)
        if status != "ok":
            continue
        if m == "m02" and contig_lengths:
            sec["images"].append(("Cumulative contig length", png_uri(
                plots.cumulative_contigs(contig_lengths, row["m02_contig_method"]))))
        elif m == "m04":
            sec["images"].append(("Telomere karyoplot",
                                  file_uri(res / "m04_telomeres" / "karyoplot.png")))
        elif m == "m06":
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
                                    "self-consistency QV.")
            sec["images"] += [("Spectra-cn", file_uri(res / "m08_merqury" / "spectra-cn.png")),
                              ("Spectra-asm",
                               file_uri(res / "m08_merqury" / "spectra-asm.png"))]
        elif m == "m11":
            rows = read_tsv(res / "m11_homopolymer" / "errors.tsv")
            sec["images"].append(("HP errors by run length", png_uri(
                plots.hp_errors(rows, "hom-alt homopolymer errors"))))
            cc = asm_only_near_errors(res)
            if cc and cc[1]:
                sec["notes"].append(f"Merqury assembly-only k-mers within 20 bp of an error: "
                                    f"{cc[0]} of {cc[1]} ({100 * cc[0] / cc[1]:.1f} %).")
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
        params=json.dumps(data["manifest"].get("parameters", {}), indent=1))
    out = res / "report.html"
    out.write_text(html)
    return out


# --- combined report -------------------------------------------------------------------
def render_combined(results: list[dict], out: Path) -> Path:
    groups = [("identity", [c for c in schema.HEADER if c[:3] not in schema.MODULES])]
    groups += [(m, [c for c in schema.HEADER if c.startswith(f"{m}_")]) for m in schema.MODULES]
    karyoplots, hp_plots, points = [], [], []
    for r in results:
        label, row, res = r["row"]["label"], r["row"], r["dir"]
        k = file_uri(res / "m04_telomeres" / "karyoplot.png")
        if k:
            karyoplots.append((label, k))
        rows = read_tsv(res / "m11_homopolymer" / "errors.tsv")
        if rows:
            hp_plots.append((label, png_uri(plots.hp_errors(rows, label))))
        x, y = row["m06_internal_stop_pct"], row["m11_hp_errors_per_mb"]
        if x != schema.NA and y != schema.NA:
            points.append((label, float(x), float(y)))
    scatter = png_uri(plots.scatter(points, "M6 internal stop codons (% of BUSCOs)",
                                    "M11 HP errors per Mb")) if points else None
    html = environment().get_template("combined.html.j2").render(
        rows=[r["row"] for r in results], groups=groups, titles=MODULE_TITLES,
        karyoplots=karyoplots, hp_plots=hp_plots, scatter=scatter, n_points=len(points))
    path = out / "report.html"
    path.write_text(html)
    return path
