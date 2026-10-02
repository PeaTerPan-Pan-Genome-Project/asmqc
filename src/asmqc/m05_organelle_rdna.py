"""Module 5: organelle and rDNA inventory (SPEC §8.5). Report only, no PASS/FAIL.

Usage: python -m asmqc.m05_organelle_rdna --paf PAF --blast TSV --fai FAI
           --outdir DIR --work DIR
"""

import argparse
from pathlib import Path

from asmqc import flags as fl
from asmqc import module
from asmqc.params import PARAMS
from asmqc.validate import CHROMOSOMES

P = PARAMS["m05"]
M = "m05"
# accessions of the organelle references (SPEC §6), as in their FASTA headers
ORGANELLES = {"NC_014057.1": "plastid", "PP555264.1": "mito"}


def merge(intervals: list[tuple[int, int]], link: int = 0) -> list[tuple[int, int]]:
    """Merge half-open intervals that overlap or lie within `link` bp."""
    out: list[list[int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1] + link:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def covered(intervals: list[tuple[int, int]]) -> int:
    return sum(e - s for s, e in merge(intervals))


# --- organelles ------------------------------------------------------------
def organelle_intervals(paf: Path, kinds: dict[str, str]) -> dict[str, dict[str, list]]:
    """query -> organelle kind -> query intervals with identity >= threshold."""
    out: dict[str, dict[str, list]] = {}
    for line in paf.read_text().splitlines():
        f = line.split("\t")
        matches, block = int(f[9]), int(f[10])
        if block == 0 or matches / block < P["organelle_min_identity"]:
            continue
        kind = kinds.get(f[5])
        if kind is None:
            raise ValueError(f"unknown organelle target {f[5]}")
        out.setdefault(f[0], {}).setdefault(kind, []).append((int(f[2]), int(f[3])))
    return out


def organelles(paf: Path, kinds: dict[str, str], lengths: dict[str, int]):
    ivs = organelle_intervals(paf, kinds)
    scaffolds, on_chrom, flags = [], [], []
    for seq, by_kind in ivs.items():
        if seq in CHROMOSOMES:
            continue
        n = lengths[seq]
        cov = {k: covered(by_kind.get(k, [])) for k in ("plastid", "mito")}
        top = max(cov, key=lambda k: cov[k])
        is_org = cov[top] >= P["organelle_scaffold_min_cov"] * n
        scaffolds.append((seq, n, f"{cov['plastid'] / n:.4f}", f"{cov['mito'] / n:.4f}",
                          top if is_org else "none"))
        if is_org:
            flags.append(fl.make(M, "organelle_scaffold", seq, 1, n, top,
                                 f"{100 * cov[top] / n:.1f} % covered by the {top} genome"))
    for c in CHROMOSOMES:
        by_kind = ivs.get(c, {})
        row = [c]
        for k in ("plastid", "mito"):
            m = merge(by_kind.get(k, []))
            row += [sum(e - s for s, e in m), max((e - s for s, e in m), default=0)]
        on_chrom.append(tuple(row))
    return scaffolds, on_chrom, flags


# --- rDNA --------------------------------------------------------------------
def rdna_hits(blast: Path) -> dict[tuple[str, str], list[tuple[int, int]]]:
    """(subject, subclass) -> subject intervals (half-open) of hits covering at
    least rdna_min_query_cov of the library subunit (outfmt "6 std qlen slen")."""
    out: dict[tuple[str, str], list] = {}
    for line in blast.read_text().splitlines():
        f = line.split("\t")
        qstart, qend, qlen = int(f[6]), int(f[7]), int(f[12])
        if (abs(qend - qstart) + 1) < P["rdna_min_query_cov"] * qlen:
            continue
        sub = f[0].split("#", 1)[1].rsplit("/", 1)[1]
        s, e = sorted((int(f[8]), int(f[9])))
        out.setdefault((f[1], sub), []).append((s - 1, e))
    return out


def rdna_arrays(blast: Path, lengths: dict[str, int]) -> list[tuple]:
    """Single-linkage clusters of copies of one family within the link distance;
    class `array` with >= rdna_min_array_copies[family] copies, else `fragment`.
    5S needs more: genomes carry many 3-4-copy clusters of 5S-like sequence."""
    copies = {k: merge(v) for k, v in rdna_hits(blast).items()}
    families = {"45S": ("18S", "5.8S", "25S"), "5S": ("5S",)}
    count_by = {"45S": "18S", "5S": "5S"}
    arrays = []
    for seq, length in lengths.items():
        for fam, subs in families.items():
            units = [iv for s in subs for iv in copies.get((seq, s), [])]
            for s, e in merge(units, P["rdna_array_link_bp"]):
                n = sum(1 for a, b in copies.get((seq, count_by[fam]), []) if a >= s and b <= e)
                cls = "array" if n >= P["rdna_min_array_copies"][fam] else "fragment"
                if seq in CHROMOSOMES:
                    context = "chromosome"
                elif cls == "array" and e - s >= P["rdna_only_min_frac"] * length:
                    context = "rdna_only_scaffold"
                else:
                    context = "unplaced"
                arrays.append((seq, s + 1, e, fam, n, cls, context))
    return arrays


def loci(arrays: list[tuple], fam: str) -> str | None:
    """Chromosomes carrying an array of the family, then 'unplaced' if any."""
    seqs = [a[0] for a in arrays if a[3] == fam and a[5] == "array"]
    names = [c for c in CHROMOSOMES if c in seqs]
    if any(s not in CHROMOSOMES for s in seqs):
        names.append("unplaced")
    return ";".join(names) if names else None


def evaluate(paf, blast, kinds, lengths):
    scaffolds, on_chrom, flags = organelles(paf, kinds, lengths)
    arrays = rdna_arrays(blast, lengths)
    org = {k: [s for s in scaffolds if s[4] == k] for k in ("plastid", "mito")}
    rdna_only = {a[0] for a in arrays if a[6] == "rdna_only_scaffold"}

    def count(fam, cls, field):
        return sum((a[4] if field == "copies" else 1) for a in arrays
                   if a[3] == fam and a[5] == cls)

    values = {
        "m05_plastid_scaffolds_n": len(org["plastid"]),
        "m05_plastid_scaffolds_bp": sum(s[1] for s in org["plastid"]),
        "m05_mito_scaffolds_n": len(org["mito"]),
        "m05_mito_scaffolds_bp": sum(s[1] for s in org["mito"]),
        "m05_chrom_plastid_like_bp": sum(r[1] for r in on_chrom),
        "m05_chrom_mito_like_bp": sum(r[3] for r in on_chrom),
        "m05_rdna45s_loci": loci(arrays, "45S"),
        "m05_rdna45s_copies": count("45S", "array", "copies"),
        "m05_rdna45s_fragments": count("45S", "fragment", "n"),
        "m05_rdna5s_loci": loci(arrays, "5S"),
        "m05_rdna5s_copies": count("5S", "array", "copies"),
        "m05_rdna5s_fragments": count("5S", "fragment", "n"),
        "m05_rdna_only_scaffolds_n": len(rdna_only),
        "m05_rdna_only_scaffolds_bp": sum(lengths[s] for s in rdna_only),
    }
    return values, flags, scaffolds, on_chrom, arrays


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", type=Path, required=True)
    ap.add_argument("--blast", type=Path, required=True)
    ap.add_argument("--fai", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    args = ap.parse_args(argv)

    kinds = ORGANELLES
    lengths = module.read_fai(args.fai)
    values, flags, scaffolds, on_chrom, arrays = evaluate(args.paf, args.blast, kinds, lengths)
    out = args.outdir
    module.write_tsv(out / "organelle_scaffolds.tsv",
                     ["seq_id", "length", "plastid_cov_frac", "mito_cov_frac", "organelle"],
                     scaffolds)
    module.write_tsv(out / "organelle_on_chromosomes.tsv",
                     ["chromosome", "plastid_like_bp", "largest_plastid_block_bp",
                      "mito_like_bp", "largest_mito_block_bp"], on_chrom)
    module.write_tsv(out / "rdna_arrays.tsv",
                     ["seq_id", "start", "end", "family", "copies", "class", "context"],
                     arrays)
    module.finish(args.work, M, values, flags)


if __name__ == "__main__":
    main()
