#!/usr/bin/env python3
"""Generate the synthetic smoke-test data set (SPEC §11.1).

Usage: make_testdata.py OUTDIR [--refs DIR] [--seed N] [--coverage X] [--no-reads]
                        [--hifi-coverage X]

Writes to OUTDIR:
  asm.fa              assembly: chr1-chr7 plus unplaced sequences
  asm.agp             matching AGP 2.0
  reads_R1.fq.gz      simulated paired reads from the "true" genome
  reads_R2.fq.gz
  reads_hifi.fq.gz    only with --hifi-coverage > 0 (development; not in §11.1)
  expected.json       the truth the smoke test asserts against

The true genome differs from the assembly by planted errors: 30 homopolymer
runs (9-12 bp) that are 1-2 nt longer, and 5 dinucleotide repeats with one
copy more or fewer. Reads are error-free apart from 0.1 % substitutions.

--refs must contain rdna_library.fasta and plastid_NC_014057.1.fa
(refs/fetch_refs.py). Output is byte-identical for a given seed and Python
version (random.Random; gzip mtime 0).
"""

import argparse
import gzip
import json
import random
import re
from dataclasses import dataclass, field
from pathlib import Path

LINE_WIDTH = 60
READ_LEN = 150
FRAG_MEAN, FRAG_SD = 400, 30
SUB_RATE = 0.001
QUAL = "F" * READ_LEN  # Q37
HIFI_MEAN, HIFI_SD = 15_000, 3_000
TELO_FWD = "TTTAGGG"  # end arm, 5'->3' towards the chromosome end
TELO_REV = "CCCTAAA"  # start arm
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")

# Chromosome lengths (bp) and arm design: capped / absent / wrong_orientation.
CHROMS = {
    "chr1": (310_000, "capped", "capped"),
    "chr2": (280_000, "capped", "capped"),
    "chr3": (260_000, "capped", "absent"),
    "chr4": (240_000, "capped", "capped"),
    "chr5": (300_000, "wrong_orientation", "capped"),
    "chr6": (270_000, "capped", "capped"),
    "chr7": (220_000, "capped", "capped"),
}
# Telomere array (repeat count, offset from the sequence end) per capped or
# wrong-orientation arm; offsets test distance_from_end_bp.
TELO_ARRAYS = {
    ("chr1", "start"): (600, 0), ("chr1", "end"): (750, 0),
    ("chr2", "start"): (420, 1_500), ("chr2", "end"): (900, 0),
    ("chr3", "start"): (500, 0),
    ("chr4", "start"): (650, 0), ("chr4", "end"): (480, 2_200),
    ("chr5", "start"): (550, 0), ("chr5", "end"): (700, 0),
    ("chr6", "start"): (800, 0), ("chr6", "end"): (450, 0),
    ("chr7", "start"): (520, 900), ("chr7", "end"): (610, 0),
}
DUP_SOURCE = ("chr5", 150_000, 50_017)  # chromosome, start (0-based), length
PLASTID_INSERT = ("chr6", 100_000, 10_000, 30_000)  # chrom, pos, cp start, cp end
PLASTID_FRAGMENT = (60_000, 90_011)  # cp coordinates of the unplaced fragment
N_HP, N_STR2 = 30, 5


@dataclass
class Event:
    """A planted feature at a fixed 0-based assembly position.

    info["run_start"] (planted errors) is 1-based, as in expected.json.
    """

    pos: int
    kind: str
    asm: str
    true: str
    info: dict = field(default_factory=dict)

    @property
    def end(self) -> int:
        return self.pos + len(self.asm)


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records, name, chunks = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks).upper()))
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    return records


class Generator:
    def __init__(self, refs: Path, seed: int):
        self.rng = random.Random(seed)
        self.seed = seed
        lib = read_fasta(refs / "rdna_library.fasta")
        self.rdna = {}
        for name, seq in lib:  # longest unit per subclass
            sub = name.split("#", 1)[1].rsplit("/", 1)[1]
            if len(seq) > len(self.rdna.get(sub, "")):
                self.rdna[sub] = seq
        (_, self.plastid), = read_fasta(refs / "plastid_NC_014057.1.fa")
        if len(self.plastid) != 122_169:
            raise SystemExit(f"unexpected plastid length {len(self.plastid)}")

    # --- sequence primitives -------------------------------------------
    def filler(self, n: int) -> str:
        return "".join(self.rng.choices("ACGT", weights=(32, 18, 18, 32), k=n))

    def base_not(self, *excluded: str) -> str:
        return self.rng.choice([b for b in "ACGT" if b not in excluded])

    # --- fixed events -------------------------------------------------
    def fixed_events(self) -> dict[str, list[Event]]:
        ev: dict[str, list[Event]] = {c: [] for c in CHROMS}

        def add(chrom, pos, kind, asm, true=None, **info):
            ev[chrom].append(Event(pos, kind, asm, asm if true is None else true, info))

        # Gap positions are not multiples of 1,000, so the only round AGP cuts
        # are the planted ones (agp_rows) and the round-length unplaced.
        for chrom, pos, length in [
            ("chr1", 90_137, 100), ("chr1", 250_311, 500), ("chr2", 200_053, 100),
            ("chr3", 130_479, 200), ("chr4", 180_221, 100), ("chr5", 100_613, 100),
            ("chr5", 230_089, 100), ("chr6", 200_167, 100), ("chr7", 110_397, 100),
        ]:
            add(chrom, pos, "gap", "N" * length, self.filler(length))
        # N-runs inside a component: 20 N is a gap for nsplit10, 5 N is not.
        add("chr6", 160_000, "n_run", "N" * 20, self.filler(20))
        add("chr6", 170_000, "n_run", "N" * 5, self.filler(5))

        soft = self.filler(10_000)
        add("chr1", 100_000, "softmask", soft.lower(), soft)
        for i, code in enumerate("RYK"):
            add("chr2", 60_000 + 1_000 * i, "iupac", code, self.filler(1))
        add("chr2", 150_000, "interstitial_telomere", TELO_FWD * 300, repeats=300)

        unit5s = self.rdna["5S"] + self.filler(230)
        add("chr1", 200_000, "rdna5s", unit5s * 10, copies=10, unit_len=len(unit5s))
        unit45s = (self.rdna["18S"] + self.filler(250) + self.rdna["5.8S"] + self.filler(220)
                   + self.rdna["25S"] + self.filler(3_000))
        add("chr4", 120_000, "rdna45s", unit45s * 3, copies=3, unit_len=len(unit45s))

        chrom, pos, cp_start, cp_end = PLASTID_INSERT
        add(chrom, pos, "plastid_insert", self.plastid[cp_start:cp_end],
            cp_start=cp_start + 1, cp_end=cp_end)
        chrom, start, length = DUP_SOURCE
        add(chrom, start, "dup_source", self.filler(length))
        return ev

    # --- planted errors -----------------------------------------------
    def planted_errors(self, ev: dict[str, list[Event]]) -> None:
        candidates = []
        for chrom, (length, _, _) in CHROMS.items():
            for pos in range(60_000, length - 60_000, 10_000):
                if all(pos + 100 < e.pos - 12_000 or pos > e.end + 12_000 for e in ev[chrom]):
                    candidates.append((chrom, pos))
        if len(candidates) < N_HP + N_STR2:
            raise SystemExit(f"only {len(candidates)} positions for planted errors")
        chosen = sorted(self.rng.sample(candidates, N_HP + N_STR2),
                        key=lambda c: (list(CHROMS).index(c[0]), c[1]))
        kinds = ["hp"] * N_HP + ["str2"] * N_STR2
        self.rng.shuffle(kinds)
        str_i = 0
        for (chrom, pos), kind in zip(chosen, kinds):
            if kind == "hp":
                b = self.rng.choice("ACGT")
                run, change = self.rng.randint(9, 12), self.rng.choice((1, 2))
                left, right = self.base_not(b), self.base_not(b)
                ev[chrom].append(Event(
                    pos, "hp", left + b * run + right, left + b * (run + change) + right,
                    {"base": b, "run_length": run, "change": change, "run_start": pos + 2}))
            else:
                unit = ["AT", "AG", "AC", "TG", "CT"][str_i]
                copies = 6 + str_i % 4
                change = 1 if str_i < 3 else -1
                left, right = self.base_not(*unit), self.base_not(*unit)
                ev[chrom].append(Event(
                    pos, "str2", left + unit * copies + right,
                    left + unit * (copies + change) + right,
                    {"unit": unit, "copies": copies, "change_copies": change,
                     "run_start": pos + 2}))
                str_i += 1

    # --- chromosome assembly ------------------------------------------
    def build_chrom(self, chrom: str, events: list[Event]) -> tuple[str, str, list, dict]:
        length, start_arm, end_arm = CHROMS[chrom]
        asm, true = [], []
        cur = 0
        arms = {}

        def put(a: str, t: str):
            nonlocal cur
            asm.append(a)
            true.append(t)
            cur += len(a)

        def telomere(arm: str, status: str) -> tuple[str, int]:
            if status == "absent":
                return "", 0
            n, offset = TELO_ARRAYS[(chrom, arm)]
            motif = TELO_REV if (arm == "start") == (status == "capped") else TELO_FWD
            return motif * n, offset

        start_seq, start_off = telomere("start", start_arm)
        end_seq, end_off = telomere("end", end_arm)
        if start_seq:
            f = self.filler(start_off)
            put(f, f)
            arms["start"] = {"status": start_arm, "start": cur + 1,
                             "end": cur + len(start_seq), "distance_from_end_bp": start_off,
                             "repeats": len(start_seq) // 7, "motif": start_seq[:7]}
            put(start_seq, start_seq)
        else:
            arms["start"] = {"status": "absent"}

        placed = []
        for e in sorted(events, key=lambda e: e.pos):
            if e.pos < cur:
                raise SystemExit(f"{chrom}: event {e.kind} at {e.pos} overlaps previous")
            f = self.filler(e.pos - cur)
            put(f, f)
            put(e.asm, e.true)
            placed.append(e)

        tail = length - cur - len(end_seq) - end_off
        if tail < 0:
            raise SystemExit(f"{chrom}: events exceed chromosome length")
        f = self.filler(tail)
        put(f, f)
        if end_seq:
            arms["end"] = {"status": end_arm, "start": cur + 1, "end": cur + len(end_seq),
                           "distance_from_end_bp": end_off, "repeats": len(end_seq) // 7,
                           "motif": end_seq[:7]}
            put(end_seq, end_seq)
            f = self.filler(end_off)
            put(f, f)
        else:
            arms["end"] = {"status": "absent"}
        return "".join(asm), "".join(true), placed, arms

    # --- unplaced -------------------------------------------------------
    def unplaced(self, chroms: dict[str, str]) -> tuple[dict[str, str], dict[str, str], dict]:
        seqs, true, classes = {}, {}, {}
        seqs["unplaced_short"] = self.filler(15)
        classes["unplaced_short"] = "short"

        body = self.filler(7_284)
        seqs["unplaced_trailing_n"] = body + "N" * 37
        true["unplaced_trailing_n"] = body
        classes["unplaced_trailing_n"] = "unique"

        chrom, start, length = DUP_SOURCE
        seqs["unplaced_dup"] = chroms[chrom][start:start + length]
        classes["unplaced_dup"] = "duplicate"

        cp_start, cp_end = PLASTID_FRAGMENT
        seqs["unplaced_plastid"] = self.plastid[cp_start:cp_end]
        classes["unplaced_plastid"] = "unique"

        for i, n in enumerate([1_000, 1_000, 1_000, 2_000, 2_000], 1):
            name = f"unplaced_round_{i}"
            seqs[name] = true[name] = self.filler(n)
            classes[name] = "unique"
        return seqs, true, classes


# --- AGP -----------------------------------------------------------------
def agp_rows(seqs: dict[str, str]) -> list[list]:
    """Components = runs between N-runs of gap events; names chosen to test §8.1."""
    rows, ctg = [], 0
    special = {("chr2", 1): "subseq_round_A", ("chr4", 2): "subseq_round_B",
               ("chr7", 1): "beg_offset"}
    for obj, seq in seqs.items():
        part, comp_i = 0, 0
        # gaps in the AGP: N-runs >= 100 bp (gap events) and terminal N-runs
        bounds = [m.span() for m in re.finditer(r"N{100,}|N+$", seq)]
        pos = 0
        for gs, ge in bounds + [(len(seq), len(seq))]:
            if gs > pos:
                ctg += 1
                comp_i += 1
                part += 1
                n = gs - pos
                kind = special.get((obj, comp_i))
                beg, end, name = 1, n, f"ctg{ctg:04d}"
                if kind == "subseq_round_A":
                    a = 25_001
                    name = f"ptg{ctg:04d}_subseq_{a}:{a + n - 1}"
                elif kind == "subseq_round_B":
                    b = (-(-n // 1000) + 7) * 1000
                    name = f"ptg{ctg:04d}_subseq_{b - n + 1}:{b}"
                elif kind == "beg_offset":
                    beg, end = 3_001, 3_000 + n
                rows.append([obj, pos + 1, gs, part, "W", name, beg, end, "+"])
            if ge > gs:
                part += 1
                rows.append([obj, gs + 1, ge, part, "N", ge - gs, "scaffold", "yes",
                             "proximity_ligation"])
            pos = ge
    return rows


def agp_cut_coordinates(rows: list[list]) -> list[int]:
    """Cut coordinates per SPEC §8.1 round-number diagnostic."""
    cuts = []
    for r in rows:
        if r[4] != "W":
            continue
        name, beg, end = r[5], r[6], r[7]
        if beg > 1:
            cuts.append(beg - 1)
        m = re.search(r"_subseq_(\d+):(\d+)$", name)
        if m:
            a, b = int(m[1]), int(m[2])
            if a > 1:
                cuts.append(a - 1)
            cuts.append(b)
        else:
            cuts.append(end)
    return cuts


# --- M1 truth --------------------------------------------------------------
def m01_truth(seqs: dict[str, str], agp: list[list], chrom_names) -> dict:
    total = sum(len(s) for s in seqs.values())
    gaps = [(n, m.start(), m.end()) for n, s in seqs.items()
            for m in re.finditer(r"[Nn]{10,}", s)]
    unplaced = [n for n in seqs if n not in chrom_names]
    cuts = agp_cut_coordinates(agp)
    lower = sum(sum(c.islower() for c in s) for s in seqs.values())
    terminal = [n for n, s in seqs.items() if s[:1] in "Nn" or s[-1:] in "Nn"]
    return {
        "n_seq": len(seqs),
        "total_bp": total,
        "n_lt20bp": sum(len(s) < 20 for s in seqs.values()),
        "n_lt200bp": sum(len(s) < 200 for s in seqs.values()),
        "n_terminal_n": len(terminal),
        "n_duplicate_names": 0,
        "n_invalid_chars": 0,
        "n_iupac": sum(len(re.findall(r"[RYSWKMBDHVryswkmbdhv]", s)) for s in seqs.values()),
        "n_gaps_ge10": len(gaps),
        "gap_bp_ge10": sum(e - s for _, s, e in gaps),
        "n_gaps_ge100": sum(e - s >= 100 for _, s, e in gaps),
        "n_seq_gt50pct_n": sum(s.upper().count("N") > 0.5 * len(s) for s in seqs.values()),
        "softmask_pct": round(100 * lower / total, 2),
        "agp_consistent": "yes",
        "n_round_lengths": sum(len(seqs[n]) % 1000 == 0 for n in unplaced),
        "round_lengths_expected": len(unplaced) / 1000,
        "n_round_agp_cuts": sum(c % 1000 == 0 for c in cuts),
        "n_agp_cuts": len(cuts),
        "ena_rules": "FAIL",
        "flags_ena_blocking": {"seq_lt_20bp": ["unplaced_short"],
                               "terminal_n": ["unplaced_trailing_n"]},
    }


# --- output ----------------------------------------------------------------
def write_fasta(path: Path, seqs: dict[str, str]) -> None:
    with path.open("w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), LINE_WIDTH):
                fh.write(seq[i:i + LINE_WIDTH] + "\n")


def simulate_reads(rng: random.Random, true: dict[str, str], coverage: float,
                   r1_path: Path, r2_path: Path) -> int:
    seqs = {n: s.upper() for n, s in true.items() if len(s) >= FRAG_MEAN + 4 * FRAG_SD}
    n_pairs = 0
    with (gzip.GzipFile(r1_path, "wb", mtime=0) as g1,
          gzip.GzipFile(r2_path, "wb", mtime=0) as g2):
        for seq in seqs.values():
            pairs = round(coverage * len(seq) / (2 * READ_LEN))
            for _ in range(pairs):
                frag = max(2 * READ_LEN // 2 + 50, round(rng.gauss(FRAG_MEAN, FRAG_SD)))
                start = rng.randrange(0, len(seq) - frag + 1)
                f = seq[start:start + frag]
                if rng.random() < 0.5:
                    f = f.translate(COMP)[::-1]
                r1, r2 = f[:READ_LEN], f[-READ_LEN:].translate(COMP)[::-1]
                n_pairs += 1
                g1.write(f"@r{n_pairs}/1\n{mutate(rng, r1)}\n+\n{QUAL}\n".encode())
                g2.write(f"@r{n_pairs}/2\n{mutate(rng, r2)}\n+\n{QUAL}\n".encode())
    return n_pairs


def simulate_hifi(rng: random.Random, true: dict[str, str], coverage: float,
                  path: Path) -> int:
    """Long accurate reads, ~15 kb, 0.1 % substitutions, both strands."""
    seqs = {n: s.upper() for n, s in true.items() if len(s) >= 2 * HIFI_MEAN}
    n_reads = 0
    with gzip.GzipFile(path, "wb", mtime=0) as gz:
        for seq in seqs.values():
            for _ in range(round(coverage * len(seq) / HIFI_MEAN)):
                length = max(5_000, round(rng.gauss(HIFI_MEAN, HIFI_SD)))
                start = rng.randrange(0, len(seq) - length + 1)
                r = seq[start:start + length]
                if rng.random() < 0.5:
                    r = r.translate(COMP)[::-1]
                n_reads += 1
                gz.write(f"@h{n_reads}\n{mutate(rng, r)}\n+\n{'F' * len(r)}\n".encode())
    return n_reads


def mutate(rng: random.Random, read: str) -> str:
    pos = int(rng.expovariate(SUB_RATE))
    if pos >= len(read):
        return read
    chars = list(read)
    while pos < len(chars):
        chars[pos] = rng.choice([b for b in "ACGT" if b != chars[pos]])
        pos += 1 + int(rng.expovariate(SUB_RATE))
    return "".join(chars)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("outdir", type=Path)
    ap.add_argument("--refs", type=Path, default=Path("/opt/asmqc/refs"))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--coverage", type=float, default=25.0)
    ap.add_argument("--no-reads", action="store_true")
    ap.add_argument("--hifi-coverage", type=float, default=0.0)
    args = ap.parse_args(argv)

    g = Generator(args.refs, args.seed)
    events = g.fixed_events()
    g.planted_errors(events)

    chroms, true, placed, arms = {}, {}, {}, {}
    for chrom in CHROMS:
        chroms[chrom], true[chrom], placed[chrom], arms[chrom] = g.build_chrom(chrom, events[chrom])
    unpl, unpl_true, m07 = g.unplaced(chroms)
    seqs = chroms | unpl
    true |= unpl_true
    agp = agp_rows(seqs)

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_fasta(args.outdir / "asm.fa", seqs)
    with (args.outdir / "asm.agp").open("w") as fh:
        fh.write("##agp-version\t2.0\n")
        for r in agp:
            fh.write("\t".join(map(str, r)) + "\n")

    def events_of(kind):
        return [{"seq": c, "start": e.pos + 1, "end": e.end, **e.info}
                for c in CHROMS for e in placed[c] if e.kind == kind]

    statuses = [a["status"] for c in CHROMS for a in arms[c].values()]
    expected = {
        "seed": args.seed,
        "coverage": None if args.no_reads else args.coverage,
        "sequences": {n: len(s) for n, s in seqs.items()},
        "m01": m01_truth(seqs, agp, CHROMS),
        "m04": {
            "arms": arms,
            "capped_arms": statuses.count("capped"),
            # T2T: both arms capped and no N-run >= 10 bp (SPEC §8.4)
            "t2t_chromosomes": sum(arms[c]["start"]["status"] == arms[c]["end"]["status"]
                                   == "capped" and not re.search("N{10,}", seqs[c].upper())
                                   for c in CHROMS),
            "wrong_orientation_arms": statuses.count("wrong_orientation"),
            "interstitial": events_of("interstitial_telomere"),
            "unplaced_with_telomere": 0,
        },
        "m05": {
            "plastid_scaffolds": ["unplaced_plastid"],
            "mito_scaffolds": [],
            "chrom_plastid_insert": events_of("plastid_insert"),
            "rdna45s": events_of("rdna45s"),
            "rdna5s": events_of("rdna5s"),
            "rdna_only_scaffolds": [],
        },
        "m07": {"classes": m07},
        "m11": {
            "hp": [{"seq": c, "run_start": e.info["run_start"],
                    "base": e.info["base"], "run_length": e.info["run_length"],
                    "change": e.info["change"], "direction": "ins"}
                   for c in CHROMS for e in placed[c] if e.kind == "hp"],
            "str2": [{"seq": c, "run_start": e.info["run_start"],
                      "unit": e.info["unit"], "copies": e.info["copies"],
                      "change_copies": e.info["change_copies"],
                      "direction": "ins" if e.info["change_copies"] > 0 else "del"}
                     for c in CHROMS for e in placed[c] if e.kind == "str2"],
        },
    }
    if not args.no_reads:
        rng = random.Random(args.seed + 1)
        expected["read_pairs"] = simulate_reads(
            rng, true, args.coverage,
            args.outdir / "reads_R1.fq.gz", args.outdir / "reads_R2.fq.gz")
    if args.hifi_coverage > 0:
        expected["hifi_coverage"] = args.hifi_coverage
        expected["hifi_reads"] = simulate_hifi(random.Random(args.seed + 2), true,
                                               args.hifi_coverage,
                                               args.outdir / "reads_hifi.fq.gz")
    (args.outdir / "expected.json").write_text(json.dumps(expected, indent=1) + "\n")


if __name__ == "__main__":
    main()
