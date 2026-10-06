"""The CRAQ script replacements in workflow/craq_patch/src give byte-identical
output to the CRAQ 1.10 originals, on seeded random depth tables and candidate
files with the edge cases the originals handle in their own way."""

import os
import random
import shutil
import subprocess
from pathlib import Path

import pytest

from asmqc import tools

PATCH = Path(__file__).parents[2] / "workflow" / "craq_patch" / "src"


def craq_src() -> Path | None:
    exe = (tools.env_bin("craq") / "craq") if tools.env_bin("craq") else shutil.which("craq")
    if not exe or not Path(exe).exists():
        return None
    src = Path(os.path.realpath(exe)).parent.parent / "src"
    return src if (src / "get_ER.pl").exists() else None


SRC = craq_src()
pytestmark = pytest.mark.skipif(SRC is None, reason="CRAQ not installed")


def depth_table(rng: random.Random) -> tuple[str, dict[str, int]]:
    """samtools depth -a style: chromosomes in blocks, positions 1..L ascending.
    Depth profiles: flat, steps, low depth, zero runs of 1-400 bp."""
    chroms = {"chrA": rng.randint(3000, 6000), "chrB": rng.randint(800, 2000),
              "one": 1, "u_7": rng.randint(150, 400), "chrC": rng.randint(2000, 4000)}
    lines = []
    for c, n in chroms.items():
        pos = 1
        while pos <= n:
            seg = rng.randint(1, 400)
            kind = rng.random()
            for p in range(pos, min(n, pos + seg - 1) + 1):
                if kind < 0.2:
                    d = 0
                elif kind < 0.35:
                    d = rng.randint(0, 5)
                elif kind < 0.7:
                    d = rng.randint(20, 40)
                else:
                    d = rng.randint(2, 12)
                lines.append(f"{c}\t{p}\t{d}\n")
            pos += seg
    return "".join(lines), chroms


def sites(rng: random.Random, chroms: dict[str, int], n: int) -> list[tuple[str, int]]:
    """Positions near ends, beyond the end, on missing and one-line sequences."""
    out = []
    names = [*chroms, "absent"]
    for _ in range(n):
        c = rng.choice(names)
        length = chroms.get(c, 500)
        r = rng.random()
        p = (rng.randint(1, 250) if r < 0.15 else rng.randint(max(1, length - 250), length + 30)
             if r < 0.3 else rng.randint(1, max(1, length)))
        out.append((c, p))
    return out


def run(script: Path, args: list[str], cwd: Path) -> bytes:
    r = subprocess.run(["perl", str(script), *args], cwd=cwd, capture_output=True, check=True)
    return r.stdout


def same(name: str, args: list[str], cwd: Path) -> bytes:
    a = run(SRC / name, args, cwd)
    b = run(PATCH / name, args, cwd)
    assert a == b, name
    return a


SEEDS = range(6)


@pytest.mark.parametrize("seed", SEEDS)
def test_get_er(tmp_path, seed):
    rng = random.Random(seed)
    dep, chroms = depth_table(rng)
    (tmp_path / "dep").write_text(dep)
    rows = []
    for c, p in sites(rng, chroms, 400):
        extra = rng.choice(["", "\t0.5\t\t", "\tx"])
        rows.append(f"{c}\t{p}\t{rng.choice('+-')}\t{rng.randint(0, 4)}\t{rng.randint(1, 40)}"
                    f"{extra}\n")
    (tmp_path / "cov").write_text("".join(rows))
    printed = 0
    for w, n, t in [(10, 20, 0.1), (10, 20, 0.5), (10, 20, 0.9), (5, 7, 0.3), (25, 4, 0.8),
                    (3, 30, 0.02), (50, 2, -1)]:
        printed += len(same("get_ER.pl", ["dep", "cov", str(w), str(n), str(t)], tmp_path))
    assert printed > 0  # the print branch is exercised


@pytest.mark.parametrize("seed", SEEDS)
def test_synthesize_lrbkdep(tmp_path, seed):
    rng = random.Random(seed)
    dep, chroms = depth_table(rng)
    (tmp_path / "dep").write_text(dep)
    rows = [f"{c}\t{p}\t{rng.choice(['+', '-', '+', '-', '.'])}\t{rng.randint(2, 30)}"
            f"{rng.choice(['', '\t7', '\t\t'])}\n" for c, p in sites(rng, chroms, 400)]
    (tmp_path / "bk").write_text("".join(rows))
    assert same("synthesize_LRbkdep_and_alldep.pl", ["bk", "dep"], tmp_path)


@pytest.mark.parametrize("seed", SEEDS)
def test_lrcoverrate_srdep_filter(tmp_path, seed):
    rng = random.Random(seed)
    dep, chroms = depth_table(rng)
    (tmp_path / "dep").write_text(dep)
    rows = [f"{c}\t{p}\t{rng.choice('+-')}\t{rng.randint(1, 30)}\t{rng.randint(1, 40)}\n"
            for c, p in sites(rng, chroms, 120)]
    (tmp_path / "cov").write_text("".join(rows))
    assert same("LRcoverRate_srdep_filter.pl", ["dep", "cov"], tmp_path)


@pytest.mark.parametrize("seed", SEEDS)
def test_synthesize_clipdicov(tmp_path, seed):
    rng = random.Random(seed)
    dep, chroms = depth_table(rng)
    (tmp_path / "dep").write_text(dep)
    rows = [f"{c}\t{p}\t{rng.choice('+-')}\t{rng.randint(0, 8)}\t{rng.randint(1, 90)}"
            f"\t{rng.choice(['D', 'I'])}\t{rng.choice([0.3, 0.6, 0.61, 1])}\n"
            for c, p in sites(rng, chroms, 600)]
    (tmp_path / "di").write_text("".join(rows))
    assert same("synthesize_clipDIcov_and_alldep.pl", ["di", "dep"], tmp_path)


@pytest.mark.parametrize("seed", SEEDS)
def test_synthesize_srbkdep(tmp_path, seed):
    rng = random.Random(seed)
    dep, chroms = depth_table(rng)
    (tmp_path / "dep").write_text(dep)
    rows = [f"{c}\t{p}\t{rng.choice(['+', '-', '+', '-', '.'])}\t{rng.randint(1, 30)}\n"
            for c, p in sites(rng, chroms, 600)]
    (tmp_path / "bk").write_text("".join(rows))
    assert same("synthesize_SRbkdep_and_alldep.pl", ["bk", "dep"], tmp_path)


@pytest.mark.parametrize("seed", SEEDS)
def test_search_dep0(tmp_path, seed):
    rng = random.Random(seed)
    dep, _ = depth_table(rng)
    (tmp_path / "dep").write_text(dep)
    assert same("search_dep0.pl", ["dep"], tmp_path)


def test_search_dep0_edges(tmp_path):
    """Runs of exactly 149 and 150 bp, single zeros, two-base runs, names that
    sort differently from file order, and a sequence split into two blocks."""
    def block(c, profile):  # profile: [(depth, length), ...] from position 1
        out, p = [], 1
        for d, n in profile:
            out += [f"{c}\t{q}\t{d}\n" for q in range(p, p + n)]
            p += n
        return out, p

    a, _ = block("z9", [(3, 5), (0, 151), (4, 2), (0, 150), (1, 1), (0, 1), (2, 3), (0, 2)])
    b, _ = block("chr10", [(0, 400), (5, 10), (0, 149)])
    c1, end = block("chr2", [(0, 100)])
    c2 = [f"chr2\t{q}\t0\n" for q in range(end, end + 80)]  # continues in a later block
    (tmp_path / "dep").write_text("".join(a + c1 + b + c2))
    out = same("search_dep0.pl", ["dep"], tmp_path).decode()
    assert "z9\t6\t+\t0" in out and "z9\t159\t+\t0" not in out


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize(("name", "args", "side"), [
    ("LReffect_size.pl", ["../dep", "100", "0.05"], "LRout/Nonmap.loc"),
    ("LReffect_size.pl", ["../dep", "30", "0.2"], "LRout/Nonmap.loc"),
    ("SReffect_size.pl", ["../dep"], "SRout/Nonmap.loc"),
])
def test_effect_size(tmp_path, seed, name, args, side):
    """stdout and the Nonmap.loc each script writes into LRout/ or SRout/."""
    rng = random.Random(seed)
    dep, _ = depth_table(rng)
    (tmp_path / "dep").write_text(dep)
    outs = []
    for src in (SRC, PATCH):
        run_dir = tmp_path / src.parent.name
        (run_dir / side).parent.mkdir(parents=True)
        outs.append((run(src / name, args, run_dir), (run_dir / side).read_bytes()))
    assert outs[0] == outs[1] and outs[0][0] and outs[0][1]


# --- the region path (asmqc_region_depth.pl) against full depth tables --------
SAMTOOLS = (tools.env_bin("craq") / "samtools") if tools.env_bin("craq") else shutil.which("samtools")


def cigar(rng: random.Random, ln: int) -> str:
    """Mostly plain matches; some soft clips at either end, small indels and
    clusters of identical clips (reads ending at the same place)."""
    r = rng.random()
    if r < 0.6:
        return f"{ln}M"
    left = f"{rng.randint(5, 60)}S" if rng.random() < 0.5 else ""
    right = f"{rng.randint(5, 60)}S" if rng.random() < 0.5 else ""
    a = rng.randint(10, ln - 10)
    mid = f"{a}M{rng.randint(3, 50)}{rng.choice('DI')}{ln - a}M" if r < 0.85 else f"{ln}M"
    return left + mid + right


def make_bam(d: Path, rng: random.Random) -> tuple[Path, dict[str, int]]:
    """Sorted, CSI-indexed BAM: covered stretches, read-free stretches, mixed
    MAPQ, a sequence without reads and one with only MAPQ < 20 reads."""
    seqs = {"s1": 30000, "s2": 12000, "noreads": 5000, "lowq": 3000, "s3": 8000}
    lines = ["@HD\tVN:1.6\tSO:unsorted"] + [f"@SQ\tSN:{s}\tLN:{n}" for s, n in seqs.items()]
    k = 0
    for s, n in seqs.items():
        if s == "noreads":
            continue
        for _ in range(n // 40):
            ln = rng.randint(80, 600)
            pos = rng.randint(1, max(1, n - ln))
            if s != "lowq" and (n // 3 < pos < n // 3 + 2000):  # read-free stretch
                continue
            mapq = 5 if s == "lowq" else rng.choice([0, 10, 20, 30, 60, 60])
            k += 1
            lines.append(f"r{k}\t0\t{s}\t{pos}\t{mapq}\t{cigar(rng, ln)}\t*\t0\t0\t*\t*")
    for s, n in seqs.items():  # clusters: indels and clips supported by several reads
        if s in ("noreads", "lowq"):
            continue
        for _ in range(15):
            ln = rng.randint(200, 500)
            pos = rng.randint(1, n - ln)
            cg = cigar(rng, ln)
            for _ in range(rng.randint(2, 6)):
                k += 1
                lines.append(f"r{k}\t0\t{s}\t{pos}\t60\t{cg}\t*\t0\t0\t*\t*")
    sam = d / "x.sam"
    sam.write_text("\n".join(lines) + "\n")
    bam = d / "x.bam"
    subprocess.run([SAMTOOLS, "sort", "-o", bam, sam], check=True, capture_output=True)
    subprocess.run([SAMTOOLS, "index", "-c", bam], check=True, capture_output=True)
    return bam, seqs


def full_table(bam: Path, mapq: int | None) -> bytes:
    if mapq is None:  # long reads: samtools depth -a BAM
        return subprocess.run([SAMTOOLS, "depth", "-a", bam], check=True,
                              capture_output=True).stdout
    view = subprocess.run([SAMTOOLS, "view", "-h", "-q", str(mapq), bam], check=True,
                          capture_output=True).stdout  # short reads: view -q | depth -a -
    return subprocess.run([SAMTOOLS, "depth", "-a", "-"], input=view, check=True,
                          capture_output=True).stdout


@pytest.mark.skipif(SAMTOOLS is None, reason="samtools not installed")
@pytest.mark.parametrize("seed", range(4))
def test_region_depth_path(tmp_path, seed):
    rng = random.Random(seed)
    bam, seqs = make_bam(tmp_path, rng)
    env = os.environ | {"PATH": f"{Path(SAMTOOLS).parent}:{os.environ['PATH']}"}
    # candidate counts keep the windows sparse (region output smaller than the table)
    for script, flank, mapq, args, n in [("get_ER.pl", 200, None, ["10", "20", "0.1"], 300),
                                         ("get_ER.pl", 200, None, ["10", "20", "0.6"], 300),
                                         ("LRcoverRate_srdep_filter.pl", 1500, 20, [], 30)]:
        full = full_table(bam, mapq)
        (tmp_path / "full").write_bytes(full)
        (tmp_path / "seqs").write_text(
            "".join(dict.fromkeys(ln.split(b"\t")[0].decode() + "\n" for ln in full.splitlines())))
        rows = [f"{c}\t{p}\t{rng.choice('+-')}\t{rng.randint(0, 4)}\t{rng.randint(1, 40)}\n"
                for c, p in sites(rng, seqs, n)]
        (tmp_path / "cand").write_text("".join(rows))
        helper = ["perl", str(PATCH / "asmqc_region_depth.pl"), "cand", str(flank), str(bam),
                  "seqs", *([str(mapq)] if mapq is not None else [])]
        region = subprocess.run(helper, cwd=tmp_path, env=env, check=True,
                                capture_output=True).stdout
        (tmp_path / "region").write_bytes(region)
        assert len(region) < len(full)
        for src in (SRC, PATCH):  # original and replacement, each on both inputs
            outs = [run(src / script, ["full", "cand", *args], tmp_path),
                    run(src / script, ["region", "cand", *args], tmp_path)]
            assert outs[0] == outs[1], (script, src)
        assert outs[0], script  # something was printed


def test_fanout(tmp_path):
    data = b"".join(f"c\t{i}\t{i % 7}\n".encode() for i in range(200000))
    ok = subprocess.run(["perl", str(PATCH / "asmqc_fanout.pl"), "cat > a", "wc -l > b",
                         "cut -f1 | uniq > c"], input=data, cwd=tmp_path, capture_output=True,
                        check=False)
    assert ok.returncode == 0
    assert (tmp_path / "a").read_bytes() == data
    assert (tmp_path / "b").read_text().strip() == "200000"
    assert (tmp_path / "c").read_text() == "c\n"
    bad = subprocess.run(["perl", str(PATCH / "asmqc_fanout.pl"), "cat > a", "head -c 10 >/dev/null",
                          "exit 3"], input=data, cwd=tmp_path, capture_output=True,
                         check=False, timeout=60)
    assert bad.returncode != 0 and b"asmqc_fanout" in bad.stderr


@pytest.mark.skipif(SAMTOOLS is None, reason="samtools not installed")
@pytest.mark.parametrize("seed", range(4))
@pytest.mark.parametrize(("flank", "mapq"), [(200, None), (1500, 20), (3, None)])
def test_region_depth_contract(tmp_path, seed, flank, mapq):
    """Output = the full table's lines within FLANK of a candidate, plus
    position 1 of each candidate's sequence, in table order."""
    rng = random.Random(seed)
    bam, seqs = make_bam(tmp_path, rng)
    full = full_table(bam, mapq).decode().splitlines(keepends=True)
    present = list(dict.fromkeys(ln.split("\t")[0] for ln in full))
    (tmp_path / "seqs").write_text("".join(s + "\n" for s in present))
    cands = sites(rng, seqs, 40)
    (tmp_path / "cand").write_text("".join(f"{c}\t{p}\t+\t1\t9\n" for c, p in cands))
    want = {}
    for c, p in cands:
        want.setdefault(c, set()).add(1)
        want[c].update(range(p - flank, p + flank + 1))
    expected = "".join(ln for ln in full
                       if int(ln.split("\t")[1]) in want.get(ln.split("\t")[0], ()))
    env = os.environ | {"PATH": f"{Path(SAMTOOLS).parent}:{os.environ['PATH']}"}
    got = subprocess.run(["perl", str(PATCH / "asmqc_region_depth.pl"), "cand", str(flank),
                          str(bam), "seqs", *([str(mapq)] if mapq is not None else [])],
                         cwd=tmp_path, env=env, check=True, capture_output=True).stdout.decode()
    assert got == expected and expected


# --- the segmented passes (asmqc_segments/par/merge) against one genome-wide pass ---
def sh(cmd: str, cwd: Path) -> bytes:
    env = os.environ | {"PATH": f"{Path(SAMTOOLS).parent}:{os.environ['PATH']}"}
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd], cwd=cwd, env=env, check=True,
                          capture_output=True).stdout


@pytest.mark.skipif(SAMTOOLS is None, reason="samtools not installed")
@pytest.mark.parametrize(("seed", "nseg"), [(0, 1), (1, 2), (2, 3), (3, 8)])
def test_segmented_passes(tmp_path, seed, nseg):
    rng = random.Random(seed)
    bam, seqs = make_bam(tmp_path, rng)
    segs = sh(f"perl {PATCH}/asmqc_segments.pl {bam} {nseg} seg", tmp_path).decode().split()
    assert 1 <= len(segs) <= max(nseg, 1)
    beds = "".join((tmp_path / s).read_text() for s in segs)
    assert [ln.split("\t")[0] for ln in beds.splitlines()] == list(seqs)  # all, in order

    def per_seg(tmpl: str) -> list[bytes]:
        return [sh(tmpl.format(s=s), tmp_path) for s in segs]

    # filtered BAM: same records in the same order
    one = sh(f"samtools view -h -q 20 -F 1796 {bam} | perl {SRC}/lrsam_cigar_filter.pl - "
             "| samtools view", tmp_path)
    parts = per_seg(f"samtools view -h -q 20 -F 1796 -M -L {{s}} {bam} "
                    f"| perl {SRC}/lrsam_cigar_filter.pl - | samtools view -S -b - -o {{s}}.f.bam")
    assert not any(parts)
    sh("samtools cat -o cat.bam " + " ".join(f"{s}.f.bam" for s in segs), tmp_path)
    assert sh("samtools view cat.bam", tmp_path) == one and one
    sh("samtools index -c cat.bam", tmp_path)

    # depth tables: long reads (samtools depth -a BAM) and short reads (view -q | depth)
    for whole, seg in [("samtools depth -a cat.bam", "samtools view -h -M -L {s} cat.bam"),
                       (f"samtools view -h -q 20 {bam} | samtools depth -a -",
                        f"samtools view -h -q 20 -M -L {{s}} {bam}")]:
        assert b"".join(per_seg(seg + " | samtools depth -a -")) == sh(whole, tmp_path)

    # clip and indel scans: same lines (CRAQ prints them in hash order); D before I
    clip = sorted(sh(f"samtools view cat.bam | perl {SRC}/caculate_breakpoint_depth.pl -",
                     tmp_path).splitlines())
    clip_seg = per_seg(f"samtools view -M -L {{s}} cat.bam | perl {SRC}/caculate_breakpoint_depth.pl -")
    assert sorted(b"".join(clip_seg).splitlines()) == clip and clip
    di = sh(f"samtools view cat.bam | perl {SRC}/caculate_clipDI_cov.pl - 3", tmp_path)
    for i, out in enumerate(per_seg(f"samtools view -M -L {{s}} cat.bam "
                                    f"| perl {SRC}/caculate_clipDI_cov.pl - 3")):
        (tmp_path / f"di{i}").write_bytes(out)
    merged = sh(f"perl {PATCH}/asmqc_merge.pl dici " + " ".join(f"di{i}" for i in range(len(segs))),
                tmp_path)
    kinds = [ln.split(b"\t")[4] for ln in merged.splitlines()]
    assert kinds == sorted(kinds, key=lambda k: k == b"I")  # all D, then all I
    assert sorted(merged.splitlines()) == sorted(di.splitlines()) and di

    # table readers per segment, merged, against one pass over the whole table
    full = sh("samtools depth -a cat.bam", tmp_path)
    (tmp_path / "full").write_bytes(full)
    (tmp_path / "LRout").mkdir()
    (tmp_path / "SRout").mkdir()
    (tmp_path / "bk").write_text("".join(f"{c}\t{p}\t{rng.choice('+-')}\t{rng.randint(2, 9)}\n"
                                         for c, p in sites(rng, seqs, 200)))
    whole = {
        "eff": sh(f"perl {SRC}/LReffect_size.pl full 100 0.05", tmp_path),
        "nonmap": (tmp_path / "LRout/Nonmap.loc").read_bytes(),
        "bk": sh(f"perl {SRC}/synthesize_LRbkdep_and_alldep.pl bk full", tmp_path),
        "srbk": sh(f"perl {SRC}/synthesize_SRbkdep_and_alldep.pl bk full", tmp_path),
        "dep0": sh(f"perl {SRC}/search_dep0.pl full", tmp_path),
        "seff": sh(f"perl {SRC}/SReffect_size.pl full", tmp_path),
    }
    for i, s in enumerate(segs):
        sh(f"samtools view -h -M -L {s} cat.bam | samtools depth -a - | perl {PATCH}/asmqc_fanout.pl "
           f"'perl {PATCH}/LReffect_size.pl /dev/stdin 100 0.05 p{i}.nonmap > p{i}.eff' "
           f"'perl {PATCH}/synthesize_LRbkdep_and_alldep.pl bk /dev/stdin > p{i}.bk' "
           f"'perl {PATCH}/synthesize_SRbkdep_and_alldep.pl bk /dev/stdin > p{i}.srbk' "
           f"'perl {PATCH}/search_dep0.pl /dev/stdin > p{i}.dep0' "
           f"'perl {PATCH}/SReffect_size.pl /dev/stdin p{i}.snonmap > p{i}.seff'", tmp_path)
    ps = range(len(segs))
    cat = lambda ext: b"".join((tmp_path / f"p{i}.{ext}").read_bytes() for i in ps)
    merge = lambda mode, ext: sh(f"perl {PATCH}/asmqc_merge.pl {mode} "
                                 + " ".join(f"p{i}.{ext}" for i in ps), tmp_path)
    assert merge("effsize", "eff") == whole["eff"] and whole["eff"]
    assert cat("nonmap") == whole["nonmap"] and whole["nonmap"]
    assert cat("bk") == whole["bk"] and whole["bk"]
    assert cat("srbk") == whole["srbk"] and whole["srbk"]
    assert merge("dep0", "dep0") == whole["dep0"]
    assert merge("effsize", "seff") == whole["seff"]


# --- issue #3: zero-depth list readers -----------------------------------------
def nonmap(rng: random.Random, dep: str, shuffle: bool = False) -> str:
    """Nonmap.loc lines (depth 0) from a depth table, optionally shuffled."""
    lines = [ln + "\n" for ln in dep.splitlines() if ln.endswith("\t0")]
    if shuffle:
        rng.shuffle(lines)
    return "".join(lines)


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize("shuffle", [False, True])
def test_get_nonmap_region(tmp_path, seed, shuffle):
    rng = random.Random(seed)
    dep, _ = depth_table(rng)
    (tmp_path / "sr").write_text(nonmap(rng, dep, shuffle))
    dep2, _ = depth_table(random.Random(seed + 100))
    (tmp_path / "lr").write_text(nonmap(rng, dep2, shuffle))
    out = same("get_nonmap_region.pl", ["sr", "sr"], tmp_path)  # as in "Search noisy error region"
    assert out
    assert same("get_nonmap_region.pl", ["sr", "lr"], tmp_path)  # as in "Create final report"


def test_get_nonmap_region_edges(tmp_path):
    """p and p+4 make one run, p and p+5 two; sequences out of string order;
    a sequence only in the second file; duplicates."""
    a = "".join(f"z\t{p}\t0\n" for p in [10, 14, 30, 35, 35, 100])
    b = "".join(f"chr10\t{p}\t0\n" for p in [1, 2, 3, 50]) + "chr2\t7\t0\n"
    (tmp_path / "f").write_text(a + b)
    (tmp_path / "g").write_text(b + "only2\t5\t0\n" + a)
    out = same("get_nonmap_region.pl", ["f", "f"], tmp_path).decode()
    assert "z\t10\t17\t0\nz\t30\t33\t0\nz\t35\t38\t0\n" in out
    same("get_nonmap_region.pl", ["f", "g"], tmp_path)
    same("get_nonmap_region.pl", ["g", "f"], tmp_path)


@pytest.mark.parametrize("seed", SEEDS)
def test_remove_ngs_normal(tmp_path, seed):
    rng = random.Random(seed)
    dep, chroms = depth_table(rng)
    (tmp_path / "nonmap").write_text(nonmap(rng, dep))
    (tmp_path / "cov").write_text("".join(
        f"{c}\t{p}\t+\t{rng.randint(0, 9)}\t{rng.choice([0, 1, 5, 20, 90])}\n"
        for c, p in sites(rng, chroms, 300)))
    (tmp_path / "cre").write_text("".join(
        f"{c}\t{p}\t+\t{rng.randint(1, 9)}\t{rng.randint(1, 40)}\n"
        for c, p in sites(rng, chroms, 300)))
    out = same("remove_ngs_normal.pl", ["cov", "nonmap", "cre"], tmp_path)
    assert out and len(out.splitlines()) < 300  # some kept, some removed


@pytest.mark.parametrize("seed", SEEDS)
def test_search_uncertain_region(tmp_path, seed):
    """Same lines as the original, which prints them in hash order."""
    rng = random.Random(seed)
    regions = []
    for c in ["chr1", "chr2", "u9"]:
        p = 1
        for _ in range(60):
            p += rng.randint(1, 3000)
            ln = rng.choice([100, 499, 500, 501, 2000, 8000])
            regions.append(f"{c}\t{p}\t{p + ln}\t0\n")
            p += ln
    regions.append(regions[5])  # duplicate
    (tmp_path / "bed").write_text("".join(regions))
    errs = [f"{c}\t{rng.randint(1, 400000)}\t+\t3\t9\n" for c in ["chr1", "chr2", "x"]
            for _ in range(150)]
    s, e = (int(x) for x in regions[3].split("\t")[1:3])
    errs += [f"chr1\t{s - 50}\t+\n", f"chr1\t{e + 51}\t+\n"]  # at the window edges
    (tmp_path / "err").write_text("".join(errs))
    a = run(SRC / "search_uncertain_region.pl", ["bed", "err"], tmp_path)
    b = run(PATCH / "search_uncertain_region.pl", ["bed", "err"], tmp_path)
    assert sorted(a.splitlines()) == sorted(b.splitlines()) and a
    assert len(set(b.splitlines())) == len(b.splitlines())


@pytest.mark.parametrize(("name", "args", "stdin"), [
    ("search_uncertain_region.pl", ["-", "err"], "bed"),  # as in runAQI.sh: piped input
    ("get_nonmap_region.pl", ["nm", "-"], "nm2"),
    ("get_ER.pl", ["-", "cov", "10", "20", "0.1"], "dep"),
    ("remove_ngs_normal.pl", ["cov", "-", "cov"], "nm"),
    ("search_dep0.pl", ["-"], "dep"),
])
def test_dash_reads_stdin(tmp_path, name, args, stdin):
    """CRAQ's scripts open their arguments with two-argument open(), where "-"
    is standard input; runAQI.sh relies on it."""
    rng = random.Random(1)
    dep, chroms = depth_table(rng)
    (tmp_path / "dep").write_text(dep)
    (tmp_path / "nm").write_text(nonmap(rng, dep))
    (tmp_path / "nm2").write_text(nonmap(rng, depth_table(random.Random(2))[0]))
    (tmp_path / "cov").write_text("".join(f"{c}\t{p}\t+\t{rng.randint(0, 4)}\t{rng.randint(1, 40)}\n"
                                          for c, p in sites(rng, chroms, 200)))
    (tmp_path / "bed").write_text("".join(f"chrA\t{p}\t{p + 900}\t0\n" for p in range(1, 6000, 1200)))
    (tmp_path / "err").write_text("chrA\t1250\t+\n")
    outs = []
    for src in (SRC, PATCH):
        with (tmp_path / stdin).open("rb") as fh:
            outs.append(subprocess.run(["perl", str(src / name), *args], stdin=fh, cwd=tmp_path,
                                       capture_output=True, check=True).stdout)
    if name == "search_uncertain_region.pl":
        outs = [sorted(o.splitlines()) for o in outs]
    assert outs[0] == outs[1] and outs[0], name
