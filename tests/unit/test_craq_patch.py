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
