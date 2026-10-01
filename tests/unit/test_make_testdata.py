"""The synthetic data set matches its own truth file (SPEC §11.1)."""

import gzip
import json
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[1]))
import make_testdata as mt


@pytest.fixture(scope="module")
def data(fake_refs, tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("td")
    mt.main([str(out), "--refs", str(fake_refs), "--coverage", "1"])
    return out


def fasta(path: Path) -> dict[str, str]:
    return dict(_records(path))


def _records(path):
    name, chunks = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                yield name, "".join(chunks)
            name, chunks = line[1:], []
        else:
            chunks.append(line)
    yield name, "".join(chunks)


def expected(data: Path) -> dict:
    return json.loads((data / "expected.json").read_text())


def test_deterministic(fake_refs, data, tmp_path):
    mt.main([str(tmp_path), "--refs", str(fake_refs), "--coverage", "1"])
    for f in ["asm.fa", "asm.agp", "expected.json", "reads_R1.fq.gz", "reads_R2.fq.gz"]:
        assert (tmp_path / f).read_bytes() == (data / f).read_bytes(), f


def test_chromosomes_and_lengths(data):
    seqs = fasta(data / "asm.fa")
    assert [n for n in seqs if n.startswith("chr")] == [f"chr{i}" for i in range(1, 8)]
    assert {n: len(s) for n, s in seqs.items()} == expected(data)["sequences"]
    for c, (length, _, _) in mt.CHROMS.items():
        assert len(seqs[c]) == length


def test_agp_covers_fasta(data):
    seqs = fasta(data / "asm.fa")
    rows = [line.split("\t") for line in (data / "asm.agp").read_text().splitlines()
            if not line.startswith("#")]
    by_obj: dict[str, list] = {}
    for r in rows:
        by_obj.setdefault(r[0], []).append(r)
    assert list(by_obj) == list(seqs)
    for obj, rs in by_obj.items():
        pos = 0
        for r in rs:
            beg, end = int(r[1]), int(r[2])
            assert beg == pos + 1
            if r[4] == "W":
                assert int(r[7]) - int(r[6]) == end - beg
                # only the planted 20 N and 5 N runs lie inside components
                runs = re.findall("N+", seqs[obj][beg - 1:end])
                assert {len(n) for n in runs} <= {5, 20}
            else:
                assert set(seqs[obj][beg - 1:end]) == {"N"}
            pos = end
        assert pos == len(seqs[obj])


def test_m01_truth(data):
    m = expected(data)["m01"]
    assert m["ena_rules"] == "FAIL"
    assert m["n_lt20bp"] == 1 and m["n_terminal_n"] == 1
    assert m["n_round_lengths"] == 5
    # 25000 (subseq A-1), 67000-type B, 3000 (component_beg-1), 5 round unplaced
    assert m["n_round_agp_cuts"] == 8
    assert m["n_iupac"] == 3


def test_telomere_arms(data):
    seqs = fasta(data / "asm.fa")
    e = expected(data)["m04"]
    assert e["capped_arms"] == 12 and e["wrong_orientation_arms"] == 1
    assert e["t2t_chromosomes"] == 5
    for chrom, arms in e["arms"].items():
        for arm, a in arms.items():
            if a["status"] == "absent":
                continue
            s = seqs[chrom][a["start"] - 1:a["end"]]
            assert s == a["motif"] * a["repeats"]
            capped_motif = mt.TELO_REV if arm == "start" else mt.TELO_FWD
            assert (a["motif"] == capped_motif) == (a["status"] == "capped")
            dist = a["start"] - 1 if arm == "start" else len(seqs[chrom]) - a["end"]
            assert dist == a["distance_from_end_bp"]


def test_planted_errors_in_assembly(data):
    seqs = fasta(data / "asm.fa")
    e = expected(data)["m11"]
    assert len(e["hp"]) == 30 and len(e["str2"]) == 5
    assert {h["base"] in "AT" for h in e["hp"]} == {True, False}
    for h in e["hp"]:
        s, i, b, n = seqs[h["seq"]], h["run_start"] - 1, h["base"], h["run_length"]
        assert s[i:i + n] == b * n and s[i - 1] != b and s[i + n] != b
    for t in e["str2"]:
        s, i, u, k = seqs[t["seq"]], t["run_start"] - 1, t["unit"], t["copies"]
        assert s[i:i + 2 * k] == u * k and s[i - 1] not in u and s[i + 2 * k] not in u
    assert {t["direction"] for t in e["str2"]} == {"ins", "del"}


def test_duplicate_is_exact_copy(data):
    seqs = fasta(data / "asm.fa")
    chrom, start, length = mt.DUP_SOURCE
    assert seqs["unplaced_dup"] == seqs[chrom][start:start + length]


def test_reads_paired(data):
    with gzip.open(data / "reads_R1.fq.gz", "rt") as f1, \
            gzip.open(data / "reads_R2.fq.gz", "rt") as f2:
        l1, l2 = f1.read().splitlines(), f2.read().splitlines()
    assert len(l1) == len(l2) == 4 * expected(data)["read_pairs"]
    assert l1[0].removesuffix("/1") == l2[0].removesuffix("/2")
    assert len(l1[1]) == mt.READ_LEN
