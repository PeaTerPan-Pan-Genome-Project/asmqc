"""M8, M9 and M11 logic on hand-made inputs."""

import pytest

from asmqc import agp
from asmqc import m08_merqury as m08
from asmqc import m09_craq as m09
from asmqc import m11_homopolymer as m11


# --- M8 ----------------------------------------------------------------------
def test_coverage_peak():
    hist = [(1, 9000), (2, 300), (3, 50), (4, 80), (5, 400), (6, 900), (7, 600), (8, 100)]
    assert m08.coverage_peak(hist) == 6
    assert m08.coverage_peak([(1, 100), (2, 50), (3, 10)]) == 0  # no trough


# --- M11 ---------------------------------------------------------------------
def call(pos, ref, alt):
    return m11.Call("s", pos, ref, alt, 50.0, "1/1", 30)


# 0-based:      0123456789012345678901234567
SEQ = "GCTAAAAAAAAACGTAGAGAGAGTTACG"


@pytest.mark.parametrize("c,want", [
    (call(3, "T", "TA"), {"class": "hp", "base": "A", "run": 9, "size": 1, "dir": "ins"}),
    (call(3, "TAA", "T"), {"class": "hp", "base": "A", "run": 9, "size": 2, "dir": "del"}),
    (call(15, "T", "TAG"), {"class": "str2", "unit": "AG", "copies": 4, "size": 1,
                            "dir": "ins"}),
    (call(15, "TAGAG", "T"), {"class": "str2", "unit": "AG", "copies": 4, "size": 2,
                              "dir": "del"}),
    (call(1, "G", "GCC"), {"class": "other_indel"}),  # run of C would be 1+2 = 3 < 4
    (call(1, "G", "GCCC"), {"class": "hp", "base": "C", "run": 1, "size": 3, "dir": "ins"}),
    (call(24, "T", "TACGG"), {"class": "other_indel"}),
    (call(13, "C", "T"), {"class": "snv"}),
])
def test_classify_call(c, want):
    assert m11.classify_call(c, SEQ) == want


def test_bins_and_inside():
    assert [m11.run_bin(n) for n in (4, 8, 9, 12, 13, 20, 21, 99)] == [
        "4-8", "4-8", "9-10", "11-12", "13-15", "16-20", ">20", ">20"]
    assert [m11.change_bin(k) for k in (1, 4, 5, 9)] == ["1", "4", ">=5", ">=5"]
    bed = {"s": ([10, 100], [50, 200])}
    assert m11.inside(bed, "s", 10, 50) and not m11.inside(bed, "s", 45, 55)
    assert not m11.inside(bed, "s", 0, 5) and not m11.inside(bed, "x", 10, 20)


def test_regions_chunks_unplaced():
    lengths = {f"chr{i}": 100 for i in range(1, 8)} | {"u1": 30_000_000, "u2": 30_000_000,
                                                      "u3": 5}
    assert m11.regions(lengths)[7:] == [["u1", "u2"], ["u3"]]


# --- M9 ----------------------------------------------------------------------
REPORT = """Short Report:
#Chr\tCovered.Rate\tLow-confident.Rate\tAvg.CRH\tAvg.CSH\tAvg.CRE(R-AQI)\tAvg.CSE(S-AQI)
Genome\t0.97\t0.01\t0\t0\t12.7659(27.8987)\t0.5(95.12)
chr1\t0.98\t0\t0\t0\t13.1(26.9)\t0(100)
"""


def test_parse_report(tmp_path):
    p = tmp_path / "r"
    p.write_text(REPORT)
    r = m09.parse_report(p)
    assert r["r_aqi"] == pytest.approx(27.8987) and r["s_aqi"] == pytest.approx(95.12)
    p.write_text(REPORT.replace("Genome", "chrX"))
    with pytest.raises(ValueError):
        m09.parse_report(p)


def test_junction_crosscheck(tmp_path):
    p = tmp_path / "a.agp"
    lines = [
        "chr1\t1\t100000\t1\tW\tc1\t1\t100000\t+",
        "chr1\t100001\t100100\t2\tN\t100\tscaffold\tyes\tproximity_ligation",
        "chr1\t100101\t200000\t3\tW\tc2\t1\t99900\t+",
        "chr1\t200001\t300000\t4\tW\tc3\t1\t100000\t+",
    ]
    p.write_text("".join(f"{line}\n" for line in lines))
    rows, _ = agp.parse(p)
    js = m09.junctions(rows)
    assert js["chr1"] == [(100001, 100100), (200000, 200001)]
    cse = [("chr1", 95_000, 95_500), ("chr1", 150_000, 150_100), ("chr1", 205_000, 205_100)]
    out = m09.junction_crosscheck(cse, js, 10_000)
    assert [r[4] for r in out] == ["yes", "no", "yes"]
    assert out[0][3] == 4501
