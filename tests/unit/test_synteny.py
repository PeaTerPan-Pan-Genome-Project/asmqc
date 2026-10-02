"""M6 synteny against the Caméor v2 BUSCO anchors (report only)."""

import pytest

from asmqc import synteny
from asmqc.validate import CHROMOSOMES


def test_committed_anchor_table():
    lengths, anchors = synteny.read_anchors(synteny.anchors_path())
    assert sorted(lengths) == CHROMOSOMES and lengths["chr5"] == 656_251_770
    assert len(anchors) == 7404
    assert {a[0] for a in anchors.values()} == set(CHROMOSOMES)
    assert {a[2] for a in anchors.values()} == {"+", "-"}


def row(bid, seq, s, e, strand, status="Complete"):
    return {"busco_id": bid, "status": status, "sequence": seq, "Gene Start": str(s),
            "Gene End": str(e), "Strand": strand}


@pytest.fixture
def anchors():
    a = {f"f{i}": ("chr1", 1000.0 * (i + 1), "+") for i in range(10)}
    a |= {f"r{i}": ("chr2", 1000.0 * (i + 1), "+") for i in range(10)}
    a |= {"t1": ("chr3", 500.0, "-")}
    return a


def test_points_and_table(anchors):
    rows = [row(f"f{i}", "chr1", 100 * i, 100 * i + 50, "+") for i in range(10)]
    # chr2 reversed: positions decrease, strands flip
    rows += [row(f"r{i}", "chr2", 10_000 - 100 * i, 10_050 - 100 * i, "-") for i in range(10)]
    rows += [row("t1", "u9", 10, 60, "-"), row("x", "chr3", 1, 2, "+"),
             row("f0", "chr1", 1, 2, "+", status="Duplicated")]
    pts = synteny.points(rows, anchors)
    assert len(pts) == 21  # only Complete, only anchored
    assert pts[0] == ("f0", "chr1", 1000, "chr1", 25, "yes")
    tab = {r[0]: r for r in synteny.table(pts, {"chr1": 2000, "chr2": 20_000, "u9": 70})}
    assert tab["chr1"][3:6] == ("chr1", "1.000", "forward")
    assert tab["chr2"][3:6] == ("chr2", "1.000", "reverse")
    assert tab["chr3"][2] == 0
    assert tab["unplaced"][1:4] == (70, 1, "chr3")


def test_spearman():
    assert synteny.spearman([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1)
    assert synteny.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1)
    assert synteny.spearman([1, 2], [1, 2]) is None
