from pathlib import Path

import pytest

from asmqc import plan as planner
from asmqc import validate as v

CHROMS = [f"chr{i}" for i in range(1, 8)]


@pytest.fixture(autouse=True)
def avx2(monkeypatch):
    monkeypatch.setattr(v, "cpu_flags", lambda: {"avx2"})


def fasta(path: Path, names: list[str]) -> Path:
    path.write_text("".join(f">{n} desc\nACGT\n" for n in names))
    return path


def opts(tmp_path, names=CHROMS, **kw) -> v.RunOptions:
    return v.RunOptions(assembly=fasta(tmp_path / "a.fa", names), label="S1",
                        outdir=tmp_path / "out", **kw)


def problems(o) -> list[str]:
    with pytest.raises(v.ValidationError) as e:
        v.validate(o)
    return e.value.problems


def test_valid(tmp_path):
    assert v.validate(opts(tmp_path, CHROMS + ["scaf1"])) == []


def test_missing_chromosome_lists_names(tmp_path):
    (p,) = problems(opts(tmp_path, CHROMS[:6] + ["Chr7"]))
    assert "chr7" in p and "Chr7" in p


def test_mapping(tmp_path):
    names = [f"LG{i}" for i in range(1, 8)]
    mapping = {f"LG{i}": f"chr{i}" for i in range(1, 8)}
    assert v.validate(opts(tmp_path, names, chromosomes=mapping)) == []


def test_mapping_source_missing_and_clash(tmp_path):
    ps = problems(opts(tmp_path, CHROMS + ["x"], chromosomes={"y": "chr1", "x": "chr2"}))
    assert any("not in the assembly: y" in p for p in ps)
    assert any("duplicate names: chr2" in p for p in ps)


def test_parse_chromosome_map():
    assert v.parse_chromosome_map("a=chr1, b=chr2") == {"a": "chr1", "b": "chr2"}
    for bad in ("a", "a=", "a=chr1,a=chr2"):
        with pytest.raises(ValueError):
            v.parse_chromosome_map(bad)


def test_reads_need_declarations(tmp_path):
    r = tmp_path / "r.fq"
    r.write_text("")
    ps = problems(opts(tmp_path, ont=[r]))
    assert any("--reads-used-in-assembly" in p for p in ps)
    assert any("--ont-chemistry" in p for p in ps)
    ps = problems(opts(tmp_path, illumina=[r], reads_used_in_assembly="no"))
    assert ps == ["--illumina: files must come in R1,R2 pairs"]


def test_label_and_avx2(tmp_path, monkeypatch):
    monkeypatch.setattr(v, "cpu_flags", set)
    o = opts(tmp_path)
    o.label = "bad label"
    ps = problems(o)
    assert any("--label" in p for p in ps) and any("AVX2" in p for p in ps)


@pytest.mark.parametrize("reads,m08,m09,m11,ngs", [
    ({}, None, None, None, False),
    ({"illumina": 2}, "illumina", None, "illumina", False),
    ({"hifi": 1}, "hifi", "hifi", "hifi", False),
    ({"ont": 1}, None, "ont_r10", None, False),
    ({"illumina": 2, "hifi": 1, "ont": 1}, "illumina", "hifi", "illumina", True),
    ({"illumina": 2, "ont": 1}, "illumina", "ont_r10", "illumina", True),
])
def test_read_type_selection(tmp_path, reads, m08, m09, m11, ngs):
    kw = {k: [tmp_path / f"{k}{i}" for i in range(n)] for k, n in reads.items()}
    p = planner.make_plan(opts(tmp_path, ont_chemistry="r10", **kw))
    for m, want in (("m08", m08), ("m09", m09), ("m11", m11)):
        assert p[m].read_type == want
        assert p[m].status == ("run" if want else "skipped_no_input")
    assert p["m09"].short_reads == ngs
    assert all(p[m].status == "run" for m in planner.ASSEMBLY_ONLY)


def test_modules_selection(tmp_path):
    p = planner.make_plan(opts(tmp_path, modules=["1", "m04", "11"]))
    assert [m for m, x in p.items() if x.status == "run"] == ["m01", "m04"]
    assert p["m11"].status == "skipped_no_input"
    assert p["m02"].status == "skipped_by_user"
    with pytest.raises(ValueError):
        planner.make_plan(opts(tmp_path, modules=["3"]))
