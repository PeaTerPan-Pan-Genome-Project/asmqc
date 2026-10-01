"""`asmqc run` end to end on a tiny assembly (Snakemake, prep, merge)."""

import gzip
import hashlib
import json
from pathlib import Path

import pytest

from asmqc import cli, prep, runner, schema
from asmqc import validate as v

WORKFLOW_RULES = Path(runner.asmqc_home()) / "workflow" / "rules"


@pytest.fixture(autouse=True)
def avx2(monkeypatch):
    monkeypatch.setattr(v, "cpu_flags", lambda: {"avx2"})


@pytest.fixture
def assembly(tmp_path) -> Path:
    text = "".join(f">LG{i} x\n{'ACGT' * 50}\n" for i in range(1, 8)) + ">s1\nACGTN\n"
    path = tmp_path / "asm.fa.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    return path


def test_prep_renames_and_hashes_decompressed(assembly, tmp_path):
    md5 = prep.prep_fasta(assembly, tmp_path / "out.fa", {"LG1": "chr1"})
    assert md5 == hashlib.md5(gzip.decompress(assembly.read_bytes())).hexdigest()
    lines = (tmp_path / "out.fa").read_text().splitlines()
    assert lines[0] == ">chr1 x" and lines[2] == ">LG2 x"


def test_command_line_hides_paths():
    argv = ["run", "--assembly", "/data/x/asm.fa", "--illumina=/d/r1.fq,/d/r2.fq",
            "--label", "S1"]
    assert runner.command_line(argv, False) == (
        "asmqc run --assembly asm.fa --illumina=r1.fq,r2.fq --label S1")
    assert "/data/x/asm.fa" in runner.command_line(argv, True)


def test_validation_failure_exits_1(assembly, tmp_path):
    rc = cli.main(["run", "--assembly", str(assembly), "--label", "S1",
                   "--outdir", str(tmp_path / "out")])
    assert rc == 1  # LG1..LG7 not mapped
    assert not (tmp_path / "out").exists()


def test_dry_run(assembly, tmp_path):
    rc = cli.main(["run", "--assembly", str(assembly), "--label", "S1", "--dry-run",
                   "--outdir", str(tmp_path / "out"),
                   "--chromosomes", ",".join(f"LG{i}=chr{i}" for i in range(1, 8))])
    assert rc == 0 and not (tmp_path / "out").exists()


def test_run(assembly, tmp_path):
    out = tmp_path / "out"
    rc = cli.main(["run", "--assembly", str(assembly), "--label", "S1",
                   "--outdir", str(out), "--threads", "2", "--mem-gb", "2",
                   "--modules", "1", "--keep-intermediates",
                   "--chromosomes", ",".join(f"LG{i}=chr{i}" for i in range(1, 8))])
    res = out / "S1"
    header, (row,) = runner.summary.read_summary(res / "qc_summary.tsv")
    assert header == schema.HEADER
    assert row["assembly_md5"] == hashlib.md5(gzip.decompress(assembly.read_bytes())).hexdigest()
    m01_implemented = (WORKFLOW_RULES / "m01.smk").exists()
    assert row["m01_status"] == ("ok" if m01_implemented else "failed")
    assert rc == (0 if m01_implemented else 2)
    assert all(row[f"{m}_status"] == "skipped_by_user" for m in schema.MODULES if m != "m01")
    assert (res / "work" / "prep" / "asm.fa").read_text().startswith(">chr1 x\n")

    manifest = json.loads((res / "run_manifest.json").read_text())
    assert manifest["inputs"]["assembly"]["name"] == "asm.fa.gz"
    assert str(tmp_path) not in manifest["command_line"]
    assert manifest["chromosome_map"]["LG7"] == "chr7"
    assert manifest["parameters"]["m02"] == {"nsplit_min_gap": 10}
