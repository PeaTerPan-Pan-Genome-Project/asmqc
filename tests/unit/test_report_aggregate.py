"""report.html and asmqc aggregate on minimal result directories."""

import gzip
import json
from pathlib import Path

import pytest

from asmqc import aggregate, report, schema, summary
from asmqc import flags as fl


def result_dir(root: Path, label: str, version: str = "1.0.0") -> Path:
    """A result directory where only M1 ran (values NA)."""
    d = root / label
    d.mkdir(parents=True)
    row = {c: schema.NA for c in schema.HEADER}
    row |= {"label": label, "asmqc_version": version, "assembly_md5": "0" * 32,
            "run_date": "2026-10-01", "ena_rules": "PASS", "n_flags_ena_blocking": "0",
            "n_flags_warning": "1", "n_flags_info": "0"}
    row |= {f"{m}_status": "skipped_by_user" for m in schema.MODULES}
    row["m01_status"] = "failed"
    summary.write_summary(d / "qc_summary.tsv", row)
    fl.write(d / "flags.tsv", [fl.make("m05", "organelle_scaffold", "u1", 1, 10, "plastid")])
    (d / "run_manifest.json").write_text(json.dumps({
        "asmqc_version": version, "git_commit": "x", "command_line": "asmqc run",
        "start": "s", "end": "e", "host": {"cpu_model": "cpu", "threads": 1, "mem_gb": 1},
        "inputs": {"reads": []}, "tools": {"samtools": "1.24"},
        "reference_data": {"plastid": {"md5": "m"}}, "parameters": {},
        "modules": {"m01": {"status": "failed", "reason": "module did not finish"}}}))
    return d


def test_render_assembly_minimal(tmp_path):
    d = result_dir(tmp_path, "S1")
    html = report.render_assembly(d).read_text()
    assert "asmqc report: S1" in html and "organelle_scaffold" in html
    assert "module did not finish" in html
    assert "http://" not in html and "https://" not in html  # self-contained


def test_aggregate_ok(tmp_path):
    dirs = [result_dir(tmp_path, "S1"), result_dir(tmp_path, "S2", "1.0.3")]
    assert aggregate.run(dirs, tmp_path / "combined") == 0
    header, rows = summary.read_summary(tmp_path / "combined" / "qc_summary.tsv")
    assert header == schema.HEADER and [r["label"] for r in rows] == ["S1", "S2"]
    assert (tmp_path / "combined" / "report.html").exists()


@pytest.mark.parametrize("problem", ["version", "label", "header", "missing"])
def test_aggregate_refuses(tmp_path, problem):
    a = result_dir(tmp_path, "S1")
    if problem == "version":
        b = result_dir(tmp_path, "S2", "1.1.0")
    elif problem == "label":
        b = result_dir(tmp_path / "other", "S1")
    elif problem == "header":
        b = result_dir(tmp_path, "S2")
        text = (b / "qc_summary.tsv").read_text().replace("m01_n_seq", "m01_nseq")
        (b / "qc_summary.tsv").write_text(text)
    else:
        b = tmp_path / "nothing"
    with pytest.raises(aggregate.AggregateError):
        aggregate.collect([a, b])
    assert aggregate.run([a, b], tmp_path / "combined") == 1


def test_crosschecks(tmp_path):
    d = tmp_path / "r"
    (d / "m06_busco").mkdir(parents=True)
    (d / "m07_redundancy").mkdir()
    (d / "m06_busco" / "busco_derived.tsv").write_text(
        "busco_id\tstatus\tn_copies\tn_on_chromosomes\tn_on_unplaced\tsequences\n"
        "b1\tDuplicated\t2\t1\t1\tchr1;u1\n"
        "b2\tDuplicated\t2\t1\t1\tchr1;u2\n"
        "b3\tComplete\t1\t0\t1\tu1\n")
    (d / "m07_redundancy" / "redundancy.tsv").write_text(
        "scaffold\tlength\tclass\n" "u1\t5000\tduplicate\n" "u2\t5000\tunique\n")
    assert report.dup_buscos_on_duplicates(d) == 1

    (d / "m08_merqury").mkdir()
    (d / "m11_homopolymer").mkdir()
    with gzip.open(d / "m08_merqury" / "asm_only_kmers.bed.gz", "wt") as fh:
        fh.write("chr1\t100\t121\nchr1\t5000\t5021\nchr2\t100\t121\n")
    with gzip.open(d / "m11_homopolymer" / "errors.bed.gz", "wt") as fh:
        fh.write("chr1\t130\t130\thp\tA\t9\t+1A\t50\t30\n")
    assert report.asm_only_near_errors(d) == (1, 3)
