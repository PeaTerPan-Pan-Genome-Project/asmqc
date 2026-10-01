"""asmqc run on the synthetic data set: M1 (and M2 when QUAST is available)
against expected.json (SPEC §11.1)."""

import json
import shutil

import pytest

from asmqc import cli, summary, tools
from asmqc import validate as v

HAS_QUAST = tools.tool_version("quast") is not None


@pytest.fixture(scope="module")
def result(testdata_noreads, tmp_path_factory):
    out = tmp_path_factory.mktemp("out")
    d = testdata_noreads
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(v, "cpu_flags", lambda: {"avx2"})
        rc = cli.main(["run", "--assembly", str(d / "asm.fa"), "--agp", str(d / "asm.agp"),
                       "--label", "T", "--outdir", str(out), "--threads", "2",
                       "--mem-gb", "2", "--modules", "1,2" if HAS_QUAST else "1"])
    _, (row,) = summary.read_summary(out / "T" / "qc_summary.tsv")
    return rc, row, json.loads((d / "expected.json").read_text()), out / "T"


def test_exit_code_and_cleanup(result):
    rc, _, _, res = result
    assert rc == 0  # an ENA FAIL is a result, not an error
    assert not (res / "work").exists()


def test_m01_matches_truth(result):
    _, row, exp, _ = result
    for key, want in exp["m01"].items():
        if key == "flags_ena_blocking":
            continue
        col = key if key == "ena_rules" else f"m01_{key}"
        if isinstance(want, float):
            assert float(row[col]) == pytest.approx(want, abs=0.005), col
        else:
            assert row[col] == str(want), col


def test_ena_blocking_flags(result):
    _, _, exp, res = result
    rows = [line.split("\t") for line in (res / "flags.tsv").read_text().splitlines()[1:]]
    blocking: dict[str, list[str]] = {}
    for module, code, severity, seq_id, *_ in rows:
        if severity == "ENA_BLOCKING":
            blocking.setdefault(code, []).append(seq_id)
    assert blocking == exp["m01"]["flags_ena_blocking"]


@pytest.mark.skipif(not HAS_QUAST, reason="QUAST not available (set ASMQC_ENV_ROOT)")
def test_m02(result):
    _, row, exp, res = result
    chrom_bp = sum(n for s, n in exp["sequences"].items() if s.startswith("chr"))
    assert row["m02_status"] == "ok" and row["m02_contig_method"] == "agp"
    assert row["m02_chrom_bp"] == str(chrom_bp)
    # the 20 N run inside a chr6 component splits it only for nsplit10
    assert row["m02_contig_n50_nsplit10"] != row["m02_contig_n50"]
    assert row["m02_n_contigs"] == str(25)
    assert (res / "m02_contiguity" / "quast_report.tsv").exists()


def test_outputs_present(result):
    *_, res = result
    for f in ["m01_integrity/integrity.tsv", "m01_integrity/sequences.tsv",
              "m01_integrity/gaps.tsv", "m01_integrity/checksums.md5"]:
        assert (res / f).exists(), f
    assert shutil.which("samtools")
