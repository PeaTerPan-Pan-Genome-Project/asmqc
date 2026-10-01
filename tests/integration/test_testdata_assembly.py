"""asmqc run on the synthetic data set, assembly-only modules, against
expected.json (SPEC §11.1). M2 runs when QUAST is available (ASMQC_ENV_ROOT);
M6 is not run, as in the smoke test (no real genes)."""

import json
import shutil

import pytest

from asmqc import cli, summary, tools
from asmqc import validate as v

HAS_QUAST = tools.tool_version("quast") is not None


@pytest.fixture(scope="module")
def result(testdata_noreads, fake_refs, tmp_path_factory):
    out = tmp_path_factory.mktemp("out")
    d = testdata_noreads
    modules = "1,2,4,5,7" if HAS_QUAST else "1,4,5,7"
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(v, "cpu_flags", lambda: {"avx2"})
        mp.setenv("ASMQC_REFS", str(fake_refs))
        rc = cli.main(["run", "--assembly", str(d / "asm.fa"), "--agp", str(d / "asm.agp"),
                       "--label", "T", "--outdir", str(out), "--threads", "2",
                       "--mem-gb", "2", "--modules", modules])
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


def test_m04(result):
    _, row, exp, res = result
    e = exp["m04"]
    for k in ("capped_arms", "t2t_chromosomes", "wrong_orientation_arms",
              "unplaced_with_telomere"):
        assert row[f"m04_{k}"] == str(e[k]), k
    assert row["m04_interstitial_arrays"] == str(len(e["interstitial"]))
    got = {(r[0], r[1]): r[2] for r in
           (line.split("\t") for line in
            (res / "m04_telomeres" / "telomeres.tsv").read_text().splitlines()[1:])}
    want = {(c, a): s["status"] for c, arms in e["arms"].items() for a, s in arms.items()}
    assert got == want
    png = res / "m04_telomeres" / "karyoplot.png"
    assert png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_m05(result):
    _, row, exp, _ = result
    e = exp["m05"]
    assert row["m05_plastid_scaffolds_n"] == str(len(e["plastid_scaffolds"]))
    assert row["m05_mito_scaffolds_n"] == "0"
    insert = e["chrom_plastid_insert"][0]
    # alignment ends may extend a few bp where a flanking base matches by chance
    want = insert["end"] - insert["start"] + 1
    assert abs(int(row["m05_chrom_plastid_like_bp"]) - want) <= 20
    assert row["m05_rdna45s_copies"] == str(e["rdna45s"][0]["copies"])
    assert row["m05_rdna5s_copies"] == str(e["rdna5s"][0]["copies"])
    assert row["m05_rdna45s_loci"] == e["rdna45s"][0]["seq"]
    assert row["m05_rdna5s_loci"] == e["rdna5s"][0]["seq"]
    assert row["m05_rdna_only_scaffolds_n"] == "0"


def test_m07(result):
    _, row, exp, res = result
    rows = [line.split("\t") for line in
            (res / "m07_redundancy" / "redundancy.tsv").read_text().splitlines()[1:]]
    assert {r[0]: r[2] for r in rows} == exp["m07"]["classes"]
    assert row["m07_n_duplicate"] == "1"


def test_report(result):
    *_, res = result
    html = (res / "report.html").read_text()
    assert "Telomere karyoplot" in html and "Unplaced bp per class" in html
    assert ("Cumulative contig length" in html) == HAS_QUAST
