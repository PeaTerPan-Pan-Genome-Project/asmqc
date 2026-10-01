"""Read-based modules on the synthetic data set (SPEC §11.1).

M11 always runs (its tools are in the dev env); M8 and M9 run when Merqury and
CRAQ are available (ASMQC_ENV_ROOT).
"""

import gzip
import json

import pytest

from asmqc import cli, summary, tools
from asmqc import validate as v

HAS_MERQURY = tools.tool_version("merqury") is not None
HAS_CRAQ = tools.tool_version("craq") is not None


@pytest.fixture(scope="module")
def result(testdata_reads, fake_refs, tmp_path_factory):
    out = tmp_path_factory.mktemp("out")
    d = testdata_reads
    modules = ",".join(["11"] + (["8"] if HAS_MERQURY else []) + (["9"] if HAS_CRAQ else []))
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(v, "cpu_flags", lambda: {"avx2"})
        mp.setenv("ASMQC_REFS", str(fake_refs))
        rc = cli.main(["run", "--assembly", str(d / "asm.fa"), "--agp", str(d / "asm.agp"),
                       "--illumina", f"{d}/reads_R1.fq.gz,{d}/reads_R2.fq.gz",
                       "--hifi", str(d / "reads_hifi.fq.gz"),
                       "--reads-used-in-assembly", "no",
                       "--label", "T", "--outdir", str(out), "--threads", "4",
                       "--mem-gb", "8", "--modules", modules])
    _, (row,) = summary.read_summary(out / "T" / "qc_summary.tsv")
    return rc, row, json.loads((d / "expected.json").read_text()), out / "T"


def test_exit_code(result):
    assert result[0] == 0


def test_m11_recovers_planted_errors(result):
    _, row, exp, res = result
    with gzip.open(res / "m11_homopolymer" / "errors.bed.gz", "rt") as fh:
        bed = [line.rstrip("\n").split("\t") for line in fh]
    hp = [b for b in bed if b[3] == "hp"]
    str2 = [b for b in bed if b[3] == "str2"]

    def found_hp(h):  # BED start is 0-based = 1-based run start - 1
        return any(b[0] == h["seq"] and int(b[1]) + 1 == h["run_start"]
                   and b[6] == f"+{h['change']}{h['base']}" and int(b[5]) == h["run_length"]
                   for b in hp)

    def found_str2(t):
        sign = "+" if t["direction"] == "ins" else "-"
        return any(b[0] == t["seq"] and int(b[1]) + 1 == t["run_start"]
                   and b[6] == f"{sign}{2 * abs(t['change_copies'])}"
                   and int(b[5]) == t["copies"] for b in str2)

    assert sum(map(found_hp, exp["m11"]["hp"])) >= 28
    assert sum(map(found_str2, exp["m11"]["str2"])) >= 4
    assert row["m11_status"] == "ok" and row["m11_read_type"] == "illumina"
    assert row["m11_reads_independent"] == "yes"


def test_m11_vcf_has_no_local_paths(result):
    *_, res = result
    with gzip.open(res / "m11_homopolymer" / "hom_calls.vcf.gz", "rt") as fh:
        header = [line for line in fh if line.startswith("##")]
    assert not any(str(res.parent) in line for line in header)


@pytest.mark.skipif(not HAS_MERQURY, reason="Merqury not available (set ASMQC_ENV_ROOT)")
def test_m08(result):
    _, row, _, res = result
    assert row["m08_status"] == "ok" and row["m08_k"] == "21"
    assert row["m08_read_type"] == "illumina" and row["m08_low_coverage"] == "no"
    assert float(row["m08_completeness_pct"]) > 95
    for f in ("spectra-cn.png", "spectra-asm.png", "asm_only_kmers.bed.gz"):
        assert (res / "m08_merqury" / f).exists()


@pytest.mark.skipif(not HAS_CRAQ, reason="CRAQ not available (set ASMQC_ENV_ROOT)")
def test_m09(result):
    _, row, _, _ = result
    assert row["m09_status"] == "ok" and row["m09_long_read_type"] == "hifi"
    assert row["m09_short_reads_used"] == "yes"
    assert row["m09_n_cse"] == "0"  # no structural error planted
    assert row["m09_cse_near_agp_junction"] == "0"


def test_no_local_paths_in_result(result, testdata_reads, fake_refs):
    """Nothing in the shared result names a local directory (SPEC §7.3)."""
    *_, res = result
    local = [str(res.parent), str(testdata_reads), str(fake_refs)]
    for f in res.rglob("*"):
        if not f.is_file() or "work" in f.relative_to(res).parts:
            continue
        if f.suffix == ".png":
            continue
        data = gzip.decompress(f.read_bytes()) if f.suffix == ".gz" else f.read_bytes()
        for p in local:
            assert p.encode() not in data, f"{p} in {f.relative_to(res)}"


def test_report(result):
    *_, res = result
    html = (res / "report.html").read_text()
    # M11 HP plot, plus the two Merqury spectra when M8 ran; M9 has no plot
    assert html.count("data:image/png;base64,") == (3 if HAS_MERQURY else 1)
    assert "src=\"http" not in html and "href=\"http" not in html
