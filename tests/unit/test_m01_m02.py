"""Scan, AGP, length statistics, M1 and M2 on small hand-made inputs."""

import json
from pathlib import Path

import pytest

from asmqc import agp, m01_integrity, m02_contiguity, scan
from asmqc.stats import nx, pieces_between

CHROMS = "".join(f">chr{i}\n" + f"{'ACGT' * 15}\n" * 4 for i in range(1, 8))  # 240 bp


def run_scan(tmp_path: Path, text: bytes) -> Path:
    fa = tmp_path / "a.fa"
    fa.write_bytes(text)
    scan.scan(fa, tmp_path / "scan")
    return tmp_path / "scan"


def codes(flags) -> list[str]:
    return sorted(f.code for f in flags)


# --- stats -------------------------------------------------------------------
def test_nx():
    assert nx([], 0.5) == (0, 0)
    assert nx([10, 2, 8], 0.5) == (10, 1)
    assert nx([5, 5, 5, 5], 0.5) == (5, 2)
    assert nx([5, 5, 5, 5], 0.9) == (5, 4)


def test_pieces_between():
    assert pieces_between(100, []) == [100]
    assert pieces_between(100, [(1, 10), (51, 60), (91, 100)]) == [40, 30]


# --- scan --------------------------------------------------------------------
def test_scan_composition(tmp_path):
    d = run_scan(tmp_path, b">s1 desc\r\nNNacgtRY\r\nACGTX\r\n>s2\n\n")
    s1, s2 = scan.read_sequences(d)
    assert (s1.length, s1.n, s1.iupac, s1.other, s1.lower) == (13, 2, 2, 1, 4)
    assert (s1.leading_n, s1.trailing_n, s1.line_width, s1.width_consistent) == (2, 0, 8, True)
    assert s2.length == 0
    facts = json.loads((d / "file.json").read_text())
    assert facts["crlf"] is True


def test_scan_gaps_and_widths(tmp_path):
    d = run_scan(tmp_path, b">s1\nACGTNNNNN\nNNNNNACG\n>s2\nACGT\nACGTACGT\nA\n")
    assert scan.read_gaps(d) == [("s1", 5, 14)]
    s1, s2 = scan.read_sequences(d)
    assert s1.width_consistent and not s2.width_consistent


# --- AGP ---------------------------------------------------------------------
def write_agp(tmp_path, lines: list[str]) -> Path:
    p = tmp_path / "a.agp"
    p.write_text("##agp-version\t2.0\n" + "\n".join(lines) + "\n")
    return p


def test_agp_consistent_and_cuts(tmp_path):
    p = write_agp(tmp_path, [
        "s1\t1\t50\t1\tW\tc1_subseq_1001:1050\t1\t50\t+",
        "s1\t51\t150\t2\tN\t100\tscaffold\tyes\tproximity_ligation",
        "s1\t151\t250\t3\tW\tc2\t2001\t2100\t-",
    ])
    rows, bad = agp.parse(p)
    assert bad == [] and agp.check(rows, {"s1": 250}) == []
    cuts = [c for _, c in agp.cut_coordinates(rows)]
    assert cuts == [1000, 1050, 2000, 2100]


def test_agp_problems(tmp_path):
    p = write_agp(tmp_path, [
        "scaffold_1\t1\t100\t1\tW\tc1\t1\t100\t+",  # name differs from FASTA
        "s2\t1\t50\t1\tW\tc2\t1\t60\t+",  # component longer than span
        "s2\t52\t60\t2\tW\tc3\t1\t9\t+",  # not contiguous
        "s3\t1\t10\t1\tW\tc4",  # malformed
    ])
    rows, bad = agp.parse(p)
    assert len(bad) == 1
    problems = agp.check(rows, {"s1": 100, "s2": 60})
    text = "\n".join(f"{o}: {m}" for o, m in problems)
    assert "scaffold_1: AGP object not in the FASTA" in text
    assert "s1: FASTA sequence not in the AGP" in text
    assert "length 60 != span 50" in text and "starts at 52, expected 51" in text


# --- M1 ----------------------------------------------------------------------
def test_m01_flags(tmp_path):
    fa = (CHROMS + ">short\nACGTACGTAC\n>tn\nACGTACGT" + "A" * 200 + "NNN\n"
          + ">mostlyN\n" + "N" * 300 + "ACGT\n>round\n" + "A" * 1000 + "\n>chr1\nACGT\n")
    d = run_scan(tmp_path, fa.encode())
    values, flags, _ = m01_integrity.evaluate(d, None)
    assert codes(flags) == ["duplicate_name", "n_fraction_gt_50pct", "round_length",
                            "seq_lt_20bp", "seq_lt_20bp", "terminal_n", "terminal_n"]
    assert values["m01_n_lt20bp"] == 2 and values["m01_n_lt200bp"] == 2  # 10 bp, 4 bp
    assert values["m01_n_terminal_n"] == 2  # 'tn' trailing, 'mostlyN' leading
    assert values["m01_n_duplicate_names"] == 1
    assert values["m01_agp_consistent"] is None and values["m01_n_agp_cuts"] is None
    assert values["m01_n_round_lengths"] == 1


def test_m01_agp_mismatch(tmp_path):
    d = run_scan(tmp_path, CHROMS.encode())
    p = write_agp(tmp_path, [f"Chr{i}\t1\t240\t1\tW\tc{i}\t1\t240\t+" for i in range(1, 8)])
    values, flags, _ = m01_integrity.evaluate(d, p)
    assert values["m01_agp_consistent"] == "no"
    assert codes(flags).count("agp_mismatch") == 14
    assert values["m01_n_agp_cuts"] == 7 and values["m01_n_round_agp_cuts"] == 0


# --- M2 ----------------------------------------------------------------------
def test_m02_method_and_values(tmp_path):
    fa = "".join(f">chr{i}\n{'A' * 100}{'N' * 20}{'C' * 80}\n" for i in range(1, 8))
    d = run_scan(tmp_path, (fa + ">u1\n" + "G" * 50 + "\n").encode())
    p = write_agp(tmp_path, [
        line for i in range(1, 8) for line in (
            f"chr{i}\t1\t100\t1\tW\tc{i}a\t1\t100\t+",
            f"chr{i}\t101\t120\t2\tN\t20\tscaffold\tyes\tproximity_ligation",
            f"chr{i}\t121\t200\t3\tW\tc{i}b\t1\t80\t+")
    ] + ["u1\t1\t50\t1\tW\tc8\t1\t50\t+"])
    v, per_chrom, _, _ = m02_contiguity.evaluate(d, p)
    assert v["m02_contig_method"] == "agp" and v["m02_n_contigs"] == 15
    assert v["m02_chrom_bp"] == 1400 and v["m02_n_unplaced"] == 1
    assert v["m02_n_gaps"] == 7 and v["m02_gap_bp"] == 140
    assert v["m02_contig_n50"] == v["m02_contig_n50_nsplit10"] == 100
    assert per_chrom[0] == ("chr1", 200, 2, 1, 20, 100)
    v2, *_ = m02_contiguity.evaluate(d, None)
    assert v2["m02_contig_method"] == "nsplit10" and v2["m02_n_contigs"] == 15


def test_m02_quast_disagreement():
    quast = {"N50": ("10", "5"), "L50": ("1", "1"), "N90": ("8", "5")}
    assert m02_contiguity.quast_disagreements(quast, [10, 8], [5]) == []
    problems = m02_contiguity.quast_disagreements(quast, [10, 8], [6])
    assert problems and "nsplit10 contig N50" in problems[0]


@pytest.mark.parametrize("bad", ["scaffold N50", "nsplit10"])
def test_m02_fails_on_disagreement(tmp_path, bad):
    d = run_scan(tmp_path, CHROMS.encode())
    report = tmp_path / "report.tsv"
    asm, broken = ("1", "240") if bad == "scaffold N50" else ("240", "1")
    report.write_text(f"Assembly\tasm\tasm_broken\nN50\t{asm}\t{broken}\n"
                      "L50\t4\t4\nN90\t240\t240\n")
    with pytest.raises(SystemExit, match="disagree"):
        m02_contiguity.main(["--scan", str(d), "--quast-report", str(report),
                             "--outdir", str(tmp_path / "o"), "--work", str(tmp_path / "w")])
