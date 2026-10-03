"""M4, M5, M6, M7 logic on hand-made tool output."""

import json

import pytest

from asmqc import m04_telomeres as m04
from asmqc import m05_organelle_rdna as m05
from asmqc import m06_busco as m06
from asmqc import m07_redundancy as m07
from asmqc.validate import CHROMOSOMES

LEN = 300_000


# --- M4 ----------------------------------------------------------------------
def windows(spec: dict[str, list[tuple[int, int, int]]]) -> list:
    """spec: seq -> [(window_end, fwd, rev)] for non-empty windows; others are 0."""
    rows = []
    for seq in CHROMOSOMES:
        nz = {w: (f, r) for w, f, r in spec.get(seq, [])}
        for end in range(10_000, LEN + 1, 10_000):
            f, r = nz.get(end, (0, 1))
            rows.append((seq, end, f, r))
    for seq, items in spec.items():
        if seq not in CHROMOSOMES:
            rows += [(seq, w, f, r) for w, f, r in items]
    return rows


def test_m04_arms_bands_and_interstitial():
    spec = {c: [(10_000, 0, 500), (LEN, 600, 0)] for c in CHROMOSOMES}
    spec["chr1"] = [(10_000, 0, 300), (20_000, 2, 200), (LEN, 400, 1)]  # band of 2 windows
    # start band 60 kb in: start arm absent, and the band is interstitial
    spec["chr2"] = [(70_000, 0, 500), (LEN, 600, 0)]
    spec["chr3"] = [(10_000, 500, 0), (LEN, 0, 500)]  # both wrong orientation
    spec["chr4"] = [(10_000, 0, 500), (150_000, 100, 0), (LEN, 600, 0)]  # interstitial
    spec["u1"] = [(10_000, 0, 30), (25_000, 0, 0)]
    lengths = dict.fromkeys(CHROMOSOMES, LEN) | {"u1": 25_000}
    values, flags, arms, inter, unpl, by_seq = m04.evaluate(windows(spec), lengths)
    status = {(a[0], a[1]): a[2] for a in arms}
    assert by_seq["chr1"][0].start == 1 and by_seq["chr1"][0].end == 20_000
    assert status[("chr2", "start")] == "absent"
    assert status[("chr3", "start")] == status[("chr3", "end")] == "wrong_orientation"
    assert values == {"m04_capped_arms": 11, "m04_t2t_chromosomes": 5,
                      "m04_wrong_orientation_arms": 2, "m04_interstitial_arrays": 2,
                      "m04_unplaced_with_telomere": 1}
    assert inter == [("chr2", 60_001, 70_000, 0, 500, 3500),
                     ("chr4", 140_001, 150_000, 100, 0, 700)]
    assert [u[2] for u in unpl] == ["start", "end"]  # one band near both ends of u1
    assert sorted(f.code for f in flags) == ["interstitial_telomere"] * 2 + [
        "wrong_orientation_telomere"] * 2


def test_m04_t2t_requires_no_gap(tmp_path):
    spec = {c: [(10_000, 0, 500), (LEN, 600, 0)] for c in CHROMOSOMES}
    lengths = dict.fromkeys(CHROMOSOMES, LEN)
    gaps_tsv = tmp_path / "gaps.tsv"
    gaps_tsv.write_text("seq_id\tstart\tend\tlength\n"
                        "chr1\t100\t199\t100\nchr1\t500\t509\t10\nchr3\t7\t16\t10\n")
    gaps = m04.gap_counts(gaps_tsv)
    assert gaps == {"chr1": 2, "chr3": 1}
    values, _, arms, *_ = m04.evaluate(windows(spec), lengths, gaps)
    assert values["m04_capped_arms"] == 14 and values["m04_t2t_chromosomes"] == 5
    assert [c for c, ok in m04.t2t(arms, gaps).items() if not ok] == ["chr1", "chr3"]


# --- M5 ----------------------------------------------------------------------
def paf(q, qlen, qs, qe, t, matches, block, mapq=60, tp="P"):
    return "\t".join(map(str, [q, qlen, qs, qe, "+", t, 1000, 0, 100, matches, block, mapq,
                               f"tp:A:{tp}"]))


def test_m05_organelles(tmp_path):
    lines = [
        paf("u1", 1000, 0, 500, "NC_014057.1", 490, 500),
        paf("u1", 1000, 400, 900, "NC_014057.1", 495, 500),  # union 900 bp: 90 %
        paf("u2", 1000, 0, 1000, "PP555264.1", 900, 1000),  # identity 0.90: ignored
        paf("chr1", LEN, 100, 2100, "PP555264.1", 2000, 2000),
        paf("chr1", LEN, 5000, 6000, "PP555264.1", 1000, 1000),
    ]
    p = tmp_path / "o.paf"
    p.write_text("\n".join(lines) + "\n")
    lengths = dict.fromkeys(CHROMOSOMES, LEN) | {"u1": 1000, "u2": 1000}
    scaffolds, on_chrom, flags = m05.organelles(p, m05.ORGANELLES, lengths)
    assert scaffolds == [("u1", 1000, "0.9000", "0.0000", "plastid")]
    assert on_chrom[0] == ("chr1", 0, 0, 3000, 2000)
    assert [f.code for f in flags] == ["organelle_scaffold"]


def blast(q, s, ss, se, qcov=1.0):
    """BLAST outfmt "6 std qlen slen" line; the query (subunit) is 100 bp."""
    return "\t".join(map(str, [q, s, 99, 100, 0, 0, 1, round(100 * qcov), ss, se, 0, 200,
                               100, 1000]))


def test_m05_rdna(tmp_path):
    q18, q25, q5 = "a#rDNA/45S_rDNA/18S", "b#rDNA/45S_rDNA/25S", "c#rDNA/5S_rDNA/5S"
    lines = []
    for i in range(3):  # 45S units every 10 kb on chr4; two 18S variants hit each unit
        lines += [blast(q18, "chr4", 1000 + i * 10_000, 2800 + i * 10_000),
                  blast("d#rDNA/45S_rDNA/18S", "chr4", 1100 + i * 10_000, 2700 + i * 10_000),
                  blast(q25, "chr4", 3500 + i * 10_000, 6900 + i * 10_000)]
    lines += [blast(q18, "chr4", 100_000, 101_800)]  # > 20 kb away: a 1-copy fragment
    lines += [blast(q5, "u1", 400 + i * 350, 280 + i * 350) for i in range(12)]  # minus strand
    lines += [blast(q5, "chr5", 9000 + i * 350, 9100 + i * 350) for i in range(4)]  # 5S < 10
    # hits covering < 50 % of the subunit are not copies
    lines += [blast(q5, "chr2", 5000 + i * 350, 5040 + i * 350, qcov=0.4) for i in range(9)]
    p = tmp_path / "b.tsv"
    p.write_text("\n".join(lines) + "\n")
    lengths = dict.fromkeys(CHROMOSOMES, LEN) | {"u1": 4800}
    arrays = m05.rdna_arrays(p, lengths)
    assert arrays == [("chr4", 1000, 26_900, "45S", 3, "array", "chromosome"),
                      ("chr4", 100_000, 101_800, "45S", 1, "fragment", "chromosome"),
                      ("chr5", 9000, 10_150, "5S", 4, "fragment", "chromosome"),
                      ("u1", 280, 4250, "5S", 12, "array", "rdna_only_scaffold")]  # 82.7 %
    assert m05.loci(arrays, "45S") == "chr4" and m05.loci(arrays, "5S") == "unplaced"


# --- M6 ----------------------------------------------------------------------
def test_m06_derived():
    rows = [
        {"busco_id": "b1", "status": "Complete", "sequence": "chr1"},
        {"busco_id": "b2", "status": "Complete", "sequence": "u1"},
        {"busco_id": "b3", "status": "Duplicated", "sequence": "chr1"},
        {"busco_id": "b3", "status": "Duplicated", "sequence": "chr2"},
        {"busco_id": "b4", "status": "Duplicated", "sequence": "chr1"},
        {"busco_id": "b4", "status": "Duplicated", "sequence": "u2"},
        {"busco_id": "b5", "status": "Fragmented", "sequence": "u3"},
        {"busco_id": "b6", "status": "Missing", "sequence": ""},
    ]
    counts, table = m06.derived(rows)
    assert counts == {"m06_complete_on_unplaced": 2, "m06_dup_both_on_chrom": 1,
                      "m06_dup_any_on_unplaced": 1}
    assert table[3] == ("b4", "Duplicated", 2, 1, 1, "chr1;u2")


def test_m06_full_table_busco6(tmp_path):
    """BUSCO 6.1.0 layout: 10 columns for found BUSCOs, 2 for missing."""
    p = tmp_path / "full_table.tsv"
    p.write_text(
        "# BUSCO version is: 6.1.0 \n"
        "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength"
        "\tOrthoDB url\tDescription\n"
        "267at72025\tComplete\tchr3\t829659\t830284\t+\t270.3\t132"
        "\thttps://v12-2.orthodb.org/?query=267at72025\tGSH-induced LITAF domain protein\n"
        "5at72025\tDuplicated\tchr1\t10\t900\t-\t500.0\t300\turl\tdesc\n"
        "5at72025\tDuplicated\tu7\t10\t900\t-\t500.0\t300\turl\tdesc\n"
        "27at72025\tMissing\n")
    rows = m06.read_full_table(p)
    assert [r["status"] for r in rows] == ["Complete", "Duplicated", "Duplicated", "Missing"]
    counts, _ = m06.derived(rows)
    assert counts == {"m06_complete_on_unplaced": 1, "m06_dup_both_on_chrom": 0,
                      "m06_dup_any_on_unplaced": 1}


def test_m06_busco_cds(tmp_path):
    d = tmp_path / "single_copy_busco_sequences"
    d.mkdir()
    (d / "10at72025.gff").write_text(
        "chr3\tminiprot\tmRNA\t100\t400\t756\t-\t.\tID=MP1\n"
        "chr3\tminiprot\tCDS\t300\t400\t756\t-\t0\tParent=MP1\n"
        "chr3\tminiprot\tCDS\t100\t200\t756\t-\t1\tParent=MP1\n"
        "chr3\tminiprot\tstop_codon\t97\t99\t0\t-\t0\tParent=MP1\n")
    (d / "9at72025.gff").write_text("chr1\tminiprot\tCDS\t1\t30\t9\t+\t0\tParent=MP2\n")
    rows = m06.cds_rows(d)
    assert rows == [("chr1", 0, 30, "9at72025", "+"), ("chr3", 99, 200, "10at72025", "-"),
                    ("chr3", 299, 400, "10at72025", "-")]
    m06.write_cds_bed(rows, tmp_path / "cds.bed.gz")
    first = (tmp_path / "cds.bed.gz").read_bytes()
    m06.write_cds_bed(rows, tmp_path / "cds.bed.gz")
    assert (tmp_path / "cds.bed.gz").read_bytes() == first  # deterministic gzip


def test_m06_lineage_checks(tmp_path):
    (tmp_path / "dataset.cfg").write_text(
        "name=fabales_odb12.2\ncreation_date=2026-05-13\nnumber_of_BUSCOs=7702\n")
    assert m06.check_lineage(tmp_path) == []
    (tmp_path / "dataset.cfg").write_text(
        "name=fabales_odb12.2\ncreation_date=2025-01-01\nnumber_of_BUSCOs=7702\n")
    assert "creation_date" in m06.check_lineage(tmp_path)[0]
    s = {"lineage_dataset": {"name": "fabales_odb12.2", "creation_date": "2026-05-13"},
         "parameters": {"gene_predictor": "metaeuk"}, "versions": {"busco": "6.1.0"},
         "results": {"n_markers": 7702}}
    assert m06.check_summary(s) == ["gene predictor metaeuk, expected miniprot"]
    json.dumps(s)


# --- M7 ----------------------------------------------------------------------
def test_m07_classes():
    lines = [
        paf("dup", 10_000, 0, 9_500, "chr1", 9_500, 9_500),  # 95 % in one chain
        paf("part", 10_000, 0, 6_000, "chr1", 6_000, 6_000),  # 60 %
        paf("rep", 10_000, 0, 3_000, "chr2", 3_000, 3_000, mapq=0),
        paf("rep", 10_000, 3_000, 6_000, "chr3", 3_000, 3_000, mapq=0, tp="S"),
        paf("lowid", 10_000, 0, 10_000, "chr1", 9_500, 10_000),  # identity 0.95
        paf("pieces", 10_000, 0, 3_000, "chr1", 3_000, 3_000),
        paf("pieces", 10_000, 5_000, 8_000, "chr2", 3_000, 3_000),  # 2 chains, 60 % union
    ]
    lengths = dict.fromkeys(CHROMOSOMES, LEN) | {
        "dup": 10_000, "part": 10_000, "rep": 10_000, "lowid": 10_000, "pieces": 10_000,
        "tiny": 900}
    classes = m07.classify(lines, lengths)
    assert {s: c[0] for s, c in classes.items()} == {
        "dup": "duplicate", "part": "partial_overlap", "rep": "repeat_like",
        "lowid": "unique", "pieces": "repeat_like", "tiny": "short"}
    assert classes["dup"][1:4] == ("chr1", 1, 100)  # paf() target span is 0-100
    v = m07.values_from(classes, lengths, total_bp=1_000_000)
    assert v["m07_duplicate_bp"] == 10_000 and v["m07_total_minus_duplicate_bp"] == 990_000
    assert v["m07_n_repeat_like"] == 2 and v["m07_short_bp"] == 900


@pytest.mark.parametrize("n", [0, 1])
def test_m07_skips_without_unplaced(tmp_path, n):
    fai = tmp_path / "a.fai"
    rows = [f"{c}\t{LEN}\t0\t60\t61" for c in CHROMOSOMES] + ["u\t500\t0\t60\t61"] * n
    fai.write_text("\n".join(rows) + "\n")
    m07.main(["lists", "--fai", str(fai), "--work", str(tmp_path / "w")])
    m07.main(["classify", "--fai", str(fai), "--paf", str(tmp_path / "none.paf"),
              "--outdir", str(tmp_path / "o"), "--work", str(tmp_path / "w")])
    s = json.loads((tmp_path / "w" / "summary.json").read_text())
    assert s["status"] == "skipped_no_input"
