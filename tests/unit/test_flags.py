import pytest

from asmqc import flags as fl


def test_severity_from_code():
    assert fl.make("m01", "seq_lt_20bp").severity == fl.ENA_BLOCKING
    assert fl.make("m08", "low_coverage").severity == fl.WARNING


def test_wrong_module_rejected():
    with pytest.raises(ValueError):
        fl.make("m04", "seq_lt_20bp")


def test_round_trip(tmp_path):
    flags = [fl.make("m01", "terminal_n", "s1", 100, 137, 37, "trailing N"),
             fl.make("m01", "iupac_present", "chr2", value=3)]
    fl.write(tmp_path / "f.tsv", flags)
    assert fl.read(tmp_path / "f.tsv") == flags


def test_read_rejects_wrong_severity(tmp_path):
    (tmp_path / "f.tsv").write_text(
        "\t".join(fl.COLUMNS) + "\nm01\tseq_lt_20bp\tINFO\ts\t\t\t\t\n")
    with pytest.raises(ValueError):
        fl.read(tmp_path / "f.tsv")
