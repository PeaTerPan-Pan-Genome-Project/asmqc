import re
from pathlib import Path

import pytest

from asmqc import schema

SPEC = Path(__file__).parents[2] / "docs" / "design" / "SPEC.md"


def spec_columns() -> list[str]:
    """Column names from SPEC §7.2, in order."""
    text = SPEC.read_text()
    section = text[text.index("### 7.2"):text.index("### 7.3")]
    section = section[section.index("**Identity and run:**"):]
    names = []
    for bullet in section.split("\n- ")[1:]:
        bullet = re.sub(r"\([^)]*\)", "", bullet)  # enumerated values, examples
        bullet = re.split(r"(?<=`): ", bullet, maxsplit=1)[0]  # "`m11_status`: each one of …"
        bullet = bullet.split("\n\n", 1)[0]  # text after the list
        names += re.findall(r"`([a-z0-9_]+)`", bullet)
    return names


def test_header_matches_spec():
    assert schema.HEADER == spec_columns()


def test_every_module_has_status_and_columns():
    for m in schema.MODULES:
        assert f"{m}_status" in schema.HEADER
        assert schema.module_columns(m)


@pytest.mark.parametrize("column,value,expected", [
    ("m01_softmask_pct", 0.5119, "0.51"),
    ("m01_softmask_pct", 100, "100.00"),
    ("m01_round_lengths_expected", 0.009, "0.009"),
    ("m11_hp_errors_per_mb", 12.34567, "12.35"),
    ("m11_snv_per_mb", 0.000123456, "0.0001235"),
    ("m08_qv", 45.678, "45.68"),
    ("m01_n_seq", 16, "16"),
    ("m01_n_seq", 16.0, "16"),
    ("m01_agp_consistent", None, "NA"),
    ("m06_lineage", "fabales_odb12.2 2026-05-13", "fabales_odb12.2 2026-05-13"),
])
def test_format(column, value, expected):
    assert schema.format_value(column, value) == expected


@pytest.mark.parametrize("column,value", [
    ("m01_softmask_pct", 101),
    ("m01_n_seq", 1.5),
    ("m01_n_seq", True),
    ("ena_rules", "pass"),
    ("label", "a\tb"),
])
def test_format_rejects(column, value):
    with pytest.raises(ValueError):
        schema.format_value(column, value)
