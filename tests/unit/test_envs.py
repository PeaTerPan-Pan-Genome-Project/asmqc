"""Environment files: exact pins, and dev.yaml consistent with core.yaml."""

from pathlib import Path

import pytest
import yaml

ENVS = Path(__file__).parents[2] / "envs"
IMAGE_ENVS = ["core", "quast", "busco", "merqury", "craq"]


def pins(name: str) -> dict[str, str]:
    deps = yaml.safe_load((ENVS / f"{name}.yaml").read_text())["dependencies"]
    return dict(d.split("=", 1) for d in deps)


@pytest.mark.parametrize("name", IMAGE_ENVS + ["dev"])
def test_every_dependency_pinned_exactly(name):
    deps = yaml.safe_load((ENVS / f"{name}.yaml").read_text())["dependencies"]
    for d in deps:
        assert isinstance(d, str) and "=" in d and not d.endswith("*"), d


@pytest.mark.parametrize("name", IMAGE_ENVS)
def test_lock_exists_and_explicit(name):
    text = (ENVS / f"{name}.lock").read_text()
    assert "@EXPLICIT" in text
    for line in text.splitlines():
        if line.startswith("https://"):
            assert "#" in line, f"missing md5: {line}"


def test_dev_matches_core():
    core, dev = pins("core"), pins("dev")
    assert {k: dev.get(k) for k in core} == core
