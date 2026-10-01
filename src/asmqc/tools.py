"""Tool environments and run-time version discovery (SPEC §5.1, §7.3).

ASMQC_ENV_ROOT=/opt/envs in the image: each tool runs from
$ASMQC_ENV_ROOT/<env>/bin. Unset (development), every tool comes from the
current environment (envs/dev.yaml) and PATH.
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

ENVS = ("core", "quast", "busco", "merqury", "craq")

# tool -> (env, command, regex on stdout+stderr); command None = read conda-meta
TOOLS: dict[str, tuple[str, list[str] | None, str]] = {
    "snakemake": ("core", ["snakemake", "--version"], r"^(\d+\.\d+\.\d+)"),
    "seqkit": ("core", ["seqkit", "version"], r"seqkit v(\S+)"),
    "samtools": ("core", ["samtools", "--version"], r"^samtools (\S+)"),
    "bcftools": ("core", ["bcftools", "--version"], r"^bcftools (\S+)"),
    "bedtools": ("core", ["bedtools", "--version"], r"bedtools v(\S+)"),
    # mm2plus is a dispatcher that also reports which binary it launches
    "mm2plus": ("core", ["mm2plus", "--version"], r"^(\d+\.\d+\S*)$"),
    "blast": ("core", ["blastn", "-version"], r"blastn: (\S+?)\+"),
    "tidk": ("core", ["tidk", "--version"], r"tidk (\S+)"),
    "quast": ("quast", ["quast.py", "--version"], r"QUAST v(\S+)"),
    "busco": ("busco", ["busco", "--version"], r"BUSCO (\d\S*)"),
    "miniprot": ("busco", ["miniprot", "--version"], r"^(\S+)"),
    "meryl": ("merqury", ["meryl", "--version"], r"meryl (\S+)"),
    "merqury": ("merqury", None, "merqury"),  # merqury.sh has no version option
    "craq": ("craq", None, "craq"),  # prints 1.0.9-alpha for 1.10 (SPEC §7.3)
}


def env_prefix(env: str) -> Path:
    root = os.environ.get("ASMQC_ENV_ROOT")
    return Path(root) / env if root else Path(sys.prefix)


def env_bin(env: str) -> Path | None:
    """bin/ to put first on PATH for a tool env; None in development."""
    return env_prefix(env) / "bin" if os.environ.get("ASMQC_ENV_ROOT") else None


def shell_prefix(env: str, tmpdir: str) -> str:
    """Shell preamble for a rule: tool env on PATH, TMPDIR under the workdir."""
    b = env_bin(env)
    path = f"export PATH={b}:$PATH; " if b else ""
    extra = ""
    if env == "merqury":  # normally set by conda activation, which rules do not use
        extra = f"export MERQURY={env_prefix('merqury') / 'share' / 'merqury'}; "
    return f"{path}{extra}mkdir -p {tmpdir}; export TMPDIR={tmpdir}; "


def conda_meta_version(env: str, package: str) -> str | None:
    meta = env_prefix(env) / "conda-meta"
    for f in sorted(meta.glob(f"{package}-*.json")):
        data = json.loads(f.read_text())
        if data.get("name") == package:
            return data["version"]
    return None


def tool_version(tool: str) -> str | None:
    env, cmd, pattern = TOOLS[tool]
    if cmd is None:
        return conda_meta_version(env, pattern)
    # The env's bin/ goes first on PATH, as in the rules: scripts such as busco
    # start with "#!/usr/bin/env python" and must find their own interpreter.
    b = env_bin(env)
    environ = os.environ | ({"PATH": f"{b}:{os.environ.get('PATH', '')}"} if b else {})
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, errors="replace", timeout=120,
                           check=False, env=environ)
    except (OSError, subprocess.TimeoutExpired):
        return None
    m = re.search(pattern, r.stdout + r.stderr, re.MULTILINE)
    return m[1] if m else None


def all_versions() -> dict[str, str | None]:
    versions = {"python": sys.version.split()[0]}
    versions.update({t: tool_version(t) for t in TOOLS})
    return versions
