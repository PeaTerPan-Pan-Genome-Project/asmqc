"""Input validation (SPEC §4.2). Every failure is collected and reported at once."""

import os
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from asmqc.seqio import fasta_names

CHROMOSOMES = [f"chr{i}" for i in range(1, 8)]
LABEL_RE = re.compile(r"[A-Za-z0-9._-]+")


@dataclass
class RunOptions:
    assembly: Path
    label: str
    outdir: Path
    agp: Path | None = None
    chromosomes: dict[str, str] = field(default_factory=dict)
    illumina: list[Path] = field(default_factory=list)
    hifi: list[Path] = field(default_factory=list)
    ont: list[Path] = field(default_factory=list)
    ont_chemistry: str | None = None
    reads_used_in_assembly: str | None = None
    threads: int = 32
    mem_gb: int = 128
    workdir: Path | None = None
    keep_intermediates: bool = False
    modules: list[str] | None = None
    manifest_full_paths: bool = False
    dry_run: bool = False

    @property
    def result_dir(self) -> Path:
        return self.outdir / self.label

    @property
    def work(self) -> Path:
        return self.workdir if self.workdir else self.result_dir / "work"


class ValidationError(Exception):
    def __init__(self, problems: list[str]):
        super().__init__("\n".join(problems))
        self.problems = problems


def parse_chromosome_map(text: str | None) -> dict[str, str]:
    """Parse 'old1=chr1,old2=chr2'. Raises ValueError on malformed input."""
    if not text:
        return {}
    mapping: dict[str, str] = {}
    for item in text.split(","):
        old, sep, new = item.strip().partition("=")
        if not sep or not old or not new:
            raise ValueError(f"--chromosomes: malformed entry {item!r}; expected old=new")
        if old in mapping:
            raise ValueError(f"--chromosomes: {old!r} mapped twice")
        mapping[old] = new
    return mapping


def cpu_flags() -> set[str]:
    try:
        text = Path("/proc/cpuinfo").read_text()
    except OSError:
        return set()
    m = re.search(r"^flags\s*:\s*(.*)$", text, re.MULTILINE)
    return set(m[1].split()) if m else set()


def check_chromosomes(names: list[str], mapping: dict[str, str]) -> list[str]:
    """Problems with chromosome naming after applying the mapping."""
    problems = []
    present = set(names)
    missing_sources = [old for old in mapping if old not in present]
    if missing_sources:
        problems.append("--chromosomes: not in the assembly: " + ", ".join(missing_sources))
    renamed = [mapping.get(n, n) for n in names]
    counts: dict[str, int] = {}
    for n in renamed:
        counts[n] = counts.get(n, 0) + 1
    missing = [c for c in CHROMOSOMES if c not in counts]
    if missing:
        shown = ", ".join(names[:20]) + (f", … ({len(names)} total)" if len(names) > 20 else "")
        problems.append(f"assembly lacks {', '.join(missing)} after --chromosomes mapping; "
                        f"names found: {shown}")
    clashes = [n for n in mapping.values() if counts.get(n, 0) > 1]
    if clashes:
        problems.append("--chromosomes: renaming creates duplicate names: "
                        + ", ".join(sorted(set(clashes))))
    return problems


def validate(opts: RunOptions) -> list[str]:
    """Raise ValidationError listing every problem; return warnings."""
    problems: list[str] = []
    warnings: list[str] = []

    if not LABEL_RE.fullmatch(opts.label):
        problems.append(f"--label {opts.label!r}: only [A-Za-z0-9._-] allowed")

    for name, path in [("--assembly", opts.assembly), ("--agp", opts.agp)]:
        if path is not None and not (path.is_file() and os.access(path, os.R_OK)):
            problems.append(f"{name} {path}: not a readable file")
    for kind in ("illumina", "hifi", "ont"):
        for path in getattr(opts, kind):
            if not (path.is_file() and os.access(path, os.R_OK)):
                problems.append(f"--{kind} {path}: not a readable file")
    if len(opts.illumina) % 2:
        problems.append("--illumina: files must come in R1,R2 pairs")

    any_reads = opts.illumina or opts.hifi or opts.ont
    if any_reads and opts.reads_used_in_assembly not in ("yes", "no"):
        problems.append("--reads-used-in-assembly yes|no is required when reads are given")
    if opts.ont and opts.ont_chemistry not in ("r9", "r10"):
        problems.append("--ont-chemistry r9|r10 is required with --ont")

    if opts.threads < 1 or opts.mem_gb < 1:
        problems.append("--threads and --mem-gb must be positive")

    if opts.assembly.is_file() and not problems:
        try:
            names = list(fasta_names(opts.assembly))
        except (OSError, EOFError, UnicodeDecodeError) as e:
            problems.append(f"--assembly {opts.assembly}: cannot read FASTA headers: {e}")
        else:
            if not names:
                problems.append(f"--assembly {opts.assembly}: no FASTA records")
            else:
                problems += check_chromosomes(names, opts.chromosomes)

    out = opts.outdir
    probe = out if out.exists() else _existing_parent(out)
    if not (probe.is_dir() and os.access(probe, os.W_OK)):
        problems.append(f"--outdir {out}: not writable")

    if "avx2" not in cpu_flags():
        problems.append("CPU lacks AVX2, required by mm2-plus")

    if problems:
        raise ValidationError(problems)

    need = estimate_work_bytes(opts)
    free = shutil.disk_usage(_existing_parent(opts.work)).free
    if free < need:
        warnings.append(f"workdir has {free / 1e9:.0f} GB free; this run may need "
                        f"~{need / 1e9:.0f} GB (SPEC §9)")
    return warnings


def estimate_work_bytes(opts: RunOptions) -> int:
    """Rough workdir need from SPEC §9, scaled to the assembly size (4.3 Gb basis)."""
    scale = max(opts.assembly.stat().st_size / 4.3e9, 0.01)
    gb = 60
    if opts.illumina:
        gb += 150 + 100
    if opts.hifi or opts.ont:
        gb += 200 + 50
    return int(gb * scale * 1e9)


def _existing_parent(path: Path) -> Path:
    path = path.absolute()
    while not path.exists():
        path = path.parent
    return path
