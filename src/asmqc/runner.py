"""`asmqc run`: validate, plan, run Snakemake, merge, clean up (SPEC §4, §5.2)."""

import datetime as dt
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

import yaml

from asmqc import __version__, summary, tools
from asmqc import flags as fl
from asmqc import plan as planner
from asmqc.params import PARAMS
from asmqc.seqio import file_md5
from asmqc.validate import RunOptions, ValidationError, cpu_flags, validate

EXIT_OK, EXIT_INVALID, EXIT_MODULE_FAILED = 0, 1, 2

# Top-level entries the workflow creates in the workdir. Cleanup removes only
# these, so a --workdir that holds other files is never wiped.
WORK_ENTRIES = ("config.yaml", "chromosome_map.json", "prep", "scan", "map", "meryl",
                "benchmarks", "tmp", ".cache", ".snakemake", *(f"m{n:02d}" for n in range(1, 12)))


def asmqc_home() -> Path:
    """Install root holding workflow/, envs/, refs/ (/opt/asmqc in the image)."""
    env = os.environ.get("ASMQC_HOME")
    return Path(env) if env else Path(__file__).resolve().parents[2]


def refs_dir() -> Path:
    """Reference data (SPEC §6): $ASMQC_REFS, else <ASMQC_HOME>/refs."""
    env = os.environ.get("ASMQC_REFS")
    return Path(env).absolute() if env else asmqc_home() / "refs"


def build_info() -> dict:
    """Version, git commit and lock-file hashes (VERSION.json in the image)."""
    home = asmqc_home()
    vfile = home / "VERSION.json"
    if vfile.exists():
        return json.loads(vfile.read_text())
    try:
        commit = subprocess.run(["git", "-C", str(home), "rev-parse", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    locks = {f.stem: hashlib.sha256(f.read_bytes()).hexdigest()
             for f in sorted((home / "envs").glob("*.lock"))}
    return {"version": __version__, "git_commit": commit, "lockfile_sha256": locks}


def reference_data() -> dict:
    """Reference md5s as verified at build time (refs/refs.tsv), plus lineage facts."""
    table = asmqc_home() / "refs" / "refs.tsv"
    out: dict[str, dict] = {}
    if table.exists():
        lines = table.read_text().splitlines()
        for line in lines[1:]:
            name, _url, md5, *_ = line.split("\t")
            out[name] = {"md5": md5}
    lineage = PARAMS["m06"]["lineage"]
    if lineage in out:
        out[lineage] |= {"date": PARAMS["m06"]["lineage_date"],
                         "markers": PARAMS["m06"]["n_markers"]}
    return out


def host_info(opts: RunOptions) -> dict:
    model = "unknown"
    try:
        m = re.search(r"^model name\s*:\s*(.*)$", Path("/proc/cpuinfo").read_text(), re.MULTILINE)
        model = m[1] if m else model
    except OSError:
        pass
    simd = sorted(f for f in cpu_flags() if f.startswith(("avx", "sse4")))
    return {"cpu_model": model, "cpu_flags": simd, "threads": opts.threads,
            "mem_gb": opts.mem_gb}


def inputs_info(opts: RunOptions, assembly_md5: str | None) -> dict:
    def name(p: Path) -> str:
        return str(p.absolute()) if opts.manifest_full_paths else p.name

    info: dict = {"assembly": {"name": name(opts.assembly), "md5": assembly_md5,
                               "bytes": opts.assembly.stat().st_size}}
    if opts.agp:
        info["agp"] = {"name": name(opts.agp), "md5": file_md5(opts.agp)}
    reads = []
    for kind in ("illumina", "hifi", "ont"):
        files = getattr(opts, kind)
        if files:
            r = {"type": kind}
            if kind == "ont":
                r["chemistry"] = opts.ont_chemistry
            r["files"] = [{"name": name(f), "bytes": f.stat().st_size} for f in files]
            r["reads_used_in_assembly"] = opts.reads_used_in_assembly
            reads.append(r)
    info["reads"] = reads
    return info


def write_config(opts: RunOptions, plan: dict[str, planner.ModulePlan]) -> Path:
    work = opts.work.absolute()
    work.mkdir(parents=True, exist_ok=True)
    cfg = {
        "label": opts.label,
        "result_dir": str(opts.result_dir.absolute()),
        "work": str(work),
        "assembly": str(opts.assembly.absolute()),
        "agp": str(opts.agp.absolute()) if opts.agp else None,
        "chromosome_map": opts.chromosomes,
        "reads": {k: [str(p.absolute()) for p in getattr(opts, k)]
                  for k in ("illumina", "hifi", "ont")},
        "reads_used_in_assembly": opts.reads_used_in_assembly,
        "ont_chemistry": opts.ont_chemistry,
        "threads": opts.threads,
        "mem_mb": opts.mem_gb * 1024,
        "plan": planner.as_dict(plan),
        "read_types": planner.read_types_needed(plan),
        "python": sys.executable,
        "refs": str(refs_dir()),
    }
    path = work / "config.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    (work / "chromosome_map.json").write_text(json.dumps(opts.chromosomes))
    return path


def run_snakemake(opts: RunOptions, config: Path) -> int:
    logs = opts.result_dir / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-m", "snakemake",
           "--snakefile", str(asmqc_home() / "workflow" / "Snakefile"),
           "--configfile", str(config),
           "--directory", str(opts.work.absolute()),
           "--cores", str(opts.threads),
           "--resources", f"mem_mb={opts.mem_gb * 1024}",
           "--rerun-incomplete", "--keep-going", "--printshellcmds"]
    with (logs / "snakemake.log").open("a") as log:
        return subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False,
                              env=os.environ | work_env(opts.work.absolute())).returncode


def work_env(work: Path) -> dict[str, str]:
    """Caches and temporary files under the workdir, never $HOME or /tmp (SPEC §5.2).

    $HOME can be read-only (image %test, some cluster nodes); Snakemake's
    source cache and matplotlib would otherwise write there.
    """
    cache, tmp = work / ".cache", work / "tmp"
    cache.mkdir(parents=True, exist_ok=True)
    tmp.mkdir(parents=True, exist_ok=True)
    return {"XDG_CACHE_HOME": str(cache), "MPLCONFIGDIR": str(cache / "matplotlib"),
            "TMPDIR": str(tmp)}


def run(opts: RunOptions, command_line: str) -> int:
    start = dt.datetime.now(dt.UTC)
    try:
        warnings = validate(opts)
        plan = planner.make_plan(opts)
    except ValidationError as e:
        for p in e.problems:
            print(f"asmqc: error: {p}", file=sys.stderr)
        return EXIT_INVALID
    except ValueError as e:
        print(f"asmqc: error: {e}", file=sys.stderr)
        return EXIT_INVALID
    for w in warnings:
        print(f"asmqc: warning: {w}", file=sys.stderr)

    print(f"asmqc {build_info()['version']}: {opts.label}\n{planner.describe(plan)}",
          file=sys.stderr)
    if opts.dry_run:
        return EXIT_OK

    opts.result_dir.mkdir(parents=True, exist_ok=True)
    config = write_config(opts, plan)
    os.environ.update(work_env(opts.work.absolute()))  # also for the report (matplotlib)
    t0 = time.monotonic()
    rc = run_snakemake(opts, config)
    work = opts.work.absolute()

    results = summary.collect(work, planner.as_dict(plan))
    md5_file = work / "prep" / "asm.md5"
    assembly_md5 = md5_file.read_text().strip() if md5_file.exists() else None
    info = build_info()
    identity = {"label": opts.label, "asmqc_version": info["version"],
                "assembly_md5": assembly_md5, "run_date": start.date().isoformat()}
    out = opts.result_dir
    summary.write_summary(out / "qc_summary.tsv", summary.summary_row(identity, results))
    fl.write(out / "flags.tsv", summary.all_flags(results))

    end = dt.datetime.now(dt.UTC)
    manifest = {
        "asmqc_version": info["version"],
        "git_commit": info["git_commit"],
        "lockfile_sha256": info["lockfile_sha256"],
        "label": opts.label,
        "command_line": command_line,
        "start": start.isoformat(timespec="seconds"),
        "end": end.isoformat(timespec="seconds"),
        "wall_seconds": summary.wall_seconds(work) | {"total": round(time.monotonic() - t0, 1)},
        "snakemake_exit_code": rc,
        "host": host_info(opts),
        "inputs": inputs_info(opts, assembly_md5),
        "chromosome_map": opts.chromosomes,
        "tools": tools.all_versions(),
        "reference_data": reference_data(),
        "parameters": PARAMS,
        "modules": {m: {"status": r.status, "read_type": plan[m].read_type,
                        "reason": r.reason} for m, r in results.items()},
    }
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    scrub_paths(out / "logs", path_placeholders(opts))
    try:
        from asmqc import report

        report.render_assembly(out, work)
    except Exception as e:  # noqa: BLE001 - the report must not hide module results
        print(f"asmqc: error: report.html not written: {e}", file=sys.stderr)
        return EXIT_MODULE_FAILED
    failed = [m for m, r in results.items() if r.status == "failed"]
    if failed:
        print(f"asmqc: module(s) failed: {', '.join(failed)}; see {out / 'logs'}",
              file=sys.stderr)
        return EXIT_MODULE_FAILED
    if not opts.keep_intermediates:
        clean_work(work)
    return EXIT_OK


def path_placeholders(opts: RunOptions) -> list[tuple[str, str]]:
    """Local directories and their placeholders, longest path first."""
    pairs = {
        str(opts.work.absolute()): "<workdir>",
        str(opts.outdir.absolute()): "<outdir>",
        str(refs_dir()): "<refs>",
        str(asmqc_home()): "<asmqc>",
        sys.prefix: "<python>",
    }
    if os.environ.get("ASMQC_ENV_ROOT"):
        pairs[os.environ["ASMQC_ENV_ROOT"]] = "<envs>"
    inputs = [opts.assembly, opts.agp, *opts.illumina, *opts.hifi, *opts.ont]
    for f in inputs:
        if f is not None:
            pairs.setdefault(str(f.absolute().parent), "<input>")
    for d in (os.environ.get("TMPDIR"), os.environ.get("HOME")):
        if d and len(d) > 1:
            pairs.setdefault(d, "<tmp>" if d == os.environ.get("TMPDIR") else "<home>")
    return sorted(pairs.items(), key=lambda kv: -len(kv[0]))


def scrub_paths(logs: Path, placeholders: list[tuple[str, str]]) -> None:
    """Replace local paths in the logs: <label>/ is shared with the consortium."""
    if not logs.is_dir():
        return
    for f in logs.iterdir():
        if not f.is_file():
            continue
        text = f.read_text(errors="replace")
        for path, ph in placeholders:
            text = text.replace(path, ph)
        f.write_text(text)


def clean_work(work: Path) -> None:
    for name in WORK_ENTRIES:
        p = work / name
        if p.is_dir() and not p.is_symlink():
            shutil.rmtree(p)
        elif p.exists() or p.is_symlink():
            p.unlink()
    try:
        work.rmdir()
    except OSError:
        pass  # not empty: holds files that are not ours


PATH_OPTIONS = ("--assembly", "--agp", "--illumina", "--hifi", "--ont", "--outdir",
                "--workdir")


def command_line(argv: list[str], full_paths: bool) -> str:
    """Quoted command line; path arguments reduced to file names (SPEC §7.3)."""
    def names(value: str) -> str:
        return ",".join(Path(p).name for p in value.split(","))

    out, path_next = [], False
    for a in argv:
        opt, eq, value = a.partition("=")
        if full_paths:
            out.append(a)
        elif path_next:
            out.append(names(a))
        elif eq and opt in PATH_OPTIONS:
            out.append(f"{opt}={names(value)}")
        else:
            out.append(a)
        path_next = a in PATH_OPTIONS
    return " ".join(shlex.quote(a) for a in ["asmqc", *out])
