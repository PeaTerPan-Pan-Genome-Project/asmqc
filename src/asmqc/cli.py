"""asmqc: comparable QC of chromosome-level pea genome assemblies (SPEC §4.1)."""

import argparse
import json
import sys
from pathlib import Path

from asmqc import runner, tools
from asmqc.validate import RunOptions, parse_chromosome_map


def _paths(text: str) -> list[Path]:
    return [Path(p) for p in text.split(",") if p]


def _not_implemented(args: argparse.Namespace) -> int:
    print(f"asmqc {args.command}: not implemented yet", file=sys.stderr)
    return 3


def _version(args: argparse.Namespace) -> int:
    info = runner.build_info()
    out = {
        "asmqc_version": info["version"],
        "git_commit": info["git_commit"],
        "lockfile_sha256": info["lockfile_sha256"],
        "tools": tools.all_versions(),
        "reference_data": runner.reference_data(),
    }
    print(json.dumps(out, indent=2))
    return 0


def _aggregate(args: argparse.Namespace) -> int:
    from asmqc import aggregate

    return aggregate.run(args.dirs, args.out)


def _run(args: argparse.Namespace, argv: list[str]) -> int:
    try:
        chromosomes = parse_chromosome_map(args.chromosomes)
    except ValueError as e:
        print(f"asmqc: error: {e}", file=sys.stderr)
        return runner.EXIT_INVALID
    opts = RunOptions(
        assembly=args.assembly, label=args.label, outdir=args.outdir, agp=args.agp,
        chromosomes=chromosomes, illumina=args.illumina, hifi=args.hifi, ont=args.ont,
        ont_chemistry=args.ont_chemistry, reads_used_in_assembly=args.reads_used_in_assembly,
        threads=args.threads, mem_gb=args.mem_gb, workdir=args.workdir,
        keep_intermediates=args.keep_intermediates,
        modules=args.modules.split(",") if args.modules else None,
        manifest_full_paths=args.manifest_full_paths, dry_run=args.dry_run)
    return runner.run(opts, runner.command_line(argv, args.manifest_full_paths))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="asmqc", description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    r = sub.add_parser("run", help="run QC on one assembly")
    r.add_argument("--assembly", type=Path, required=True, help="FASTA, gzip allowed")
    r.add_argument("--label", required=True, help="[A-Za-z0-9._-]+; names all outputs")
    r.add_argument("--outdir", type=Path, required=True, help="results go to OUTDIR/LABEL/")
    r.add_argument("--agp", type=Path)
    r.add_argument("--chromosomes", metavar="OLD=chrN,...",
                   help="map assembly names to chr1..chr7 (FASTA and AGP)")
    r.add_argument("--illumina", type=_paths, default=[], metavar="R1,R2[,R1b,R2b...]")
    r.add_argument("--hifi", type=_paths, default=[], metavar="FILE[,FILE...]")
    r.add_argument("--ont", type=_paths, default=[], metavar="FILE[,FILE...]")
    r.add_argument("--ont-chemistry", choices=("r9", "r10"))
    r.add_argument("--reads-used-in-assembly", choices=("yes", "no"),
                   help="required when any reads are given")
    r.add_argument("--threads", type=int, default=32)
    r.add_argument("--mem-gb", type=int, default=128)
    r.add_argument("--workdir", type=Path, help="intermediates; default OUTDIR/LABEL/work")
    r.add_argument("--keep-intermediates", action="store_true")
    r.add_argument("--modules", help="comma list, e.g. 1,2,4; default: all applicable")
    r.add_argument("--manifest-full-paths", action="store_true",
                   help="record full input paths in run_manifest.json, not only names")
    r.add_argument("--dry-run", action="store_true",
                   help="validate inputs, print the plan, run nothing")
    r.set_defaults(func=_run)

    sub.add_parser("version", help="print versions as JSON").set_defaults(func=_version)
    a = sub.add_parser("aggregate", help="merge result directories into one table and report")
    a.add_argument("dirs", type=Path, nargs="+", metavar="DIR", help="<outdir>/<label>/ dirs")
    a.add_argument("--out", type=Path, required=True)
    a.set_defaults(func=_aggregate)
    sub.add_parser("test", help="run the built-in smoke test").set_defaults(
        func=_not_implemented)
    return parser


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    args = build_parser().parse_args(argv)
    if args.func is _run:
        return _run(args, argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
