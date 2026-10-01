"""Command-line entry point (SPEC §4.1).

Subcommands are registered here; their implementations land in later
milestones (IMPLEMENTATION_PLAN §3).
"""

import argparse
import sys

from asmqc import __version__


def _not_implemented(args: argparse.Namespace) -> int:
    print(f"asmqc {args.command}: not implemented yet", file=sys.stderr)
    return 3


def _version(args: argparse.Namespace) -> int:
    print(f"asmqc {__version__}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="asmqc", description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("run", help="run QC on one assembly").set_defaults(func=_not_implemented)
    sub.add_parser("version", help="print versions").set_defaults(func=_version)
    sub.add_parser("aggregate", help="merge result directories").set_defaults(
        func=_not_implemented
    )
    sub.add_parser("test", help="run the built-in smoke test").set_defaults(
        func=_not_implemented
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
