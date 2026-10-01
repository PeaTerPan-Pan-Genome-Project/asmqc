"""Shared `map_<readtype>` stage (SPEC §5.2): mm2-plus, sort, index.

One mm2plus run per read file (per R1/R2 pair for Illumina), each with its
own read group, then merged into one coordinate-sorted, indexed BAM.

Usage: python -m asmqc.mapping --fa FA --read-type T --label L --threads N
           --mem-mb M --out BAM UNIT [UNIT ...]
       UNIT is one file, or R1,R2 for Illumina.
"""

import argparse
import shlex
import subprocess
from pathlib import Path

PRESETS = {"illumina": "sr", "hifi": "map-hifi", "ont_r10": "lr:hq", "ont_r9": "map-ont"}
PLATFORMS = {"illumina": "ILLUMINA", "hifi": "PACBIO", "ont_r10": "ONT", "ont_r9": "ONT"}


def sh(cmd: str) -> None:
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fa", type=Path, required=True)
    ap.add_argument("--read-type", choices=PRESETS, required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--threads", type=int, required=True)
    ap.add_argument("--mem-mb", type=int, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("units", nargs="+")
    args = ap.parse_args(argv)

    rt = args.read_type
    tmp = args.out.parent / f"{rt}_parts"
    tmp.mkdir(parents=True, exist_ok=True)
    sort_threads = max(1, min(8, args.threads // 4))
    sort_mem = max(256, min(2048, args.mem_mb // (2 * sort_threads)))
    parts = []
    for i, unit in enumerate(args.units, 1):
        files = unit.split(",")
        rg = f"@RG\\tID:{rt}{i}\\tSM:{args.label}\\tPL:{PLATFORMS[rt]}"
        part = tmp / f"part{i}.bam"
        sh(f"mm2plus -t {args.threads} -ax {PRESETS[rt]} --max-chain-skip=1000000"
           f" -R {shlex.quote(rg)} {shlex.quote(str(args.fa))}"
           f" {' '.join(shlex.quote(f) for f in files)}"
           f" | samtools sort -@ {sort_threads} -m {sort_mem}M -T {tmp}/sort{i}"
           f" -o {part} -")
        parts.append(part)
    if len(parts) == 1:
        parts[0].rename(args.out)
    else:
        sh(f"samtools merge -@ {args.threads} -f -o {args.out} "
           + " ".join(str(p) for p in parts))
        for p in parts:
            p.unlink()
    sh(f"samtools index -@ {args.threads} {args.out}")
    tmp.rmdir()


if __name__ == "__main__":
    main()
