#!/usr/bin/env python3
"""Build the M6 synteny anchor table from an asmqc M6 run on the reference.

Anchors are the Complete (single-copy) BUSCOs on chr1-chr7 of the reference
assembly (Caméor v2, GCA_977071245.1), with their positions and strands.
M6 joins them with the tested assembly's BUSCOs for the synteny dotplot;
the reference sequence itself is not needed at run time.

Usage: make_busco_anchors.py RESULT_DIR FAI OUT.tsv
  RESULT_DIR  asmqc result of the reference (<outdir>/<label>/), module 6 run
              with --chromosomes mapping to chr1-chr7
  FAI         .fai of the reference, sequence names as in the source FASTA
  OUT.tsv     anchor table; header lines carry provenance and chr lengths

The committed table was built by tests/realdata/run_cameor_busco.sh with
asmqc 0.1.0rc2 (BUSCO 6.1.0, fabales_odb12.2).
"""

import json
import sys
from pathlib import Path

CHROMOSOMES = [f"chr{i}" for i in range(1, 8)]


def main(argv: list[str]) -> None:
    if len(argv) != 3:
        sys.exit(__doc__)
    res, fai, out = Path(argv[0]), Path(argv[1]), Path(argv[2])
    manifest = json.loads((res / "run_manifest.json").read_text())
    to_chr = manifest["chromosome_map"]
    lengths = {}
    for line in fai.read_text().splitlines():
        name, length = line.split("\t")[:2]
        if name in to_chr:
            lengths[to_chr[name]] = int(length)
    rows, header = [], None
    for line in (res / "m06_busco" / "full_table.tsv").read_text().splitlines():
        if line.startswith("# Busco id"):
            header = line[2:].split("\t")
        if line.startswith("#") or not line.strip():
            continue
        r = dict(zip(header, line.split("\t")))
        if r["Status"] == "Complete" and r["Sequence"] in CHROMOSOMES:
            s, e = sorted((int(r["Gene Start"]), int(r["Gene End"])))
            rows.append((r["Busco id"], r["Sequence"], s, e, r["Strand"]))
    rows.sort(key=lambda x: (CHROMOSOMES.index(x[1]), x[2], x[0]))
    with out.open("w") as fh:
        fh.write("# M6 synteny anchors: Complete single-copy BUSCOs on chr1-chr7\n")
        fh.write("# reference\tCameor v2\tGCA_977071245.1\n")
        fh.write(f"# asmqc\t{manifest['asmqc_version']}\tBUSCO\t{manifest['tools']['busco']}"
                 f"\tlineage\t{manifest['reference_data'].get('fabales_odb12.2', {}).get('date', '')}"
                 "\tfabales_odb12.2\n")
        for c in CHROMOSOMES:
            fh.write(f"# length\t{c}\t{lengths[c]}\n")
        fh.write("busco_id\tchromosome\tstart\tend\tstrand\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    print(f"{len(rows)} anchors -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1:])
