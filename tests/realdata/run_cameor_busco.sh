#!/usr/bin/env bash
# BUSCO on the full Cameor v2 assembly, for the M6 synteny anchors.
#
# Runs asmqc module 6 only, so BUSCO version, lineage (fabales_odb12.2) and
# parameters are those of M6 on any tested assembly, and sequences are
# renamed to chr1-chr7 by the same --chromosomes mapping. The anchor table
# is then built from results/cameor_v2/m06_busco/full_table.tsv.
#
# Expects in the directory of this script: asmqc.sif, FULL.fa, chromosomes.txt
# (symlinks are fine; the parent directory is bound into the container).
#
# Usage: [THREADS=64] [MEM_GB=128] ./run_cameor_busco.sh
# Run detached: nohup ./run_cameor_busco.sh > run.log 2>&1 &
set -euo pipefail

DIR=$(cd "$(dirname "$0")" && pwd)
cd "$DIR"
THREADS=${THREADS:-64}
MEM_GB=${MEM_GB:-128}
SING=$(command -v singularity || command -v apptainer) || { echo "no singularity/apptainer" >&2; exit 1; }
for f in asmqc.sif FULL.fa chromosomes.txt; do
    [[ -s $f ]] || { echo "missing $DIR/$f" >&2; exit 1; }
done
BIND=$(dirname "$DIR")   # covers symlink targets in sibling directories
export TMPDIR="$DIR/tmp"
mkdir -p "$TMPDIR"

echo "$(date +%F\ %T) BUSCO on Cameor v2, $THREADS threads"
"$SING" run -B "$BIND" asmqc.sif run --assembly "$DIR/FULL.fa" --label cameor_v2 \
    --chromosomes "$(cat chromosomes.txt)" --modules 6 \
    --threads "$THREADS" --mem-gb "$MEM_GB" \
    --outdir "$DIR/results" --workdir "$DIR/work" --keep-intermediates
echo "$(date +%F\ %T) done: $DIR/results/cameor_v2/m06_busco/full_table.tsv"
