#!/usr/bin/env bash
# Run the Cameor v2 fixture on one large host, without a scheduler.
#
# Stage "reads" extracts the Illumina, HiFi and ONT read subsets in parallel
# (extract_reads.sh, THREADS/3 threads each; mm2plus peaks at ~18 GB each).
# Stage "fixture" runs asmqc on the fixture with all three read types, then
# M9 once more on ONT alone (the ont_r9 path). Every stage resumes: finished
# runs and finished asmqc rules are skipped on a rerun.
#
# Expects in the directory of this script: asmqc.sif, FULL.fa,
# fixture/{fixture.fa,fixture.agp,regions.bed,chromosomes.txt}, extract_reads.sh.
#
# Usage: [THREADS=96] [MEM_GB=200] ./run_cameor_fixture.sh [reads|fixture|all]
# Run detached, e.g.: nohup ./run_cameor_fixture.sh all > run.log 2>&1 &
set -euo pipefail

DIR=$(cd "$(dirname "$0")" && pwd)
cd "$DIR"
STAGE=${1:-all}
THREADS=${THREADS:-96}
MEM_GB=${MEM_GB:-200}
PER=$((THREADS / 3))
SING=$(command -v singularity || command -v apptainer) || { echo "no singularity/apptainer" >&2; exit 1; }
export TMPDIR="$DIR/tmp"
mkdir -p "$TMPDIR" logs

for f in asmqc.sif FULL.fa fixture/fixture.fa fixture/fixture.agp fixture/regions.bed \
         fixture/chromosomes.txt extract_reads.sh; do
    [[ -s $f ]] || { echo "missing $DIR/$f" >&2; exit 1; }
done

reads() {
    local pids=() names=() rc=0
    declare -A RUNS=(
        [illumina]="sr ERR3014756 ERR3014760 ERR3014765"
        [hifi]="map-hifi ERR9972527 ERR9972528 ERR9972529 ERR9972530"
        [ont]="map-ont ERR9980780 ERR9980790 ERR9980792 ERR9980791 ERR9980797 ERR9980793"
    )
    for t in illumina hifi ont; do
        read -r preset runs <<<"${RUNS[$t]}"
        echo "$(date +%F\ %T) reads: $t started ($PER threads), log logs/reads_$t.log"
        # shellcheck disable=SC2086
        bash extract_reads.sh asmqc.sif FULL.fa fixture/regions.bed "$preset" \
            "reads_$t" "$PER" $runs > "logs/reads_$t.log" 2>&1 &
        pids+=($!); names+=("$t")
    done
    for i in "${!pids[@]}"; do
        if wait "${pids[$i]}"; then
            echo "$(date +%F\ %T) reads: ${names[$i]} done"
        else
            echo "$(date +%F\ %T) reads: ${names[$i]} FAILED, see logs/reads_${names[$i]}.log" >&2
            rc=1
        fi
    done
    return $rc
}

fixture() {
    local chrom common
    chrom=$(cat fixture/chromosomes.txt)
    common=(--assembly "$DIR/fixture/fixture.fa" --agp "$DIR/fixture/fixture.agp"
            --chromosomes "$chrom" --reads-used-in-assembly yes
            --threads "$THREADS" --mem-gb "$MEM_GB" --outdir "$DIR/results")
    echo "$(date +%F\ %T) fixture: all modules, Illumina + HiFi + ONT"
    "$SING" run -B "$DIR" asmqc.sif run "${common[@]}" --label cameor_fixture \
        --workdir "$DIR/work/cameor_fixture" --keep-intermediates \
        --illumina "$DIR/reads_illumina/subset_R1.fq.gz,$DIR/reads_illumina/subset_R2.fq.gz" \
        --hifi "$DIR/reads_hifi/subset.fq.gz" \
        --ont "$DIR/reads_ont/subset.fq.gz" --ont-chemistry r9 \
        > logs/fixture_all.log 2>&1 || echo "fixture: asmqc exited $? (see logs/fixture_all.log)"
    echo "$(date +%F\ %T) fixture: M9 on ONT R9.4.1 only"
    "$SING" run -B "$DIR" asmqc.sif run "${common[@]}" --label cameor_fixture_ont \
        --workdir "$DIR/work/cameor_fixture_ont" --keep-intermediates --modules 9 \
        --ont "$DIR/reads_ont/subset.fq.gz" --ont-chemistry r9 \
        > logs/fixture_ont.log 2>&1 || echo "fixture: asmqc exited $? (see logs/fixture_ont.log)"
    echo "$(date +%F\ %T) fixture: results in $DIR/results"
}

case $STAGE in
    reads) reads ;;
    fixture) fixture ;;
    all) reads && fixture ;;
    *) echo "stage must be reads, fixture or all" >&2; exit 1 ;;
esac
