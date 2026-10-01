#!/usr/bin/env bash
# Extract the reads of a real-data fixture from public ENA runs.
#
# Each run is downloaded, md5-checked, mapped with mm2-plus to the FULL source
# assembly, and only reads whose primary alignment overlaps the fixture
# regions (regions.bed from make_cameor_fixture.py) are kept. Mapping to the
# full assembly, not the fixture, keeps reads of repeats from elsewhere in the
# genome off the fixture. One run is on disk at a time.
#
# Tools (mm2plus, samtools, seqkit) come from the asmqc image, so the subset
# is reproducible with the pinned versions. Needs curl and singularity or
# apptainer on the host.
#
# Usage: extract_reads.sh SIF ASSEMBLY REGIONS_BED PRESET OUTDIR THREADS RUN [RUN ...]
#   PRESET  sr (Illumina pairs), map-hifi (PacBio HiFi) or map-ont (ONT R9.4.1)
#   RUN     ENA run accession, e.g. ERR9972527
# Writes per run OUTDIR/<RUN>.subset.fq.gz (or <RUN>_1/_2 for pairs), all runs
# together as OUTDIR/subset.fq.gz (or subset_R1/_R2.fq.gz), and
# OUTDIR/subset.stats.tsv. A pair is kept when either mate's primary
# alignment overlaps the regions, so R1 and R2 stay in step.
set -euo pipefail

[[ $# -ge 7 ]] || { sed -n '2,20p' "$0"; exit 1; }
SIF=$(realpath "$1"); ASM=$(realpath "$2"); REGIONS=$(realpath "$3")
PRESET=$4; OUT=$(realpath -m "$5"); THREADS=$6
shift 6
case $PRESET in sr|map-hifi|map-ont) ;; *) echo "PRESET must be sr, map-hifi or map-ont" >&2; exit 1;; esac

mkdir -p "$OUT"
WORK="$OUT/work"
mkdir -p "$WORK"
SING=$(command -v singularity || command -v apptainer)
binds="$(dirname "$ASM"),$(dirname "$REGIONS"),$OUT"
tool() {
    "$SING" exec -B "$binds" "$SIF" env PATH=/opt/envs/core/bin:/usr/bin:/bin "$@"
}

# One index for all runs; -I larger than the genome keeps it in one part,
# so primary alignments are not split across index parts.
IDX="$WORK/asm.$PRESET.mmi"
if [[ ! -s $IDX ]]; then
    tool mm2plus -x "$PRESET" -I 16G -t "$THREADS" -d "$IDX" "$ASM"
fi

for RUN in "$@"; do
    if [[ -s "$OUT/$RUN.subset.fq.gz" || -s "$OUT/${RUN}_2.subset.fq.gz" ]]; then
        echo "$RUN: done before, skipped"; continue
    fi
    meta=$(curl -sf "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$RUN&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" | tail -n 1)
    IFS=';' read -r -a urls <<<"$(cut -f2 <<<"$meta")"
    IFS=';' read -r -a md5s <<<"$(cut -f3 <<<"$meta")"
    [[ ${#urls[@]} -gt 0 && -n ${urls[0]} ]] || { echo "$RUN: no FASTQ at ENA" >&2; exit 1; }
    if [[ $PRESET == sr ]]; then
        [[ ${#urls[@]} -eq 2 ]] || { echo "$RUN: expected 2 FASTQ files for pairs" >&2; exit 1; }
        names=(_1 _2)
    else
        urls=("${urls[0]}"); md5s=("${md5s[0]}"); names=("")
    fi
    fqs=()
    for i in "${!urls[@]}"; do
        fq="$WORK/$RUN${names[$i]}.fastq.gz"
        echo "$RUN: downloading ${urls[$i]}"
        curl -sf --retry 5 -C - -o "$fq" "https://${urls[$i]}"
        echo "${md5s[$i]}  $fq" | md5sum -c --quiet
        fqs+=("$fq")
    done

    echo "$RUN: mapping"
    tool bash -o pipefail -c "mm2plus -ax $PRESET --secondary=no -t $THREADS '$IDX' ${fqs[*]} \
        | samtools view -F 0x904 -L '$REGIONS' - | cut -f1 | sort -u > '$WORK/$RUN.names'"
    for i in "${!fqs[@]}"; do
        sub="$OUT/$RUN${names[$i]}.subset.fq.gz"
        tool seqkit grep -j "$THREADS" -f "$WORK/$RUN.names" "${fqs[$i]}" -o "$sub.part"
        mv "$sub.part" "$sub"
    done
    echo "$RUN: kept $(wc -l < "$WORK/$RUN.names") reads or pairs"
    rm -f "${fqs[@]}"
done

if [[ $PRESET == sr ]]; then
    cat "$OUT"/ERR*_1.subset.fq.gz > "$OUT/subset_R1.fq.gz"
    cat "$OUT"/ERR*_2.subset.fq.gz > "$OUT/subset_R2.fq.gz"
    tool seqkit stats -T -a "$OUT/subset_R1.fq.gz" "$OUT/subset_R2.fq.gz" > "$OUT/subset.stats.tsv"
else
    cat "$OUT"/ERR*.subset.fq.gz > "$OUT/subset.fq.gz"
    tool seqkit stats -T -a "$OUT/subset.fq.gz" > "$OUT/subset.stats.tsv"
fi
cat "$OUT/subset.stats.tsv"
