#!/usr/bin/env bash
# Fetch the inputs of cameor_reads.pbs into the current directory:
# FULL.fa (Cameor v2, GCA_977071245.1, md5-checked, decompressed) and
# asmqc.sif. Run once on a frontend before qsub; needs curl and singularity.
#
# Usage: get_inputs.sh [asmqc image version, default 0.1.0rc1]
set -euo pipefail
VERSION=${1:-0.1.0rc1}
BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/977/071/245/GCA_977071245.1_pisum_sativum_cameor_v2
GZ=GCA_977071245.1_pisum_sativum_cameor_v2_genomic.fna.gz

if [[ ! -s FULL.fa ]]; then
    curl -sf --retry 5 -O "$BASE/md5checksums.txt"
    curl -sf --retry 5 -C - -O "$BASE/$GZ"
    grep " ./$GZ\$" md5checksums.txt | sed 's| \./| |' | md5sum -c
    gzip -dc "$GZ" > FULL.fa.part && mv FULL.fa.part FULL.fa && rm "$GZ"
fi
if [[ ! -s asmqc.sif ]]; then
    singularity pull asmqc.sif "oras://ghcr.io/peaterpan-pan-genome-project/asmqc/sif:$VERSION"
fi
ls -l FULL.fa asmqc.sif
