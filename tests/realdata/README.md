# Real-data fixture: Cameor v2 subset

A ~50 Mb fixture cut from the published *Pisum sativum* 'Caméor' v2 assembly
(GCA_977071245.1, ENA PRJEB78861; Kreplak et al., Scientific Data) with the
matching subsets of its public HiFi and ONT reads and of the Illumina reads
of the Caméor v1 project (PRJEB30482, same genotype). It fills the gap between
the synthetic smoke-test data (exact truth, no real repeats or genes) and
full-genome runs (hours to days).

The scripts are committed; the data are not (they live under the untracked
`validation/`).

## Contents

* **Chromosomes**: `chr1`–`chr7` (source records `CDSBUU010000001.1` …
  `CDSBUU010000007.1`), each built from both 3 Mb ends (telomeres) and
  blocks around rDNA arrays and interstitial telomeric arrays found by asmqc
  on the full assembly, joined by 100 bp N gaps.
* **Unplaced scaffolds**: organelle-like, rDNA-only and telomere-carrying
  scaffolds, duplicates whose source lies inside the fixture, then a
  deterministic sample, up to 5 Mb.
* **AGP**: components are named `<source>_subseq_<start>:<end>` with source
  coordinates; the joins are `contig` gaps without linkage.

The joins are artificial. No read spans them, so M9 (CRAQ) may report
breakpoints there; they count as CSEs near AGP junctions.

## Building the fixture

1. **Source assembly**: download `GCA_977071245.1_pisum_sativum_cameor_v2_genomic.fna.gz`
   from the NCBI FTP site, check its md5 against `md5checksums.txt` and
   decompress it to `FULL.fa`.
2. **Full-assembly run**: run asmqc modules 1, 4, 5 and 7 on it:

   ```
   M=CDSBUU010000001.1=chr1,CDSBUU010000002.1=chr2,CDSBUU010000003.1=chr3,CDSBUU010000004.1=chr4,CDSBUU010000005.1=chr5,CDSBUU010000006.1=chr6,CDSBUU010000007.1=chr7
   singularity run asmqc.sif run --assembly FULL.fa --label cameor_v2 \
       --outdir full_run --chromosomes "$M" --modules 1,4,5,7 --threads 16 --mem-gb 64
   ```
3. **Fixture**: cut it from `FULL.fa` with that result:

   ```
   python tests/realdata/make_cameor_fixture.py --fasta FULL.fa \
       --result full_run/cameor_v2 --chromosomes "$M" --out fixture
   ```

   The output is `fixture/fixture.fa`, `fixture.agp`, `regions.bed`,
   `chromosomes.txt` and `fixture.json` (blocks, reasons, totals).

## Extracting the reads (MetaCentrum)

The reads are mapped to the full assembly, not the fixture, so that reads of
repeats elsewhere in the genome do not pile up on it. `extract_reads.sh`
downloads one ENA run at a time, checks its md5, maps it with mm2-plus from
the asmqc image, keeps reads whose primary alignment overlaps `regions.bed`
and deletes the run.

Put `fixture/regions.bed`, `extract_reads.sh`, `cameor_reads.pbs` and
`get_inputs.sh` in one directory on storage. `get_inputs.sh` downloads and
checks `FULL.fa` and pulls `asmqc.sif` there (run it once on a frontend).
Then submit from that directory:

```
qsub -v READS=illumina cameor_reads.pbs
qsub -v READS=hifi     cameor_reads.pbs
qsub -v READS=ont      cameor_reads.pbs
```

* `illumina`: three 2×150 bp paired-end runs of Caméor v1, ERR3014756,
  ERR3014760, ERR3014765 (365–477 bp inserts, ~120 Gb, about 31×). A pair is
  kept when either mate maps primarily into the regions.
* `hifi`: all four HiFi runs, ERR9972527–ERR9972530 (73.1 Gb, about 19× of
  the genome).
* `ont`: the six largest ONT R9.4.1 runs (about 130 Gb, about 33×); all 20
  runs are ERR9980778–ERR9980797.

Results: `reads_<type>/subset.fq.gz` (Illumina: `subset_R1.fq.gz`,
`subset_R2.fq.gz`) and `subset.stats.tsv`. The job
requests 32 CPUs, 96 GB and 200 GB scratch for 48 h; these are first
estimates.

## Using the fixture

All three read types went into Cameor v2 (ONT for assembly, HiFi and short
reads for polishing), so runs declare `--reads-used-in-assembly yes`:

```
singularity run asmqc.sif run --assembly fixture/fixture.fa --agp fixture/fixture.agp \
    --chromosomes "$(cat fixture/chromosomes.txt)" --label cameor_fixture \
    --illumina reads_illumina/subset_R1.fq.gz,reads_illumina/subset_R2.fq.gz \
    --hifi reads_hifi/subset.fq.gz --ont reads_ont/subset.fq.gz --ont-chemistry r9 \
    --reads-used-in-assembly yes --outdir fixture_run
```

M8 and M11 use Illumina; M9 uses HiFi with Illumina as `-ngs`. Leave out
`--illumina` to run M8 and M11 on HiFi, or `--hifi` to run M9 on ONT. HiFi
coverage is about 19×, below the 20× `low_coverage` threshold, so M9 (and M8
on HiFi) raise that warning.

There is no ground truth. Uses:

* regression: a run on the same image must reproduce `qc_summary.tsv`
  (except `run_date`);
* plausibility: telomere arm statuses and rDNA loci agree with the full run;
  BUSCO completeness is proportional to the gene-rich share;
* resources: per-module times and peak memory to extrapolate SPEC §9.
