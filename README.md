# asmqc

asmqc runs a fixed set of quality-control checks on chromosome-level genome
assemblies of pea (*Pisum*; *Lathyrus oleraceus* and *L. fulvus*) and writes
one fixed-schema result per assembly. It is one Singularity image: every
group in the pangenome consortium runs identical tools, versions, parameters
and reference data, so the results merge into one comparable table.

asmqc measures; it never modifies an assembly. The design specification is
[`docs/design/SPEC.md`](docs/design/SPEC.md).

## How it works

A Snakemake workflow inside the image switches modules on according to the
inputs supplied. It runs on a single machine; wrap it in your own scheduler.

| Module | Checks | Runs when | Tools |
|---|---|---|---|
| M1 | integrity, format, ENA submission rules, AGP consistency | always | Python, samtools |
| M2 | contiguity, anchoring to `chr1`–`chr7` | always | QUAST, Python |
| M4 | telomeres (TTTAGGG) at chromosome ends | always | tidk |
| M5 | organelle-derived sequence, rDNA arrays | always | mm2-plus, BLAST+ |
| M6 | gene-space completeness, `fabales_odb12.2`; synteny dotplot against Caméor v2 (report only) | always | BUSCO + miniprot |
| M7 | redundancy of unplaced scaffolds | unplaced ≥ 1 kb exist | mm2-plus |
| M8 | QV, k-mer completeness, spectra-cn | Illumina or HiFi given | meryl, Merqury |
| M9 | read-back structural validation | HiFi or ONT given | CRAQ |
| M11 | homopolymer and dinucleotide-repeat length errors | Illumina or HiFi given | bcftools |

**Read mapping**: each read type is mapped once with mm2-plus and shared by
the modules. M8 and M11 use Illumina reads if given, otherwise HiFi; ONT is
never used for them. M9 uses HiFi if given, otherwise ONT, plus Illumina when
given.

**CRAQ**: M9 runs CRAQ 1.10 from a copy in the work directory with patched
drivers and scripts (`workflow/craq_patch/`). They stream data that CRAQ
loads into memory and run its steps per chromosome in parallel; the output is
byte-identical to unpatched CRAQ.

**Determinism**: the same inputs on the same image give identical values in
`qc_summary.tsv` (except `run_date`). Results are comparable only between
runs of the same `MAJOR.MINOR` image version.

## Requirements

* **CPU**: x86-64 with AVX2 (required by mm2-plus); asmqc stops with an
  error otherwise. 64 threads recommended; the default is 32.
* **Memory**: 256 GB (`--mem-gb 256`) is tested for runs with reads; the
  largest single step used 72 GB. 128 GB, the default, is expected to work
  but has not been measured. Assembly-only runs: 64 GB.
* **Disk**: about 400 GB free for `--workdir` with reads, on local disk;
  the run peaked at 333 GB. The result directory stays below 50 MB.
* **Software**: Singularity or Apptainer to run the image; tested with
  Apptainer 1.4. Nothing else is installed on the host.
* **Network**: needed only to pull the image. Runs are offline; all
  reference data are inside the image.

The figures are for one 4.3 Gb pea assembly; see [Resources](#resources).

## Installation

Pull the image and check its sha256 against the value published with the
release:

```
singularity pull asmqc_0.1.5.sif oras://ghcr.io/peaterpan-pan-genome-project/asmqc/sif:0.1.5
sha256sum asmqc_0.1.5.sif
```

Releases are listed on the
[GitHub releases page](https://github.com/PeaTerPan-Pan-Genome-Project/asmqc/releases);
use the same `MAJOR.MINOR` version as the other groups.

To build the image yourself, run `./build.sh` in a clone of this repository
(Apptainer 1.4, network access during the build). The build downloads every
reference file, verifies its md5 and runs the built-in test.

## Usage

```
singularity run [-B <bind paths>] asmqc_<version>.sif run \
    --assembly ASM.fa[.gz] --label SAMPLE01 --outdir results/ [options]
```

Bind every directory that holds inputs or outputs (`-B /data,/scratch`).

**Example 1, assembly only** (M1, M2, M4–M7; about 30 min for 4.3 Gb):

```
singularity run -B /data asmqc_0.1.5.sif run \
    --assembly /data/PS01.fa.gz --agp /data/PS01.agp --label PS01 \
    --outdir /data/qc --threads 32 --mem-gb 64
```

**Example 2, with Illumina reads** (adds M8 and M11; mapping dominates the run time):

```
singularity run -B /data,/scratch asmqc_0.1.5.sif run \
    --assembly /data/PS01.fa.gz --agp /data/PS01.agp --label PS01 \
    --illumina /data/PS01_R1.fq.gz,/data/PS01_R2.fq.gz \
    --reads-used-in-assembly no \
    --outdir /data/qc --workdir /scratch/asmqc_PS01 --threads 32 --mem-gb 128
```

**Example 3, with HiFi and ONT reads** (adds M8, M9 and M11, all from HiFi).
ONT reads are mapped only when no HiFi reads are given; with both, the ONT
files are recorded in the manifest but not used:

```
singularity run -B /data,/scratch asmqc_0.1.5.sif run \
    --assembly /data/PS01.fa.gz --label PS01 \
    --chromosomes Chr1=chr1,Chr2=chr2,Chr3=chr3,Chr4=chr4,Chr5=chr5,Chr6=chr6,Chr7=chr7 \
    --hifi /data/PS01_hifi.fq.gz --ont /data/PS01_ont.fq.gz --ont-chemistry r10 \
    --reads-used-in-assembly yes \
    --outdir /data/qc --workdir /scratch/asmqc_PS01 --threads 32 --mem-gb 128
```

Options:

* `--assembly FILE` FASTA, gzip allowed. Required.
* `--label NAME` names all outputs; letters, digits, `.`, `_`, `-`. Required.
* `--outdir DIR` results go to `DIR/<label>/`. Required.
* `--agp FILE` AGP 2.0 of the assembly. Contigs are then the AGP components.
* `--chromosomes OLD=chrN,...` maps assembly names to `chr1`–`chr7`, in the FASTA and the AGP.
* `--illumina R1,R2[,R1b,R2b...]` paired reads, in order.
* `--hifi FILE[,FILE...]` HiFi reads, FASTQ or FASTA, gzip allowed.
* `--ont FILE[,FILE...]` ONT reads; needs `--ont-chemistry r9|r10`.
* `--reads-used-in-assembly yes|no` required when any reads are given.
* `--threads N` default (32). `--mem-gb N` default (128).
* `--workdir DIR` intermediates. Default (`<outdir>/<label>/work`).
* `--keep-intermediates` keep BAMs, k-mer databases and work files after success.
* `--modules 1,2,4,...` run only these modules. Default (all applicable).
* `--manifest-full-paths` record full input paths in the manifest, not only file names.
* `--dry-run` validate the inputs, print the plan and stop.

Re-running the same command resumes from the workdir.

Other subcommands:

* `asmqc version` prints the image version, git commit, tool versions and reference md5s.
* `asmqc aggregate DIR1 DIR2 ... --out combined/` merges result directories into one table and report.
* `asmqc test` runs the built-in smoke test on a synthetic assembly (about 1 min).

### Required declarations

**Chromosome names**: after `--chromosomes` mapping, the FASTA must contain
exactly `chr1`…`chr7`. Every other sequence is treated as unplaced. There is
no auto-detection: it would miscount silently on an assembly with a broken
chromosome.

**`--reads-used-in-assembly`**: when the reads were used to build the
assembly, the M8 QV measures self-consistency, not accuracy, and M11 finds
fewer errors. The declaration is recorded in every output. QVs are
comparable only between assemblies with the same declaration.

**`--ont-chemistry`**: selects the mapping preset (`r10` → `lr:hq`, `r9` →
`map-ont`). The R10 preset on R9 data under-aligns and creates false
breakpoints in M9.

### Exit codes

* `0` every applicable module finished. An ENA rule failure is a result, not an error.
* `1` input validation failed; nothing ran.
* `2` at least one module failed. The other results and `qc_summary.tsv` are still written; the failed module's columns are `NA`.

## Resources

Measured with asmqc 0.1.5 on a 4.3 Gb chromosome-level pea assembly
(7 chromosomes up to 752 Mb, 56 unplaced scaffolds) on an AMD EPYC 9654
host, workdir on local SSD:

| Run type | Threads, memory | Wall time | CPU time | Peak memory (largest step) | Work disk (peak) |
|---|---|---|---|---|---|
| Assembly only (M1, M2, M4–M7) | 32, 128 GB | 28 min | not measured | 23 GB (M7) | not measured |
| Illumina 58 Gb (~13×) + HiFi (4 cells), all modules | 64, 256 GB | 3 h 11 min | ~110 CPU-hours | 72 GB (meryl) | 333 GB |

Where the time goes in the run with reads:

| Stage | Wall time | Note |
|---|---|---|
| Read mapping (Illumina and HiFi at the same time) | 1 h 41 min | ~90 of the 110 CPU-hours; HiFi 41 GB, Illumina 33 GB memory |
| CRAQ (M9) | 1 h 09 min | at most 16 threads, 4 GB memory |
| meryl k-mer counting | 16 min | runs before mapping |
| All other modules | under 32 min each | run alongside CRAQ |

Run time grows mainly with read volume (mapping). Runs with ONT reads,
Illumina only or HiFi only, and deeper coverage have not been timed. The
assembly-only time was measured on 0.1.0rc3; BUSCO (20 GB) and M7 (23 GB)
can run at the same time.

The same run took 14 h 39 min with 0.1.1, with CRAQ at 12 h 21 min and
114 GiB. 0.1.2 to 0.1.5 replaced the CRAQ steps that held per-base depth in
memory; every release gave byte-identical results.

Each run records per-rule wall time and peak memory in `logs/benchmarks/`.
In `run_manifest.json`, `wall_seconds` gives each stage's elapsed time and
`rule_seconds` the sum of its rules' times (larger when a stage's rules run
at the same time). A free-space shortfall in `--workdir` is reported as a
warning before the run starts.

## Output

`<outdir>/<label>/`:

* `qc_summary.tsv` one row, fixed column order; `NA` where a module did not run.
* `flags.tsv` every flag (`ENA_BLOCKING`, `WARNING`, `INFO`) with sequence, position and message.
* `run_manifest.json` versions, inputs (file names and md5), parameters, reference data, host, wall times.
* `report.html` self-contained report; opens offline.
* `m01_integrity/` … `m11_homopolymer/` per-module tables and plots; a module that did not run has no directory.
* `logs/` tool logs; local paths are replaced with placeholders such as `<workdir>`.

After a successful run the intermediate files are deleted from the workdir,
unless `--keep-intermediates` is given. Other files in a user-supplied
`--workdir` are never touched.

### What to send back

Pack the result directory and send the archive to the consortium:

```
tar czf PS01_asmqc.tgz -C /data/qc --exclude PS01/work PS01
```

The archive is below 50 MB. It contains input file names but no local
paths. `--exclude` matters only with `--keep-intermediates` and the default
workdir, which lies inside the result directory.

## Reading the report

* The **header** shows the ENA result (PASS/FAIL), flag counts and the read declarations.
* **At a glance** and **per chromosome** summarise the key values.
* **Flags** lists ENA-blocking flags and warnings with the meaning of each code; information flags are folded.
* Each **module** section explains what is measured, then lists every value with a plain label, unit,
  definition and its `qc_summary.tsv` column, followed by tables and plots.
* **Glossary**, **output files** and **provenance** (command line, tool versions, reference md5s, parameters) close the report.

Hovering over a label shows its definition.

For chromosome-level assemblies the scaffold N50 is about one chromosome;
the contig N50 is the discriminating number. BUSCO completeness saturates
near 100 %; the duplicated percentage and the placement of duplicated
BUSCOs carry the signal. In M11, insertions relative to the assembly mean
that assembled homopolymer runs are too short.

## Development

```
mamba env create -f envs/dev.yaml && mamba activate asmqc-dev
pip install --no-deps --no-build-isolation -e .
python refs/fetch_refs.py .refs
ASMQC_REFS=.refs pytest
```

`docs/design/IMPLEMENTATION_PLAN.md` describes the code layout and
development decisions.

## Credits

* The rDNA subunit library is `data/rdna_library.fasta` from
  [kavonrtep/CARP](https://github.com/kavonrtep/CARP) 1.9.0 (GPL-3.0).
* The M6 synteny anchors are the BUSCO genes of *Pisum sativum* 'Caméor' v2
  (GCA_977071245.1; Kreplak et al., Scientific Data).
* The telomere karyoplot is adapted from `scripts/telomere_karyoplot.py` in
  [kavonrtep/ont_genome_assembly_pipeline](https://github.com/kavonrtep/ont_genome_assembly_pipeline)
  (GPL-3.0); origin and changes are in the header of `src/asmqc/karyoplot.py`.
* `workflow/craq_patch/src/` holds patched drivers and streaming script
  replacements for [CRAQ](https://github.com/JiaoLaboratory/CRAQ) 1.10 (MIT,
  `workflow/craq_patch/LICENSE.CRAQ`) with byte-identical output; origin in
  each file header.
* Tools: mm2-plus, minimap2, QUAST, tidk, BLAST+, BUSCO, miniprot, meryl,
  Merqury, CRAQ, samtools, bcftools, bedtools, seqkit, Snakemake.
* Reference data: BUSCO lineage `fabales_odb12.2` (OrthoDB 12.2),
  *Pisum sativum* chloroplast NC_014057.1, *L. oleraceus* mitochondrion
  PP555264.1.

## License

GPL-3.0; see [`LICENSE`](LICENSE).
