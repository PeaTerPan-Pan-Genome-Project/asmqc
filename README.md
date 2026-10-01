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
| M6 | gene-space completeness, `fabales_odb12.2` | always | BUSCO + miniprot |
| M7 | redundancy of unplaced scaffolds | unplaced ≥ 1 kb exist | mm2-plus |
| M8 | QV, k-mer completeness, spectra-cn | Illumina or HiFi given | meryl, Merqury |
| M9 | read-back structural validation | HiFi or ONT given | CRAQ |
| M11 | homopolymer and dinucleotide-repeat length errors | Illumina or HiFi given | bcftools |

**Read mapping**: each read type is mapped once with mm2-plus and shared by
the modules. M8 and M11 use Illumina reads if given, otherwise HiFi; ONT is
never used for them. M9 uses HiFi if given, otherwise ONT, plus Illumina when
given.

**Determinism**: the same inputs on the same image give identical values in
`qc_summary.tsv` (except `run_date`). Results are comparable only between
runs of the same `MAJOR.MINOR` image version.

## Installation

Pull the image and check its sha256 against the value published with the
release:

```
singularity pull asmqc_1.0.0.sif oras://ghcr.io/peaterpan-pan-genome-project/asmqc/sif:1.0.0
sha256sum asmqc_1.0.0.sif
```

The image needs no network at run time. The CPU must support AVX2
(required by mm2-plus); asmqc stops with an error otherwise.

To build the image yourself, run `./build.sh` in a clone of this repository
(Apptainer 1.4, network access during the build). The build downloads every
reference file, verifies its md5 and runs the built-in test.

## Usage

```
singularity run [-B <bind paths>] asmqc_<version>.sif run \
    --assembly ASM.fa[.gz] --label SAMPLE01 --outdir results/ [options]
```

Bind every directory that holds inputs or outputs (`-B /data,/scratch`).

**Example 1, assembly only** (M1, M2, M4–M7; about 6–12 h for 4.3 Gb):

```
singularity run -B /data asmqc_1.0.0.sif run \
    --assembly /data/PS01.fa.gz --agp /data/PS01.agp --label PS01 \
    --outdir /data/qc --threads 32 --mem-gb 64
```

**Example 2, with Illumina reads** (adds M8 and M11):

```
singularity run -B /data,/scratch asmqc_1.0.0.sif run \
    --assembly /data/PS01.fa.gz --agp /data/PS01.agp --label PS01 \
    --illumina /data/PS01_R1.fq.gz,/data/PS01_R2.fq.gz \
    --reads-used-in-assembly no \
    --outdir /data/qc --workdir /scratch/asmqc_PS01 --threads 32 --mem-gb 128
```

**Example 3, with HiFi and ONT reads** (adds M8, M9 and M11, all from HiFi).
ONT reads are mapped only when no HiFi reads are given; with both, the ONT
files are recorded in the manifest but not used:

```
singularity run -B /data,/scratch asmqc_1.0.0.sif run \
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

Estimates for a 4.3 Gb assembly on 32 threads; measured values will replace
them before v1.0.

| Run type | Time | RAM | Workdir |
|---|---|---|---|
| Assembly only (M1–M7) | ~6–12 h | 64 GB | ~60 GB |
| With Illumina and long reads | ~1–2 days | 128 GB | ~500 GB |

Put `--workdir` on fast scratch with enough space. A free-space shortfall is
reported as a warning before the run starts.

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
* **Flags** lists ENA-blocking flags and warnings; information flags are folded.
* Each **module** section gives its numbers, a short paragraph on how to read them, and plots.
* **Provenance** lists the command line, tool versions, reference md5s and parameters.

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
* The telomere karyoplot is adapted from `scripts/telomere_karyoplot.py` in
  [kavonrtep/ont_genome_assembly_pipeline](https://github.com/kavonrtep/ont_genome_assembly_pipeline)
  (GPL-3.0); origin and changes are in the header of `src/asmqc/karyoplot.py`.
* Tools: mm2-plus, minimap2, QUAST, tidk, BLAST+, BUSCO, miniprot, meryl,
  Merqury, CRAQ, samtools, bcftools, bedtools, seqkit, Snakemake.
* Reference data: BUSCO lineage `fabales_odb12.2` (OrthoDB 12.2),
  *Pisum sativum* chloroplast NC_014057.1, *L. oleraceus* mitochondrion
  PP555264.1.

## License

GPL-3.0; see [`LICENSE`](LICENSE).
