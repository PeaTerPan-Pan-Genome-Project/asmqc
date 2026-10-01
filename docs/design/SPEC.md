# asmqc: specification

**asmqc** is a single Singularity image that runs a fixed, comparable set of
quality-control checks on chromosome-level genome assemblies of **pea**
(*Pisum*; *Lathyrus oleraceus* and *L. fulvus*). It is built for a pangenome
consortium in which every group runs QC on its own assemblies.

Because every group uses identical tools, versions, parameters and reference
data, and writes one fixed output schema, the results can be merged into one
comparable table.

This document is the design specification for version 1.0. Decisions recorded
here are project decisions. Points marked **[default]** are implementation
defaults that may change; §13 lists the open points that go to the
maintainers.

---

## 1. Scope

**In scope:** QC of a chromosome-level genome assembly: 7 pseudomolecules named
`chr1`–`chr7`, plus unplaced scaffolds.

**Out of scope:**
- Gene annotation QC.
- Contamination screening (NCBI FCS-GX / FCS-adaptor; run separately).
- Genome size against flow cytometry.
- Centromere and satellite validation.
- LAI.
- Any modification of the assembly. asmqc **measures**; it never writes a
  corrected assembly.

**Species:** pea only. Reference data and motifs are fixed for pea; there are
no per-species presets.

---

## 2. Design principles

| | Principle |
|---|---|
| A | **All reference data inside the image.** One image digest means identical tools and identical data. |
| B | **One command:** `singularity run asmqc.sif …`. A Snakemake workflow inside switches modules on according to the inputs supplied. It runs on a single machine; there is no scheduler integration (users wrap it in their own scheduler). |
| C | **Machine-readable output plus a self-contained HTML report**: a fixed-schema TSV and JSON. |
| D | **Pea only.** |
| E | **Comparability boundary = minor version.** Results are comparable only between runs of the same `MAJOR.MINOR` image. |

---

## 3. Modules

| # | Module | Runs when | Programs |
|---|---|---|---|
| 1 | Integrity, format and ENA rules | always | Python, `seqkit`, `samtools faidx` |
| 2 | Contiguity and anchoring | always | QUAST 5.3.0, Python (AGP) |
| 4 | Telomeres | always | tidk 0.2.65, Python |
| 5 | Organelle and rDNA inventory | always | mm2-plus 1.3, BLAST+ 2.17.0 |
| 6 | Gene-space completeness | always | BUSCO 6.1.0 + miniprot 0.18 |
| 7 | Redundancy of unplaced scaffolds | always (if any unplaced ≥ 1 kb) | mm2-plus 1.3 |
| 8 | QV, k-mer completeness, spectra-cn | Illumina or HiFi given | meryl 1.4.2, Merqury 1.4.1 |
| 9 | Read-back structural validation | HiFi or ONT given | CRAQ 1.10 on mm2-plus BAMs |
| 11 | Homopolymer and short-STR errors | Illumina or HiFi given | bcftools 1.24 on the mm2-plus BAM |

Module numbers 3 and 10 are intentionally unused; they were taken out of scope.
The numbering is kept because it appears in column prefixes.

### 3.1 Read mapper: mm2-plus

Every minimap2-type step uses **mm2-plus**
([at-cg/mm2-plus](https://github.com/at-cg/mm2-plus); bioconda `mm2plus` 1.3).
- It is a drop-in replacement built on minimap2 2.31, using AVX2/AVX-512
  alignment and parallel sorting, with 1.6–7.2× speedup on whole-genome
  alignment.
- Read mapping is identical to minimap2 2.31 with `--max-chain-skip=1000000`.
  Genome-to-genome alignment shows negligible differences, which do not affect
  comparability because everyone runs the same binary.
- Confirm the binary name in the package. Record the host CPU flags in the
  manifest.

**GPU-accelerated mm2-gb was considered and rejected:**
- users' hardware cannot be assumed;
- the CPU side is single-threaded;
- it is built on minimap2 2.24;
- it has no package.

A GPU path would split users into two code paths.

---

## 4. User interface

### 4.1 Command

```
singularity run [-B <bind paths>] asmqc_<version>.sif run \
    --assembly  ASM.fa[.gz]            # required
    --label     SAMPLE01               # required; [A-Za-z0-9._-]+, used in all output names
    --outdir    results/               # required; results go to <outdir>/<label>/
    [--agp      ASM.agp]
    [--chromosomes old1=chr1,old2=chr2,...]   # map names to chr1..chr7 (FASTA and AGP)
    [--illumina R1.fq.gz,R2.fq.gz[,R1b.fq.gz,R2b.fq.gz ...]]   # pairs, in order
    [--hifi     reads1.fq.gz[,reads2...]]     # FASTQ/FASTA, gz ok
    [--ont      reads1.fq.gz[,...] --ont-chemistry r9|r10]
    [--reads-used-in-assembly yes|no]         # REQUIRED if any reads are given
    [--threads 32] [--mem-gb 128]
    [--workdir /scratch/asmqc_SAMPLE01]       # intermediates; default <outdir>/<label>/work
    [--keep-intermediates]                    # default: delete BAMs, meryl DBs, work files on success
    [--modules 1,2,4,5,6,7,8,9,11]            # default: all applicable
    [--dry-run]                               # validate inputs, print the plan, run nothing
```

Other subcommands:

| Subcommand | Purpose |
|---|---|
| `asmqc version` | image version, git commit, every tool version, reference-data md5s |
| `asmqc aggregate <dir1> <dir2> … --out combined/` | merge several `<label>/` result dirs into `combined/qc_summary.tsv` and `combined/report.html` (one row per assembly) |
| `asmqc test` | runs the built-in smoke test (§11.1) |

### 4.2 Input validation

Validation runs before anything else, and also under `--dry-run`. Any failure
exits with **code 1** and a clear message.

- The assembly exists and is readable. gzip is fine; it is decompressed once
  into the workdir.
- **Chromosome naming:** after `--chromosomes` mapping, the FASTA must contain
  **exactly `chr1`…`chr7`**; there is no auto-detection. Everything else is
  "unplaced". If any of `chr1`–`chr7` is missing, exit 1, listing the names
  found. Auto-detection would silently miscount on an assembly with a broken
  chromosome.
- If any reads are given, `--reads-used-in-assembly` is required.
- If `--ont` is given, `--ont-chemistry` is required.
- The output dir is writable. Free workdir space below the §9 estimate is a
  warning, not a failure.
- **If AVX2 is absent, exit 1** (mm2-plus).

### 4.3 Read-type selection (fixed rules)

| Purpose | Rule |
|---|---|
| Module 8 (QV) and module 11 | **Illumina if given, else HiFi.** ONT is never used for these. |
| Module 9 (CRAQ) long reads | **HiFi if given, else ONT.** Illumina is additionally passed as `-ngs` if given. |

The read type actually used is recorded in each module's output.

### 4.4 Exit codes

| Code | Meaning |
|---|---|
| 0 | every applicable module finished; results written |
| 1 | input validation failed; nothing run |
| 2 | at least one module failed. The other modules' results are still written and `qc_summary.tsv` is produced, with that module's columns `NA` and `mNN_status = failed` |

**An ENA rule FAIL is a result, not an error:** exit code 0.

---

## 5. Architecture

### 5.1 Image layout

- Base: `docker://mambaorg/micromamba:2.9.0` **[default]**.
- Separate micromamba environments, because some tools have conflicting
  dependencies. Every environment is built from a **committed explicit lock
  file**.

| Environment | Pinned contents |
|---|---|
| `/opt/envs/core` | python 3.12, snakemake 9.27.0, pandas, matplotlib, jinja2, seqkit 2.14.0, samtools 1.24, bcftools 1.24, bedtools 2.31.1, mm2plus 1.3, blast 2.17.0, tidk 0.2.65 |
| `/opt/envs/quast` | quast 5.3.0 |
| `/opt/envs/busco` | busco 6.1.0, miniprot 0.18 (plus BUSCO's dependencies) |
| `/opt/envs/merqury` | merqury 1.4.1, meryl 1.4.2 (plus R as packaged) |
| `/opt/envs/craq` | craq 1.10 |

- The lock files are the pins for any library not listed here.
- Reference data in `/opt/asmqc/refs/` (§6).
- Workflow and scripts in `/opt/asmqc/`.
- `%runscript` → `/opt/asmqc/bin/asmqc "$@"`.
- **No network at run time.** BUSCO runs `--offline`, nothing is downloaded,
  and the image must work on an offline node.

### 5.2 Workflow

- Snakemake, called by the `asmqc` wrapper with `--cores`,
  `--resources mem_mb`, `--rerun-incomplete` and `--keep-going`.
- Re-running the same command resumes from the workdir.
- Each rule puts its environment's `bin/` on PATH explicitly. There is no
  conda activation at run time.
- **Shared stages**, computed once:
  - `prep`: decompress and index the assembly; apply `--chromosomes` to a
    working copy. The sequence content is unchanged, and the renaming map goes
    into the manifest.
  - `map_<readtype>`: one mm2-plus mapping per read type, then
    `samtools sort` + `index`, with read groups and
    `--max-chain-skip=1000000`:

    | Reads | Preset |
    |---|---|
    | HiFi | `-ax map-hifi` |
    | ONT, `r10` | `-ax lr:hq` |
    | ONT, `r9` | `-ax map-ont` |
    | Illumina | `-ax sr` (paired) |

  - `meryl_reads`: the read k-mer database (k = 21).
- All temporary files go under `--workdir`, never `/tmp`: set `TMPDIR` for
  every rule.
- On success, without `--keep-intermediates`, delete the BAMs, meryl DBs and
  scratch. Keep only `<outdir>/<label>/`.

### 5.3 Determinism

- The same inputs on the same image give identical `qc_summary.tsv` values.
  Wall times belong in the manifest only, never in the summary.
- Seeds are fixed wherever a tool takes one.
- Tool outputs whose order varies with threading are sorted before parsing.

---

## 6. Reference data (inside the image)

The build downloads each item and **verifies its md5; the build fails on any
mismatch**. All md5s are written to the run manifest.

| Item | Source | md5 | Notes |
|---|---|---|---|
| BUSCO lineage `fabales_odb12.2` | `https://busco-data.ezlab.org/v6/data/lineages/fabales_odb12.2.2026-05-13.tar.gz` (387,854,545 B) | `5505de6096f169497b957fd56475d434` | `dataset.cfg`: creation 2026-05-13, **7,702 BUSCOs**, 11 species, OrthoDB 12.2. Unpack under `/opt/asmqc/refs/busco_downloads/lineages/fabales_odb12.2`. |
| Plastid genome | NCBI **NC_014057.1** (*Pisum sativum* chloroplast, complete genome, 122,169 bp) | `959a1c9be9e5c9e21a47a7abd31fd295` (efetch FASTA, 2026-10-01) | |
| Mitochondrial genome | NCBI **PP555264.1** (*L. oleraceus* subsp. *oleraceus* voucher JI_281 mitochondrion, complete genome, 363,796 bp) | `3b6f306756e911a3ed2ff7c5005d3a03` | |
| rDNA subunit library | [`kavonrtep/CARP`](https://github.com/kavonrtep/CARP) tag `1.9.0`, `data/rdna_library.fasta` | `f09e64b09f99f83a44653a6a9153373b` | GPL-3.0. 117 sequences: 23 × 18S, 18 × 5.8S, 17 × 25S, 59 × 5S. Headers are `>id#rDNA/45S_rDNA/18S`. |
| Telomere motif | constant `TTTAGGG` | — | |
| Karyoplot code | [`kavonrtep/ont_genome_assembly_pipeline`](https://github.com/kavonrtep/ont_genome_assembly_pipeline) `scripts/telomere_karyoplot.py` @ `2144dc8b467b7c606c52b6e8667b8ce24640b077` | — | GPL-3.0; adapted (§8.4), with its origin in the file header |

If NCBI's FASTA formatting ever changes the file md5, verify the md5 of the
**sequence** and record both. Never accept a different sequence silently.

---

## 7. Output contract

### 7.1 Results directory

```
<outdir>/<label>/
  qc_summary.tsv            # ONE row, fixed schema (§7.2)
  run_manifest.json         # §7.3
  report.html               # self-contained (§7.4)
  flags.tsv                 # every flag from every module (§7.5)
  m01_integrity/   integrity.tsv  sequences.tsv  gaps.tsv  checksums.md5
  m02_contiguity/  contiguity.tsv  per_chromosome.tsv  quast_report.tsv
  m04_telomeres/   telomeres.tsv  telomeres_unplaced.tsv  interstitial.tsv  tidk_windows.tsv  karyoplot.png
  m05_organelle_rdna/ organelle_scaffolds.tsv  organelle_on_chromosomes.tsv  rdna_arrays.tsv
  m06_busco/       short_summary.json  full_table.tsv  busco_derived.tsv
  m07_redundancy/  redundancy.tsv
  m08_merqury/     merqury.qv  per_chromosome_qv.tsv  completeness.stats  asm_only_kmers.bed.gz  spectra-cn.png  spectra-asm.png
  m09_craq/        out_final.Report  CRE.bed  CSE.bed  craq_derived.tsv
  m11_homopolymer/ errors.tsv  errors.bed.gz  hom_calls.vcf.gz
  logs/
```

- A module that did not run has no directory.
- QUAST's own HTML/PDF/Icarus output is not included.
- BAMs, meryl DBs and BUSCO raw output stay in the workdir.
- Target size: **< 50 MB**. This directory is what users send to the
  consortium.

### 7.2 `qc_summary.tsv`: schema

- One header line plus one row; tab-separated; UTF-8.
- **Column order is fixed** as listed.
- `NA` for a module that did not run or a metric that does not apply.
- Percentages are 0–100 with 2 decimals. Counts and bp are integers. Rates have
  4 significant digits.
- `asmqc aggregate` refuses inputs whose headers differ.

**Identity and run:**
- `label`, `asmqc_version`, `assembly_md5` (of the decompressed bytes),
  `run_date`
- `ena_rules` (`PASS`/`FAIL`)
- `n_flags_ena_blocking`, `n_flags_warning`, `n_flags_info`
- `m01_status`, `m02_status`, `m04_status`, `m05_status`, `m06_status`,
  `m07_status`, `m08_status`, `m09_status`, `m11_status`: each one of `ok` ·
  `skipped_no_input` · `skipped_by_user` · `failed`

**Module 1:**
- `m01_n_seq`, `m01_total_bp`, `m01_n_lt20bp`, `m01_n_lt200bp`,
  `m01_n_terminal_n`, `m01_n_duplicate_names`, `m01_n_invalid_chars`,
  `m01_n_iupac`, `m01_n_gaps_ge10`, `m01_gap_bp_ge10`, `m01_n_gaps_ge100`,
  `m01_n_seq_gt50pct_n`, `m01_softmask_pct`
- `m01_agp_consistent` (`yes`/`no`/`NA`)
- `m01_n_round_lengths`, `m01_round_lengths_expected`
- `m01_n_round_agp_cuts`, `m01_n_agp_cuts`

**Module 2:**
- `m02_scaffold_n50`, `m02_scaffold_l50`, `m02_scaffold_n90`, `m02_longest_bp`
- `m02_contig_method` (`agp`/`nsplit10`), `m02_contig_n50`, `m02_contig_l50`,
  `m02_contig_n90`, `m02_n_contigs`, `m02_contig_n50_nsplit10`
- `m02_chrom_bp`, `m02_anchored_pct`, `m02_n_unplaced`, `m02_unplaced_bp`,
  `m02_n_gaps`, `m02_gap_bp`

**Module 4:**
- `m04_capped_arms` (0–14), `m04_t2t_chromosomes` (0–7)
- `m04_wrong_orientation_arms`, `m04_interstitial_arrays`,
  `m04_unplaced_with_telomere`

**Module 5:**
- `m05_plastid_scaffolds_n`, `m05_plastid_scaffolds_bp`,
  `m05_mito_scaffolds_n`, `m05_mito_scaffolds_bp`
- `m05_chrom_plastid_like_bp`, `m05_chrom_mito_like_bp`
- `m05_rdna45s_loci` (semicolon list), `m05_rdna45s_copies`,
  `m05_rdna5s_loci`, `m05_rdna5s_copies`
- `m05_rdna_only_scaffolds_n`, `m05_rdna_only_scaffolds_bp`

**Module 6:**
- `m06_complete_pct`, `m06_single_pct`, `m06_duplicated_pct`,
  `m06_fragmented_pct`, `m06_missing_pct`, `m06_n_markers`
- `m06_internal_stop_pct`
- `m06_complete_on_unplaced`
- `m06_dup_both_on_chrom`, `m06_dup_any_on_unplaced`
- `m06_lineage` (`fabales_odb12.2 2026-05-13`)

**Module 7:**
- `m07_n_duplicate`, `m07_duplicate_bp`
- `m07_n_partial_overlap`, `m07_partial_overlap_bp`
- `m07_n_repeat_like`, `m07_repeat_like_bp`
- `m07_n_unique`, `m07_unique_bp`
- `m07_n_short`, `m07_short_bp`
- `m07_total_minus_duplicate_bp`

**Module 8:**
- `m08_read_type` (`illumina`/`hifi`), `m08_reads_independent` (`yes`/`no`)
- `m08_k` (21), `m08_kmer_coverage`, `m08_low_coverage` (`yes`/`no`)
- `m08_qv`, `m08_error_rate`, `m08_completeness_pct`

**Module 9:**
- `m09_long_read_type` (`hifi`/`ont_r9`/`ont_r10`), `m09_short_reads_used`
  (`yes`/`no`), `m09_coverage`, `m09_low_coverage`
- `m09_aqi`, `m09_r_aqi`, `m09_s_aqi`, `m09_n_cre`, `m09_n_cse`
- `m09_cse_near_agp_junction`, `m09_cse_inside_contig`

**Module 11:**
- `m11_read_type`, `m11_reads_independent`, `m11_callable_bp`
- `m11_hp_errors`, `m11_hp_errors_per_mb`, `m11_hp_errors_per_10k_runs`
- `m11_dinuc_errors`, `m11_dinuc_errors_per_mb`
- `m11_other_indels_per_mb`, `m11_snv_per_mb`
- `m11_hp_pct_of_errors`, `m11_hp_ins_del_ratio`, `m11_hp_at_pct`
- `m11_het_calls`

### 7.3 `run_manifest.json`

```json
{
  "asmqc_version": "1.0.0", "git_commit": "…", "lockfile_sha256": {"core": "…"},
  "label": "…", "command_line": "…", "start": "…", "end": "…",
  "wall_seconds": {"m01": 0},
  "host": {"cpu_model": "…", "cpu_flags": ["avx2", "avx512f"], "threads": 32, "mem_gb": 128},
  "inputs": {
    "assembly": {"name": "…", "md5": "…", "bytes": 0},
    "agp": {"name": "…", "md5": "…"},
    "reads": [{"type": "illumina", "files": [{"name": "…", "bytes": 0}], "reads_used_in_assembly": "no"},
              {"type": "ont", "chemistry": "r10", "files": [], "reads_used_in_assembly": "yes"}]
  },
  "chromosome_map": {"old1": "chr1"},
  "tools": {"mm2plus": "1.3", "busco": "6.1.0"},
  "reference_data": {"fabales_odb12.2": {"date": "2026-05-13", "markers": 7702, "md5": "5505de60…"}},
  "parameters": {"m02": {"nsplit_min_gap": 10}, "m04": {"window": 10000, "min_repeats": 25, "terminal_bp": 50000}},
  "modules": {"m08": {"status": "ok", "read_type": "illumina", "reason": null}}
}
```

- Tool versions are read from the tools at run time. Exception: CRAQ, whose
  1.10 package prints `1.0.9-alpha`; its version is read from
  `/opt/envs/craq/conda-meta`.
- Every parameter named in §8 appears under `parameters`.
- Input file **names**, not full paths, are recorded **[default]**, so a shared
  result does not expose the user's directory layout. `--manifest-full-paths`
  opts in to full paths.

### 7.4 `report.html`

**One self-contained file**: inline CSS, embedded PNGs, no external requests.
Sections:

1. **Header card:** label, assembly md5, version, date, `ena_rules` badge,
   flag counts, and the read types with their `reads_independent`
   declarations.
2. **Flags:** every `ENA_BLOCKING` and `WARNING` flag (`INFO` collapsed).
3. **Per module:** key numbers, a short "how to read this" paragraph, and
   plots:
   - M2: cumulative contig length
   - M4: telomere karyoplot
   - M7: bp per class
   - M8: spectra-cn and spectra-asm
   - M11: errors by homopolymer length, stacked A/T vs G/C
4. **Provenance:** versions, reference data, parameters and the command line.

The combined report from `asmqc aggregate`:
- one table with assemblies as rows and the columns grouped by module;
- small multiples of the M4 karyoplots and the M11 plots;
- M6 internal stop-codon % plotted against the M11 HP-error rate.

### 7.5 Flags: `flags.tsv`

Columns: `module`, `code`, `severity` (`ENA_BLOCKING` | `WARNING` | `INFO`),
`seq_id`, `start`, `end`, `value`, `message`.

| Code | Severity | Module |
|---|---|---|
| `seq_lt_20bp` | ENA_BLOCKING | M1 |
| `terminal_n` | ENA_BLOCKING | M1 |
| `duplicate_name` | ENA_BLOCKING | M1 |
| `invalid_char` (non-IUPAC) | ENA_BLOCKING | M1 |
| `empty_sequence` | ENA_BLOCKING | M1 |
| `seq_lt_200bp` | WARNING | M1 |
| `n_fraction_gt_50pct` | WARNING | M1 |
| `agp_mismatch` | WARNING | M1 |
| `crlf_line_endings` | WARNING | M1 |
| `iupac_present` | INFO | M1 |
| `inconsistent_line_width` | INFO | M1 |
| `round_length` | INFO | M1 |
| `round_agp_cut` | INFO | M1 |
| `organelle_scaffold` | WARNING (ENA accepts declared organelle sequence) | M5 |
| `wrong_orientation_telomere` | INFO | M4 |
| `interstitial_telomere` | INFO | M4 |
| `low_coverage` | WARNING | M8, M9 |
| `reads_not_independent` | INFO | M8, M11 |

`ena_rules` is FAIL if and only if there is at least one `ENA_BLOCKING` flag.
Round-number diagnostics never affect it.

---

## 8. Module specifications

Throughout this section:
- "chromosome" means `chr1`–`chr7` after mapping;
- "unplaced" means every other sequence;
- identity is matches ÷ alignment block length (PAF columns 10 ÷ 11, with
  `-c`).

### 8.1 Module 1: integrity, format and ENA rules

**Per sequence (`sequences.tsv`):**
- name; length
- counts of A, C, G, T, N (upper + lower case), IUPAC ambiguity codes, and
  other characters
- leading and trailing N-run lengths
- N fraction; soft-masked fraction

**`gaps.tsv`:** every N-run of ≥ 10 bp, with seq, start, end and length.

**Whole file:**
- duplicate names; name characters (printable, no whitespace)
- CRLF line endings; line-width consistency
- md5 of the FASTA (and of the AGP) in `checksums.md5`

**AGP check (when `--agp` is given):**
1. Parse the AGP 2.0 `W` and `N`/`U` lines.
2. Map object names through `--chromosomes`.
3. Every FASTA sequence must appear as an object of the same length, object
   coordinates must be contiguous from 1, and the components must cover every
   base.
4. On any mismatch, raise `agp_mismatch` with the details, set
   `m01_agp_consistent = no`, and have module 2 use `nsplit10`. One common
   real-world case is an AGP whose object names differ from the FASTA headers.

**Round-number diagnostic** (INFO only; never affects `ena_rules`).
Manual curation of Hi-C contact maps can place breaks on a fixed coordinate
grid (e.g. 1 kb) rather than at the true junction. That leaves round-length
pieces and very short stubs behind.
- Count the unplaced sequences whose length is an exact multiple of 1,000,
  against the expectation n_unplaced / 1000.
- With an AGP (even one that failed the FASTA check), collect for every `W`
  line:
  1. `component_beg − 1`, if `component_beg > 1`;
  2. if the component name contains `_subseq_A:B` (a piece cut from a parent
     contig, coordinates on the parent): `A − 1` if `A > 1`, and `B`;
  3. otherwise `component_end`.

  Report how many of the collected coordinates are multiples of 1,000
  (`m01_n_round_agp_cuts`) out of all collected (`m01_n_agp_cuts`). The chance
  rate is ~0.1 %.

**Rules:**

| Rule | Severity |
|---|---|
| < 20 bp | ENA_BLOCKING |
| < 200 bp | WARNING |
| leading or trailing N | ENA_BLOCKING |
| > 50 % N | WARNING |

### 8.2 Module 2: contiguity and anchoring

```
quast.py --large --min-contig 0 --split-scaffolds --no-icarus --no-plots \
         --threads N -o work/quast ASM.fa
```
- **No reference**: comparing against another accession would mix biology
  into QC.
- No `--conserved-genes-finding` and no `--gene-finding`: module 6 covers
  completeness.
- `--split-scaffolds` splits at ≥ 10 N; record `nsplit_min_gap: 10`.
- Keep `report.tsv` as `quast_report.tsv`. The Python-computed scaffold N50
  and contig N50 from the N-split must agree with QUAST's; a disagreement is a
  bug.

**Contigs:**
- If the AGP is consistent, contigs are its `W` components (`agp`).
- Otherwise contigs are the pieces between N-runs of ≥ 10 bp (`nsplit10`).
- `m02_contig_n50_nsplit10` is always reported. A difference from the AGP
  value (N-runs inside components) is informative.

**Anchoring:**
- `m02_chrom_bp` = bp in `chr1`–`chr7`
- `m02_anchored_pct` = `m02_chrom_bp` ÷ total × 100

**`per_chromosome.tsv`:** chromosome, length, contigs, gaps, gap bp, contig
N50.

For chromosome-level assemblies scaffold N50 is about one chromosome; contig
N50 is the discriminating number. The report says so.

### 8.4 Module 4: telomeres

1. `tidk search --string TTTAGGG --window 10000 --output <label> --dir work/tidk ASM.fa`
   (tidk 0.2.65). Copy the output to `tidk_windows.tsv`.
2. **Telomeric window:** forward + reverse ≥ **25** repeats. Merge adjacent
   telomeric windows into bands, recording for each band its start, end,
   Σforward and Σreverse.
3. Classify each chromosome arm:
   - **start arm:** the band nearest position 0.
     - Within **50,000 bp** of the start and dominated by `CCCTAAA` →
       `capped`.
     - Within 50 kb but dominated by the other orientation →
       `wrong_orientation`.
     - No band within 50 kb → `absent`.
   - **end arm:** the mirror image; dominated by `TTTAGGG`.
   - Confirm tidk's forward/reverse convention on real data (§13).
4. **T2T chromosome:** both arms `capped`.
5. **Approximate array length:** (Σforward + Σreverse) × 7 bp of the terminal
   band. **Distance from the end:** band start, or sequence length − band end.
6. **Interstitial arrays:** bands more than 50 kb from both ends →
   `interstitial.tsv` and INFO. Pea may carry genuine interstitial telomeric
   repeats.
7. **Unplaced sequences with a band within 50 kb of an end** →
   `telomeres_unplaced.tsv`, listed and not plotted. They are chromosome ends
   that exist but were not anchored.

**Karyoplot** (adapted from the source in §6):
- draw `chr1`–`chr7` only;
- draw every band, coloured by status;
- use the 50 kb rule; no reorientation.

`telomeres.tsv` columns: chromosome, arm, status, distance_from_end_bp,
approx_array_bp, fwd_repeats, rev_repeats.

### 8.5 Module 5: organelles and rDNA

**Organelles:**
1. `mm2plus -x asm20 -c -t N refs/organelles.fa ASM.fa > org.paf`, with the
   assembly as query.
2. Keep records with identity ≥ 0.95.
3. Per query, compute the union of covered bp, separately for plastid and
   mito.
4. **Unplaced sequence ≥ 80 % covered by one organelle:**
   `organelle_scaffold` (WARNING), typed by the larger coverage.
5. **Chromosomes:** plastid-like bp, mito-like bp and the largest block per
   chromosome (`organelle_on_chromosomes.tsv`). **Report only**: nuclear
   insertions of organelle DNA are real biology, but a very large block can
   indicate a misjoin.
6. No *P. fulvum* organelle references exist; the 95 % threshold accommodates
   *L. fulvus*.

**rDNA:**
1. Search the assembly with BLAST:
   ```
   makeblastdb -dbtype nucl -in ASM.fa
   blastn -task blastn -query refs/rdna_library.fasta -db ASM -evalue 1e-10 \
          -outfmt "6 std qlen slen" -max_target_seqs 1000000 -max_hsps 1000000 -num_threads N
   ```
   The 5S unit (~120 bp) is too short for minimap2.
2. Take the subclass (18S/5.8S/25S/5S) from the query header after `#`.
3. **Copies:** merge overlapping hits of the same subclass into one copy.
   45S copies = 18S copies; 5S copies = 5S copies.
4. **Arrays:** single-linkage clustering of copies of the same family within
   **20 kb**.
5. `rdna_arrays.tsv`: seq_id, start, end, family, copies, and context:
   `chromosome`, `unplaced`, or `rdna_only_scaffold` (the array spans ≥ 80 % of
   an unplaced sequence).
6. **Report only:**
   - The 45S NORs are expected on chr4 and chr7 in pea; never PASS/FAIL.
   - Assembled copy number measures how much of the array was captured, not
     the plant.

### 8.6 Module 6: BUSCO

```
busco --in ASM.fa --mode genome --lineage_dataset fabales_odb12.2 --offline \
      --opt-out-run-stats \
      --download_path /opt/asmqc/refs/busco_downloads --cpu N --out busco --out_path work/ -f
```
- BUSCO 6.1.0 with the eukaryote default predictor, **miniprot**; assert it.
- `--opt-out-run-stats`: BUSCO 6 otherwise sends anonymous usage data, which
  breaks "no network at run time" (§5.1).
- **Assert the lineage `fabales_odb12.2`, dated 2026-05-13, with 7,702
  markers**; the module fails otherwise.
- From `short_summary.json`: C / S / D / F / M %, number of markers, and the
  internal stop-codon %.

**Derived values, from `full_table.tsv`** (no extra run):
- `m06_complete_on_unplaced`: Complete BUSCOs with a hit on an unplaced
  sequence.
- Duplicated BUSCOs: all copies on chromosomes (`m06_dup_both_on_chrom`) vs at
  least one copy unplaced (`m06_dup_any_on_unplaced`, a likely false
  duplication).
- Per-BUSCO detail in `busco_derived.tsv`; `full_table.tsv` is kept.

**Report text:** in high-quality assemblies C % saturates near 100, so D % and
the derived numbers carry the signal.

### 8.7 Module 7: redundancy of unplaced scaffolds

**This is a measurement, never a purge.**
1. `mm2plus -x asm5 -c --secondary=yes -N 5 -t N chroms.fa unplaced_ge1kb.fa > red.paf`.
2. Classify each query:

   | Class | Rule |
   |---|---|
   | `duplicate` | a **single primary PAF record** (one collinear chain to one locus) spanning ≥ 90 % of the query, identity ≥ 0.99, MAPQ ≥ 20 |
   | `partial_overlap` | the same criteria, spanning 50–90 % |
   | `repeat_like` | no qualifying single record ≥ 50 %, but the union of query coverage over **all** records (identity ≥ 0.99, any MAPQ) is ≥ 50 % |
   | `unique` | otherwise |
   | `short` | < 1 kb: not aligned, counted only |

**Why one chain and not summed coverage.** The pea genome is about 85 % repeat,
with young, near-identical 5–15 kb LTR retrotransposons. Summing coverage over
all alignments would call purely repetitive scaffolds duplicates.

3. `redundancy.tsv`: scaffold, length, class, chromosome, target start,
   target end, identity, MAPQ, coverage_frac.
4. `m07_total_minus_duplicate_bp` = total − `m07_duplicate_bp`.
5. **Cross-check (report only):** how many `m06_dup_any_on_unplaced` BUSCOs
   sit on `duplicate` scaffolds.

### 8.8 Module 8: QV, k-mer completeness, spectra-cn

**Reads:** Illumina, else HiFi. **Never ONT**: its residual errors are
concentrated in homopolymers, the same class that persists in ONT consensus
sequence, so the QV would measure the reads as much as the assembly.

**Steps:**
1. Count k-mers per read file and merge:
   ```
   meryl count k=21 memory=<mem> threads=N output <f>.meryl <file>
   meryl union-sum output reads.meryl *.meryl
   ```
   **k = 21 is fixed** for all users.
2. In the module work dir: `merqury.sh reads.meryl ASM.fa <label>`.

**Parse:**
- QV and error rate
- per-chromosome QV
- completeness %
- the spectra-cn and spectra-asm plots
- the assembly-only k-mer BED (gzipped)

**Coverage:** from the first major peak of `meryl histogram` after the error
trough. Below 20× → `low_coverage` (WARNING), and the module still runs.

**Declaration:** `m08_reads_independent` is the inverse of
`--reads-used-in-assembly`.
- If the reads were used to build the assembly, the QV is labelled
  **"self-consistency QV"**, with INFO `reads_not_independent`.
- QVs are comparable only between assemblies with the same declaration.

### 8.9 Module 9: CRAQ

Inspector was considered and not included:
- it accepts only raw reads and maps them itself, which would be a second,
  inconsistent mapping;
- its QV duplicates module 8 and its structural errors duplicate CRAQ.

**Inputs:**
- `-sms`: the shared long-read BAM (HiFi, else ONT).
- `-ngs`: the shared Illumina BAM, if Illumina was given.
- ONT preset from `--ont-chemistry`: `r10` → `lr:hq`, `r9` → `map-ont`.
  `lr:hq` on R9 data would under-align and invent breakpoints.

**Command:** `craq -g ASM.fa -sms long.sorted.bam [-ngs short.sorted.bam] -t N`
with CRAQ defaults.

**Parse `out_final.Report`:** AQI, R-AQI, S-AQI, CRE count, CSE count. Keep
the CRE and CSE BEDs.

**Coverage:** median long-read depth over chromosomes. Below 20× →
`low_coverage`.

**AGP-junction cross-check** (diagnostic):
- CSEs within 10 kb of an AGP gap or W–W junction, against CSEs inside a
  component.
- A breakpoint at a Hi-C join suggests a scaffolding error; one inside a contig
  suggests an assembler error.

### 8.11 Module 11: homopolymer and short-STR errors

**Background.**
- Assemblies built from Oxford Nanopore reads can systematically misreport the
  length of homopolymer runs; typically the run in the assembly is shorter
  than accurate reads indicate.
- Read-based variant calling against such an assembly would be swamped by these
  false indels. Module 11 quantifies them.
- Similar length errors occur in dinucleotide repeats, so these are counted
  too.
- Heterozygous calls in an inbred line mostly reflect mismapping in repeats,
  so they are never counted as errors.
- The expected errors are short length changes in runs that 150-bp reads span,
  so bcftools is adequate.

**Pipeline:**
1. **Reads:** the Illumina BAM, else the HiFi BAM.
2. **Callable region:** positions with depth in [0.5 × median, 2 × median],
   where the median depth over chromosomes comes from `samtools depth -a`
   sampled every 1 kb. Write `callable.bed`.
3. **Call**, per chromosome and per unplaced chunk in parallel:
   ```
   bcftools mpileup -f ASM.fa -a AD,DP -q 20 -Q 20 -d 500 [HiFi: -X pacbio-ccs] \
       -T callable.bed -r <region> reads.bam \
     | bcftools call -m -v --ploidy 2 -Oz
   ```
   Concatenate, then `bcftools norm -f ASM.fa -m -both` (left-aligned).
4. **Error calls:** GT 1/1, QUAL ≥ 30, inside `callable.bed`. GT 0/1 is counted
   in `m11_het_calls` only.
5. **Classify each hom-alt indel** against the assembly right after the anchor
   base:
   - **HP:** the inserted or deleted bases are all one base *b*; L is the
     assembly run length of *b* after the anchor. It counts when
     **max(L, L + change) ≥ 4**.
   - **STR2:** an even-length indel equal to whole copies of a 2-bp unit *u*
     (u[0] ≠ u[1]), with ≥ 2 tandem copies of *u* in the assembly.
   - **Other indel:** the rest.
   - **SNV/MNV:** counted separately.
6. **Denominators:**
   - `m11_callable_bp`
   - homopolymer runs ≥ 4 bp lying entirely inside `callable.bed`
7. **`errors.tsv`:** HP counts by
   - assembly run length: 4–8, 9–10, 11–12, 13–15, 16–20, > 20
   - change size: 1, 2, 3, 4, ≥ 5
   - base: A/T vs G/C
   - direction: ins vs del

   plus STR2 counts by unit and direction.
8. **`errors.bed.gz`:** one line per error: seq, start, end, class, base or
   unit, run length, change (e.g. `+2A`), QUAL, depth. **`hom_calls.vcf.gz`:**
   the filtered calls.
9. **Report:**
   - a stacked bar chart of HP errors by run-length bin;
   - text: insertions relative to the assembly mean the assembly run is too
     short;
   - cross-check (report only): the share of Merqury assembly-only k-mers
     within 20 bp of an error.

---

## 9. Resource requirements

Estimates for a 4.3 Gb assembly on 32 threads. **To be replaced by measured
values** before v1.0.

| Stage | Wall time | Peak RAM | Work disk |
|---|---|---|---|
| prep + M1 | < 30 min | 8 GB | 10 GB |
| M2 QUAST | 1–3 h | 32 GB | 10 GB |
| M4 tidk | < 30 min | 4 GB | < 1 GB |
| M5 | 1–2 h | 16 GB | 15 GB |
| M6 BUSCO | 3–8 h | 32–64 GB | 20 GB |
| M7 | ~1 h | 30–40 GB | 5 GB |
| map Illumina 30× | 4–8 h | 30 GB | 100–150 GB |
| map long reads 30× | 4–10 h | 30 GB | 100–200 GB |
| M8 | 3–6 h | 64–128 GB | 50–100 GB |
| M9 CRAQ | 4–12 h | tbd | 50 GB |
| M11 | 2–6 h | 8 GB | 5 GB |

| Run type | Time | RAM | Workdir |
|---|---|---|---|
| Full run with Illumina and long reads | ~1–2 days | 128 GB | ~500 GB |
| Assembly-only run (M1–M7) | ~6–12 h | 64 GB | ~60 GB |

---

## 10. Repository, build and release

- **Repository:** this repository, public, **GPL-3.0**. It reuses GPL-3.0 code
  and data from `kavonrtep/ont_genome_assembly_pipeline` and `kavonrtep/CARP`;
  both are credited in the README and in the vendored files' headers.
- **Never commit:**
  - unpublished data of any kind (assemblies, AGPs, reads, BAMs, excerpts);
  - results from real assemblies;
  - host names or internal paths.

  Validation against consortium data lives in an untracked `validation/`
  directory (listed in `.gitignore`) that takes its paths as arguments. Public
  reference data are downloaded at build time, not committed.
- **Layout:** `Singularity`, `envs/*.yaml` + `envs/*.lock`, `envs/dev.yaml`
  (development and unit tests outside the image), `pyproject.toml`,
  `src/asmqc/` (Python package: CLI, validation, parsers, classifiers,
  summary writer), `workflow/` (Snakemake; rules stay thin and call
  `python -m asmqc.<module>`), `templates/`, `tests/`, `docs/design/`,
  `README.md` (user guide), `CHANGELOG.md`, `LICENSE`.
- **Pinning:** every conda package exactly as in §5.1, built from committed
  explicit lock files. No unpinned dependency anywhere.
- **The def file:**
  1. create the environments from the locks;
  2. download and md5-verify the reference data;
  3. install the workflow;
  4. `%test` runs `asmqc test`.
- **Versioning:** semantic versioning; image `asmqc_<version>.sif`.
  `/opt/asmqc/VERSION.json` holds the version, git commit and lock-file
  hashes. Any change to a tool version, parameter or reference file is a minor
  or major bump (design principle E).
- **First task: release-path canary.**
  1. A `workflow_dispatch` workflow (`permissions: {contents: read,
     packages: write}`) installs Apptainer (e.g.
     `eWaterCycle/setup-apptainer`) and ORAS (`oras-project/setup-oras`).
  2. It builds a ~5 MB test SIF (`From: alpine`) and pushes it with
     `GITHUB_TOKEN` to `ghcr.io/peaterpan-pan-genome-project/asmqc-canary:test`.
  3. It pulls the image back and compares the sha256.
  - On failure, report the step and error to the maintainers; the cause is
    usually an organization setting (allowed actions, package creation).
  - On success, delete the canary package and keep the workflow as the template
    for `release.yml`.
- **Release:** build on a pushed `v*` tag. Publish the SIF as a GitHub release
  asset and to
  `oras://ghcr.io/peaterpan-pan-genome-project/asmqc/sif:<version>` (public
  package, linked to this repository), with its sha256. The first versioned
  release needs the maintainers' go-ahead.
- **README (user guide):**
  - installation (`singularity pull`, sha256 check);
  - three worked examples: assembly-only, +Illumina, +HiFi+ONT;
  - the required declarations and why they matter;
  - the resource table;
  - what to send back: `tar czf <label>_asmqc.tgz -C <outdir> <label>`;
  - how to read the report.

---

## 11. Testing and acceptance

### 11.1 Built-in smoke test (`asmqc test`, run in `%test`)

A synthetic data set under 20 MB, generated by a committed, seeded script
(`tests/make_testdata.py`). It contains:
- 7 "chromosomes" of 200–400 kb, `chr1`–`chr7`, with:
  - telomeric arrays: 12 capped arms, 1 arm with no array, 1 arm in the wrong
    orientation;
  - one interstitial array;
  - a 45S-like array (3 copies) and a 5S-like array (10 copies) built from the
    library;
  - a 20 kb plastid segment inside one chromosome.
- Unplaced sequences:
  - one 15 bp (→ ENA FAIL)
  - one with a trailing N
  - one 50 kb exact copy of a chromosome segment (→ `duplicate`)
  - one 30 kb plastid fragment (→ `organelle_scaffold`)
  - several of exactly 1,000 and 2,000 bp (round lengths)
- A matching AGP, in which some components are named `…_subseq_A:B` with
  round cut coordinates.
- Simulated paired short reads, ~20× from a "true" genome in which 30
  homopolymers of 9–12 bp are 1–2 nt longer than in the assembly, plus 5
  dinucleotide differences.

The BUSCO run is skipped (no real genes), but the lineage is asserted present:
`fabales_odb12.2`, 2026-05-13, 7,702 markers, and the md5.

**Assertions:**
- exact M1 counts (including the round-cut counts), `ena_rules = FAIL`
- the M4 arm statuses
- the M5 classes and rDNA copy counts
- the M7 classes
- in M11, **≥ 28 of the 30 planted HP errors and ≥ 4 of the 5 STR2 errors
  recovered**, with the correct direction and size

### 11.2 Validation on real assemblies

Before v1.0 the image is validated on consortium assemblies with known
properties, from the untracked `validation/` directory. That covers:
- integrity flags on an assembly with known ENA-blocking sequences, round
  lengths and an AGP name mismatch;
- a clean assembly (`ena_rules = PASS`);
- BUSCO sanity (7,702 markers, the expected range) and determinism (two runs
  identical);
- rDNA-only scaffolds recognised;
- telomere arm statuses on carefully curated assemblies (this also confirms
  the orientation convention);
- the Python contiguity numbers agree with QUAST;
- module 11 against an independent GATK call set on the same Illumina reads
  (**≥ 85 %** of the GATK hom-alt homopolymer calls recovered; total HP count
  within ±20 %);
- module 8 runs end to end;
- measured resources replace the §9 estimates.

The concrete data and expected values are kept privately by the maintainers.

### 11.3 Acceptance for v1.0

- `asmqc test` passes in the build.
- The §11.2 validation passes, or deviations are documented and accepted.
- The image runs with no network.
- `asmqc aggregate` works on the validation outputs.
- `report.html` opens offline and is complete.
- The README covers §10.

---

## 12. Decision log

| Date | Decision |
|---|---|
| 2026-10-01 | Modules 1, 2, 4, 5, 6, 7, 8, 9, 11; out of scope: genome size vs flow cytometry; centromere/satellite validation |
| 2026-10-01 | Reference data in the image; `singularity run` + Snakemake on a single machine; TSV/JSON + HTML; pea only |
| 2026-10-01 | M1: ENA rules as PASS/FAIL; round-number diagnostic report-only; no sample-sheet check |
| 2026-10-01 | mm2-plus everywhere; mm2-gb rejected |
| 2026-10-01 | M2: `chr1`–`chr7` required (or a mapping); N-split at ≥ 10; QUAST HTML not bundled |
| 2026-10-01 | M4: capped = band within 50 kb; wrong orientation diagnostic; unplaced telomere scaffolds listed, not plotted |
| 2026-10-01 | M5: plastid NC_014057.1, mito PP555264.1; ≥ 80 % at ≥ 95 %; organelle scaffold = WARNING; rDNA library from CARP 1.9.0 |
| 2026-10-01 | M6: BUSCO 6.1.0 + miniprot, `fabales_odb12.2`; derived unplaced and duplicate-placement numbers; `full_table.tsv` kept |
| 2026-10-01 | M7: one-chain rule (≥ 90 %, ≥ 99 % identity, MAPQ ≥ 20); scaffolds ≥ 1 kb |
| 2026-10-01 | M8: k = 21; `--reads-used-in-assembly` required; Illumina preferred over HiFi; run and flag below 20× |
| 2026-10-01 | M9: CRAQ only; ONT preset from the declared chemistry; AGP-junction cross-check |
| 2026-10-01 | M11: bcftools; homopolymer ≥ 4 bp; dinucleotide class; error BED kept |
| 2026-10-01 | Public repository, GPL-3.0; no unpublished data committed |
| 2026-10-01 | Python code as package `src/asmqc/`; thin Snakemake rules; `envs/dev.yaml` + pytest for development outside the image. Build order in `IMPLEMENTATION_PLAN.md` |
| 2026-10-01 | CRAQ version taken from `conda-meta` (its own output reports `1.0.9-alpha`) |
| 2026-10-01 | BUSCO runs with `--opt-out-run-stats` (no network at run time) |

## 13. Open points (to the maintainers before deciding)

1. A pinned version that cannot be installed, even in a separate environment.
2. CRAQ 1.10's exact output paths and `out_final.Report` format: show one
   parsed example.
3. The tidk forward/reverse orientation convention, confirmed on real data.
4. If M11 misses the 85 % target: the per-bin concordance first.
5. The canary outcome; the go-ahead for the first versioned release.
6. Any measured resource figure more than 2× the estimate.
7. BUSCO 6.1.0 must accept `--lineage_dataset fabales_odb12.2 --offline`.
