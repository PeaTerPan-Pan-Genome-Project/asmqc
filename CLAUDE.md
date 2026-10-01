# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Status

Pre-implementation. The repository holds only `SPEC.md`, the v1.0 design
specification. There are no build, lint or test commands yet; add them here
once they exist. `SPEC.md` is the design document: read the relevant section
before changing anything, and stop and ask if a change diverges from it.
Points marked **[default]** may change; §13 lists open points that go to the
maintainers, not to be decided unilaterally.

## What it is

asmqc is one Singularity image that runs a fixed QC suite on chromosome-level
pea assemblies (`chr1`–`chr7` + unplaced) for a pangenome consortium. The
purpose is cross-group comparability: identical tools, versions, parameters,
reference data and output schema. It measures; it never modifies an assembly.

## Planned layout (SPEC §10)

`Singularity` (def file), `envs/*.yaml` + `envs/*.lock` (explicit lock files),
`workflow/` (Snakemake), `scripts/`, `templates/` (HTML report), `tests/`
(incl. seeded `tests/make_testdata.py`), `README.md`, `CHANGELOG.md`, `LICENSE`
(GPL-3.0). Untracked `validation/` (in `.gitignore`) holds runs against
consortium data.

## Architecture (SPEC §5)

- `asmqc` wrapper (`/opt/asmqc/bin/asmqc`, the `%runscript`) with subcommands
  `run`, `version`, `aggregate`, `test`. `run` validates inputs (exit 1 on
  failure, also under `--dry-run`), then calls Snakemake with `--cores`,
  `--resources mem_mb`, `--rerun-incomplete`, `--keep-going`. Single machine,
  no scheduler integration.
- Five separate micromamba envs under `/opt/envs/` (`core`, `quast`, `busco`,
  `merqury`, `craq`). Each rule puts its env's `bin/` on PATH explicitly; no
  conda activation at run time.
- Shared stages computed once and reused by modules: `prep` (decompress,
  index, rename via `--chromosomes`), `map_<readtype>` (one mm2-plus BAM per
  read type), `meryl_reads` (k = 21).
- Read-type selection is fixed: M8/M11 use Illumina else HiFi (never ONT);
  M9 uses HiFi else ONT, plus Illumina as `-ngs`.
- Module numbers 1, 2, 4, 5, 6, 7, 8, 9, 11. Numbers 3 and 10 are
  intentionally unused; the numbering is baked into `mNN_` column prefixes and
  output directory names.

## Invariants that span many files

- **Output contract (§7) is the interface.** `qc_summary.tsv` has a fixed
  column order; `asmqc aggregate` refuses differing headers. Adding, removing
  or reordering a column, or changing any tool version, parameter or reference
  file, is a minor/major version bump (results are comparable only within one
  `MAJOR.MINOR`).
- **Determinism:** same inputs + same image → identical `qc_summary.tsv`.
  Fixed seeds; sort thread-order-dependent tool output before parsing; wall
  times go only in `run_manifest.json`.
- Every parameter named in §8 must appear under `parameters` in the manifest.
  Tool versions are read from the tools at run time.
- Exit codes: 0 all applicable modules ran (an ENA FAIL is still 0); 1
  validation failure; 2 a module failed (other results still written, its
  columns `NA`, `mNN_status = failed`).
- `ena_rules = FAIL` iff at least one `ENA_BLOCKING` flag.
- All temporaries under `--workdir` (set `TMPDIR` per rule), never `/tmp`.
- No network at run time (BUSCO `--offline`). Reference data are downloaded
  and md5-verified at image build time; the build fails on mismatch.
- mm2-plus (`mm2plus`) for every minimap2-type step; AVX2 required.

## Repository rules

- Never commit unpublished data (assemblies, AGPs, reads, BAMs, excerpts),
  results from real assemblies, host names or internal paths.
- No unpinned dependency anywhere; envs are built from committed explicit lock
  files.
- Vendored GPL-3.0 code/data (`kavonrtep/ont_genome_assembly_pipeline`
  karyoplot script, `kavonrtep/CARP` rDNA library) must carry their origin in
  the file header and be credited in the README.
- First implementation task per §10: the release-path canary GitHub workflow
  (build a tiny SIF, push to ghcr.io via ORAS, pull back, compare sha256).
