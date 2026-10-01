# asmqc: implementation plan

Build order, code layout and test strategy for v1.0. The behaviour is defined
in [SPEC.md](SPEC.md); this document defines how and in what order it is
built. Section references (§) point to SPEC.md.

---

## 1. Code layout

| Path | Contents |
|---|---|
| `pyproject.toml` | package metadata; `asmqc` console entry point |
| `src/asmqc/cli.py` | subcommands `run`, `version`, `aggregate`, `test` |
| `src/asmqc/validate.py` | input validation (§4.2) |
| `src/asmqc/plan.py` | module switching and read-type selection (§4.3); writes the Snakemake config |
| `src/asmqc/schema.py` | ordered `qc_summary.tsv` columns with type and format (§7.2) |
| `src/asmqc/flags.py` | flag codes and severities (§7.5) |
| `src/asmqc/summary.py` | merges module results into `qc_summary.tsv`, `flags.tsv`, `run_manifest.json` |
| `src/asmqc/mNN_*.py` | one module each: parsers, classifiers, `main()` called by the rule |
| `workflow/Snakefile`, `workflow/rules/*.smk` | thin rules; shell out to tools and `python -m asmqc.mNN_*` |
| `templates/` | Jinja2 templates for `report.html` and the combined report |
| `envs/*.yaml`, `envs/*.lock` | image environments (§5.1) and explicit locks |
| `envs/dev.yaml` | development environment (python, pytest, snakemake, core tools) |
| `tests/unit/` | parser and classifier tests on small fixtures |
| `tests/make_testdata.py` | seeded synthetic data set (§11.1) |

## 2. Cross-cutting design

- **Schema as code.** `schema.py` is the single source of the column order and
  formatting: percentages 2 dp, counts integer, rates 4 significant digits,
  enums checked. A unit test compares it with the §7.2 list.
- **Module output contract.** Each module writes its files under
  `<label>/mNN_*/`, plus `work/mNN/summary.json` (its own columns) and
  `work/mNN/flags.tsv`.
- **Merging outside Snakemake.** The wrapper merges after Snakemake returns.
  With `--keep-going` a failed module would block a merge rule; merging in the
  wrapper guarantees `qc_summary.tsv`, `flags.tsv`, `ena_rules` and exit code
  2 even when a module failed (§4.4).
- **Planning in the wrapper.** `asmqc run` validates, decides each module's
  status (`run`, `skipped_no_input`, `skipped_by_user`) and read types, and
  writes a Snakemake config YAML. `--dry-run` prints the plan and stops.
- **Environment switching.** Rules prefix
  `PATH=$ASMQC_ENV_ROOT/<env>/bin:$PATH`. `ASMQC_ENV_ROOT=/opt/envs` in the
  image; in development every env name resolves to the dev environment.
  `TMPDIR` is set per rule under the workdir.
- **Pure parsers.** Parsers and classifiers take files or records and return
  data structures, with no tool calls, so they are unit-testable. Fixtures are
  synthetic or derived from public data only.
- **Determinism.** Thread-order-dependent tool output is sorted before
  parsing; no timestamps outside the manifest; fixed seeds.

## 3. Milestones

Development is local, with commits directly on `main`; pushing is a separate decision.

| # | Milestone | Content |
|---|---|---|
| 0 | Scaffold | `pyproject.toml`, `src/asmqc/` skeleton, `envs/dev.yaml`, `.gitignore` (`validation/`, `*.sif`, work dirs), GPL-3.0 `LICENSE`, pytest CI workflow |
| 1 | Environments and locks | `envs/{core,quast,busco,merqury,craq}.yaml` with exact pins; `envs/make_locks.sh` writes explicit locks (`conda list --explicit --md5`). Confirms every pin resolves (§13.1) and the mm2plus binary name before code depends on them |
| 2 | Release canary | `.github/workflows/canary.yml` per §10. Pushing and triggering require maintainer approval; outcome reported (§13.5) |
| 3 | Synthetic test data | `tests/make_testdata.py`: assembly, AGP, Illumina reads and `expected.json` with the §11.1 truth. Built early because every module asserts against it |
| 4 | Core runner | CLI, validation, `prep`, Snakefile with module switching, summary/flags/manifest merge, exit codes 0/1/2 |
| 5 | M1, M2 | integrity, ENA rules, AGP check, round-number diagnostic, N-split contiguity; QUAST run and agreement check |
| 6 | Assembly-only modules | M4 (tidk bands, arm classification, karyoplot from the vendored script), M5 (organelle PAF, rDNA copies and arrays), M7 (one-chain classification), M6 (BUSCO, lineage assertion, derived values; skipped on test data) |
| 7 | Read-based modules | shared `map_<readtype>` and `meryl_reads`; M8 (Merqury, coverage peak); M11 (callable BED, parallel bcftools, HP/STR2 classifier); M9 (CRAQ; parser waits for §13.2) |
| 8 | Report and aggregate | `report.html` with embedded PNGs and per-module plots; `asmqc aggregate` with header check and combined report |
| 9 | Image and release | `Singularity` def: envs from locks, md5-verified references, `VERSION.json`, `%test` = `asmqc test`; `release.yml` on `v*` tags |
| 10 | Validation | runs on consortium data from the untracked `validation/` directory (§11.2); measured resources replace §9; v1.0 acceptance (§11.3) |

## 4. Testing

- **Unit:** per-module parsers and classifiers on fixtures (`tests/unit/`).
- **Schema:** column list against §7.2.
- **Integration:** `asmqc run` on the synthetic data in the dev environment
  for modules whose tools it contains.
- **Smoke:** full `asmqc test` inside the SIF (§11.1).
- **Determinism:** two runs on the same input give identical
  `qc_summary.tsv`.

## 5. Risks and unknowns

- Pinned versions may not resolve, alone or together (§13.1). Milestone 1
  surfaces this first.
- CRAQ output paths and format (§13.2) and the tidk orientation convention
  (§13.3) need real output.
- BUSCO 6.1.0 offline with `fabales_odb12.2` (§13.7).
- Merqury and its R dependencies inside the image.
- `%test` must run offline, on < 20 MB of data, in acceptable build time.
- Read simulation for the test data: a seeded, pinned simulator or a small
  Python implementation.

## 6. Findings

Milestone 1 (2026-10-01). All five image environments install from their
locks; every pinned tool version in §5.1 exists on bioconda.

- `pandas` is pinned to 2.3.3 in `core` and `dev`: snakemake 9.27.0 requires
  pandas < 3.
- `quast` env pins Python 3.11.16: QUAST 5.3.0 imports `distutils`, removed in
  Python 3.12, and fails under the otherwise-resolved Python 3.13.
- mm2plus binary is `mm2plus`, a dispatcher that launches
  `mm2plus.{avx2,avx512,…}` and prints the choice to stderr. Version parsing
  must take the last stdout line.
- CRAQ 1.10 (bioconda `craq-1.10`) prints `CRAQ Version: 1.0.9-alpha` in its
  usage text; upstream did not update the string. "Tool versions are read
  from the tools at run time" (§7.3) cannot hold for CRAQ: read the version
  from the environment's `conda-meta` instead (decided 2026-10-01; SPEC §7.3).
- `merqury.sh` needs `MERQURY=/opt/envs/merqury/share/merqury`, normally set by
  conda activation. Rules do not activate (§5.2), so the rule sets it.
- `craq` and `quast` envs pull minimap2 2.31 as a dependency. asmqc never calls
  it: CRAQ receives the mm2-plus BAMs and QUAST runs without a reference.

Milestone 3 (2026-10-01). `tests/make_testdata.py` writes `asm.fa`,
`asm.agp`, paired reads (25×, 157,861 pairs) and `expected.json`; 17 MB,
17 s, byte-identical per seed. Design choices beyond §11.1:
- Gap positions are not multiples of 1,000, so the round AGP cuts are only
  the planted ones: 8 of 28 (subseq `A − 1`, subseq `B`,
  `component_beg − 1`, five round-length unplaced).
- chr6 carries a 20 N run (a gap for `nsplit10`, not in the AGP) and a 5 N run
  (below the gap threshold), so `m02_contig_n50_nsplit10` differs from the AGP
  value.
- 25× read coverage instead of ~20×: k-mer coverage at k = 21 is ~22×, above
  the 20× `low_coverage` threshold.
- Unplaced names are descriptive (`unplaced_dup`, `unplaced_plastid`, …).
- Prototype M11 (mm2plus `-ax sr`, the §8.11 bcftools pipeline) on data built
  with random stand-in references: 30/30 HP and 5/5 STR2 errors recovered as
  hom-alt, no other hom-alt calls.

With the real references (md5s match §6): BLAST of the rDNA library finds
only the planted arrays (merged copies: chr4 3 × 18S/5.8S/25S, chr1 10 × 5S);
the plastid segments give no rDNA hits; mm2plus finds the plastid segments at
identity 1.0 and no mito hits.

Milestone 4 (2026-10-01). Core runner in `cli.py`, `validate.py`, `plan.py`,
`runner.py`, `summary.py`, `schema.py`, `flags.py`, `params.py`, `tools.py`,
`prep.py`, `workflow/Snakefile`, `workflow/rules/prep.smk`. Points the spec
leaves open, settled as **[default]**:
- Number format `dec2` (2 decimals) for values that are neither percentages,
  counts nor rates: `m08_kmer_coverage`, `m08_qv`, `m09_coverage`, `m09_aqi`,
  `m09_r_aqi`, `m09_s_aqi`.
- `ena_rules` is `NA` when M1 did not finish; PASS would be unfounded.
- `run_date` differs between runs, so the determinism check (§5.3) compares
  every column except `run_date`.
- Manifest `inputs.assembly.md5` is the decompressed md5 (= `assembly_md5`);
  `bytes` is the size of the file as given.
- `command_line` in the manifest reduces path arguments to file names, like
  `inputs` (§7.3); `--manifest-full-paths` keeps them.
- Manifest extras: `snakemake_exit_code`, `wall_seconds.total`.
- Cleanup removes only the entries the workflow creates in the workdir, never
  a user-supplied `--workdir` wholesale.
- A module decides `skipped_no_input` itself when it depends on assembly
  content (M7: no unplaced ≥ 1 kb), via `status` in its `summary.json`.
- A planned module without a rule file (not yet implemented) ends as
  `failed`; the run exits 2.
- Tool versions run each tool with its env's `bin/` first on PATH: scripts
  such as `busco` use `#!/usr/bin/env python`.

Milestone 5 (2026-10-01). M1 and M2, with a shared `scan` stage (one pass over
the FASTA; per-sequence composition, N-runs ≥ 10, file facts) that both read.
Scan of a 500 Mb sequence: 18 s, 2 GB RSS. On the synthetic data every M1
value equals `expected.json`; M2 agrees with QUAST. **[default]** choices:
- `seq_lt_200bp` is raised for 20–199 bp only (a < 20 bp sequence carries
  `seq_lt_20bp`); `m01_n_lt200bp` counts every sequence < 200 bp. An empty
  sequence carries `empty_sequence` only.
- `terminal_n`: one flag per affected end; `m01_n_terminal_n` counts
  sequences.
- `m01_n_duplicate_names` counts distinct names that occur more than once.
- An empty name or one with non-printable characters raises `invalid_char`
  with the name as `seq_id`; §7.5 has no separate code.
- `m01_n_invalid_chars` and `m01_n_iupac` count characters, not sequences.
- `m02_n_gaps` / `m02_gap_bp` follow the contig method: AGP gap lines when
  `agp`, N-runs ≥ 10 bp when `nsplit10`.
- QUAST agreement compares N50, L50 and N90 of the scaffolds (QUAST column 1)
  and of the N-split contigs (`_broken` column). QUAST's `# contigs` for the
  broken assembly is not compared: it differs from its own `>= 0 bp` count.
