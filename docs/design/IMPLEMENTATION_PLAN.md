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

Milestone 6 (2026-10-01). M4, M5, M6, M7. On the synthetic data with the real
references every M4, M5 and M7 value equals `expected.json`; M6 runs (0 %
complete, no genes) and its lineage and predictor assertions pass. Resolved:
- §13.3: tidk `forward_repeat_number` counts the search string (TTTAGGG),
  `reverse_repeat_number` its reverse complement; `window` is the window end,
  capped at the sequence length. Confirmed on the planted arms; real curated
  assemblies (§11.2) remain the biological check.
- §13.7: BUSCO 6.1.0 runs `--offline` with `fabales_odb12.2` (25 s on 2 Mb).
  BUSCO 6 sends usage statistics unless `--opt-out-run-stats`; added to §8.6.
**[default]** choices:
- M4 distances are at window resolution: band start − 1, or length − band
  end. Dominance is a strict majority; a tie is `wrong_orientation`. A band
  > 50 kb from both ends is interstitial even when it is the band nearest an
  end (that arm is then `absent`). `m04_unplaced_with_telomere` counts
  sequences.
- M5 `m05_rdna45s_loci` / `m05_rdna5s_loci`: chromosomes carrying an array of
  the family, in chr1–chr7 order, then `unplaced` if any array lies on an
  unplaced sequence (e.g. `chr4;chr7;unplaced`). Array coordinates span the
  first to the last subunit hit. Largest organelle block = largest merged
  interval. `organelle_scaffolds.tsv` lists every unplaced sequence with an
  organelle hit; the `organelle` column says whether it passed 80 %.
- M6 `m06_complete_on_unplaced` counts Complete (single-copy and duplicated)
  BUSCO ids with at least one copy on an unplaced sequence.
- M7 total for `m07_total_minus_duplicate_bp` is the whole assembly. The
  §8.7 cross-check with M6 is computed in the report (milestone 8).
M4 `karyoplot.png`: `src/asmqc/karyoplot.py`, adapted from the vendored
`telomere_karyoplot.py` (§6, origin and changes in its header). Bands are
coloured by status (status palette, each colour with a legend label); absent
arms carry a marker.

Milestone 7 (2026-10-01). Shared `map_<readtype>` (`asmqc.mapping`: one
mm2plus run per file or R1/R2 pair with its own read group, sort, merge,
index) and `meryl_reads`; M8, M11, M9. Full run with Illumina 25× and HiFi
25× on the synthetic data: 75 s, all nine modules `ok`; M11 recovers 30/30
HP and 5/5 STR2 errors with correct direction and size and no other hom-alt
call. The test generator gained `--hifi-coverage` (development only; the
§11.1 data set is unchanged).
- `gawk` 5.4.1 added to `core`/`dev`: the callable BED is
  `samtools depth -a | gawk` per sequence in parallel (Python is too slow
  for per-base depth on 4 Gb). The median depth samples `samtools depth -a -b`
  at 1-kb positions over chr1–chr7.
- §13.2 (CRAQ 1.10 format), from a run on BAM input: results under
  `output/runAQI_out/`; `out_final.Report` has a `Genome` row and one row per
  sequence with `Covered.Rate`, `Low-confident.Rate`, `Avg.CRH`, `Avg.CSH`,
  `Avg.CRE(R-AQI)`, `Avg.CSE(S-AQI)`; CRE/CSE beds in
  `locER_out/out_final.CRE.bed` and `strER_out/out_final.CSE.bed`. CRAQ refuses
  an existing `output/`, so the rule declares a sentinel and clears it first.
  **There is no single AQI column**: `m09_aqi` is `NA` until the maintainers
  define it. CRAQ also reports CREs at sequence ends and at AGP gaps.
- `hom_calls.vcf.gz` holds the hom-alt error calls with `DP` only; header
  lines carrying local paths (`##reference=`, `##bcftools*`) are dropped.
  BUSCO's `short_summary.json` path parameters are reduced to their last
  component. Both files are part of the shared result.
- M8 on the synthetic data: k-mer peak at 20×, exactly the `low_coverage`
  threshold (not raised). The smoke test has no margin there.
**[default]**: M11 counts het calls at QUAL ≥ 30 like errors; `errors.bed.gz`
`change` is `+2A`/`-1T` for HP and `+2`/`-4` (bp) for STR2;
`m11_hp_ins_del_ratio` is `NA` without deletions.

Milestone 8 (2026-10-01). `report.py`, `plots.py`, `aggregate.py`,
`templates/`. `report.html` is written by the wrapper after the merge (a
report failure exits 2); it embeds every plot as a base64 PNG and has no
external request. **[default]** choices:
- Plots are matplotlib PNGs with one colour scheme (categorical slot 1 for a
  single series, slots 1–2 for A/T vs G/C, status colours in the karyoplot).
  The HTML follows light/dark; PNGs sit on a light card in both.
- Each module section lists all its `qc_summary.tsv` columns as key numbers,
  a fixed "how to read" paragraph and its plots. The §8.7 and §8.11
  cross-checks appear as notes in the M6 and M11 sections.
- The M2 cumulative plot uses contig lengths from the workdir
  (`m02/contig_lengths.txt`); the result directory does not carry them.
- `asmqc aggregate` also refuses mixed asmqc MAJOR.MINOR versions (principle
  E) and duplicate labels, besides differing headers (§7.2). It writes
  `combined/qc_summary.tsv` and `combined/report.html` (grouped table, M4
  karyoplots and M11 plots as small multiples, M6 internal stop % against
  M11 HP errors per Mb).

Milestone 9 (2026-10-01). `asmqc test` (`smoketest.py`), `Singularity`,
`build.sh`, `README.md`, `CHANGELOG.md`, `release.yml`. Local build with
Apptainer 1.4.5: 4 min, `asmqc_0.1.0.dev0.sif` 2.9 GB, `%test` passes.
Checked with the image:
- `asmqc test` in a network namespace with only a downed loopback: passes,
  51 s (no network at run time, §11.3).
- full run, all nine modules incl. BUSCO and CRAQ, offline: 70 s, every
  module `ok`; `qc_summary.tsv` identical to the dev-env run apart from the
  identity columns.
- tool versions as pinned; no local path in the result directory.
Findings:
- Snakemake writes a source cache under `$XDG_CACHE_HOME` (default
  `~/.cache`), read-only in `%test`. The runner now puts `XDG_CACHE_HOME`,
  `MPLCONFIGDIR` and `TMPDIR` under the workdir.
- Image size 2.9 GB exceeds the 2 GiB GitHub release-asset limit (§13.9);
  `release.yml` attaches the SIF only if it fits.
- `release.yml` and the canary are untested; they need a pushed tag.

Release path (2026-10-01). Canary run 36895113235 passed (§13.5): the
organisation allows the actions, `GITHUB_TOKEN` creates and pushes a ghcr.io
package, the pulled-back SIF has the same sha256. The `asmqc-canary` package
must be deleted by a user with the `delete:packages` scope.
`release.yml` builds as root, verifies the pushed image by pulling it back,
and marks a/b/rc/dev versions as pre-releases; first test with `v0.1.0rc1`.
Release run 36897267740 for `v0.1.0rc1` passed in 13 min 51 s: build as
root, `asmqc test` passed on the runner, image pushed to
`oras://ghcr.io/peaterpan-pan-genome-project/asmqc/sif:0.1.0rc1`
(sha256 `18dfbbad…9239`), pulled back and verified, pre-release created with
the sha256 asset and the pull command. The package is private after the
first push; anonymous pulls fail until its visibility is set to public.
Public packages had to be enabled at organisation level (Settings →
Packages → Package creation) before `asmqc/sif` could be made public. After
that (2026-10-01) an anonymous `apptainer pull` takes 2.5 min, the sha256
matches the release asset, and the pulled image passes `asmqc test` offline.

Real data, Cameor v2 (GCA_977071245.1), 2026-10-01. Full assembly (3.9 Gb),
modules 1, 4, 5, 7 with the `0.1.0rc1` image: 9 min on 14 threads (scan 2.2
min, M5 2.2 min, M7 4.1 min); ENA PASS. Findings:
- **Bug, fixed**: BUSCO 6.1.0 `full_table.tsv` has 10 columns for found
  BUSCOs (`OrthoDB url`, `Description`) and 2 for missing ones; the parser
  assumed 8 and M6 failed. The synthetic data have no genes, so only Missing
  rows were ever seen. The parser now takes column names from the
  `# Busco id` line. `0.1.0rc1` has the bug.
- **M5 rDNA, open**: every BLAST hit counts as a copy, so dispersed subunit
  fragments make "arrays" on all chromosomes; `m05_rdna45s_loci` and
  `m05_rdna5s_loci` list every chromosome. 22 45S "arrays" contain no 18S
  copy. Real arrays stand out: 45S on chr4 (16 copies) and chr7 (11), 5S on
  chr1 (403), chr3 (127), chr2 (45). A minimum hit length per copy and a
  minimum copy number per array are needed; to be decided.
- M4: 8 of 14 arms capped, 2 T2T; 17 of 19 interstitial arrays form one
  cluster on chr4 at 74.0-75.3 Mb.
Real-data fixture (`tests/realdata/`): ~50 Mb cut from Cameor v2 with
`make_cameor_fixture.py` (arrays with >= 10 copies, block ends moved off
1-kb multiples so no round AGP cut is an artefact of the fixture). Modules
1, 2, 4, 5, 6, 7 on it: 4.3 min (BUSCO 4 min), all ok; BUSCO C 3.1 %,
internal stops 4.2 %.
M5 rDNA decided (2026-10-02, SPEC §8.5): copies need ≥ 50 % subunit
coverage; arrays need ≥ 3 (45S) or ≥ 10 (5S) copies, smaller clusters are
fragments with their own columns. With ≥ 3 for 5S, hundreds of 3-4-copy
clusters of 5S-like sequence (within ~1 kb) still made 5S "arrays" on every
chromosome. Cameor v2 with the final rule: 45S loci chr3;chr4;chr7;unplaced
(6 arrays, 60 copies), 5S loci chr1;chr2;chr3;unplaced (8 arrays, 2,644
copies). Released as `0.1.0rc2` with the M6 fix.

M6 synteny (2026-10-02). Prototype on the phylogeny BUSCO records
(`fabales_odb12`, 22 pangenome assemblies, scratch only): 3 s for all
dotplots; shows chr1/chr5 and chr3/chr5 translocations relative to Caméor in
several accessions. Implemented with anchors from asmqc M6 on full Caméor v2
(`0.1.0rc2`, 64 threads: 9.5 min; C 99.5 %, 7,404 single-copy anchors on
chr1-chr7). On the 50 Mb fixture: 215 shared BUSCOs, every chromosome
matches its own Caméor chromosome, forward. §9 note: M6 on 3.9 Gb took 9.5
min on 64 threads, far below the 3-8 h estimate.

Large chromosomes (2026-10-02). The JI1006 run with `0.1.0rc3` failed at
`samtools index`: BAI and TBI store positions only up to 2^29
(536,870,912 bp); pea chromosomes reach ~650 Mb. Audit of every index:
- our BAMs (`mapping.py`): now `samtools index -c`; consumers
  (`samtools depth -r/-b`, `bcftools mpileup -R`) read CSI;
- M11 VCF: `bcftools index` was CSI by default, now explicit `-c`; no tabix;
- CRAQ 1.10: `bin/craq` dies unless `<bam>.bai` exists (existence test
  only); `runLR.sh`/`runSR.sh` run plain `samtools index` without `set -e`,
  so a failed BAI was silent and the index is never read (all later steps
  stream). Fix: inputs linked into `m09/craq/inputs/` with `.bai` -> `.csi`
  links, and `workflow/bin/craq_shim/samtools` (CRAQ rule only) turns
  `index` into `index -c`. Synthetic M9 results unchanged.
- No index involved: mm2-plus, faidx, seqkit, BLAST, tidk, QUAST, BUSCO,
  meryl/Merqury, `samtools depth | gawk`; the full 656 Mb Caméor chr5 went
  through M1, M4, M5, M6, M7 without trouble.
Test: `tests/unit/test_large_chromosomes.py` (a 700 Mb header, a read at
600 Mb: BAI fails, our CSI path and the shim work).
Released as `0.1.0rc4`.

### Explicit report and T2T rule (2026-10-03)

After the first full JI1006 run the report was made readable without the
specification: `src/asmqc/docs.py` holds a label, unit and definition for
every `qc_summary.tsv` column (a test enforces completeness), the module
explanations, flag meanings, a glossary and an output file guide. New
report-only tables: per chromosome, largest rDNA arrays with an rDNA
karyoplot, CSE positions, and M11 errors inside the coding exons of the
Complete single-copy BUSCOs (M6 now writes `busco_cds.bed.gz` from BUSCO's
per-gene miniprot GFF). On JI1006, 2 of 30,540 homopolymer errors fall in
11.3 Mb of BUSCO exons and 102 within 100 bp of one; the errors are in A/T
runs of ≥ 9 bp that coding sequence rarely contains. M4 T2T now also
requires no gap (SPEC §8.4). Released as `0.1.0rc5`; `0.1.0rc6` adds the
template row missed in rc5. `0.1.0` (2026-10-05) is the rc6 code after a
full JI1006 run with reads passed.

### CRAQ script replacements (2026-10-06, 0.1.2)

Issue #2: after 0.1.1, CRAQ was 12.4 h of a 14.7 h JI1006 run (mean load
0.81, max_rss 116 GB): about 6.2 h for the concurrent read passes and 6.2 h
for the AQI stage. Six CRAQ scripts load per-base depth tables (75 GB) into
Perl hashes or regex every line; `workflow/craq_patch/src/` replaces them with
streaming versions (two from the issue, four new). CRAQ is copied into
`work/m09/craq/sw/` and run with the replacements, so the env stays pristine
and development and image behave alike. `tests/unit/test_craq_patch.py`
compares each original and replacement byte for byte on seeded tables
(several sequences, one-line and missing sequences, sites at ends, zero runs,
seven get_ER parameter sets, run lengths 149/150); mutations of the
replacements are caught. On a 30 M-line table: 1.2–3.3× faster, hash-based
scripts at under 5 % of the memory. M9 on synthetic data: identical results to
0.1.1; CRAQ intermediates identical up to row order (CRAQ's own
`caculate_breakpoint_depth.pl` writes in hash order).

### CRAQ without depth tables (2026-10-06, 0.1.3)

Stage 1 and 2 of the no-table plan. Readers of `LR_sort.depth`: effective
size (+ `LRout/Nonmap.loc`), `synthesize_LRbkdep`, `synthesize_clipDIcov`,
`get_ER` (AQI). Readers of `SR_sort.depth`: effective size (+
`SRout/Nonmap.loc`), `synthesize_SRbkdep`, `search_dep0`,
`LRcoverRate_srdep_filter` (AQI). Clip and indel sites need only the BAM, so
the patched drivers extract them first and stream one depth table per BAM to
all within-pass readers (`asmqc_fanout.pl`: block copy to child pipes; fails
if a child fails, unlike `tee` with FIFOs, which hangs, or process
substitution, which loses exit codes). The AQI lookups cross read types and
run after both passes; they use region queries. Checked on synthetic BAMs:
region output equals the full table (order, read-free stretches, regions past
the end, unsorted and overlapping BED); `samtools depth -a -Q 20` equals
CRAQ's `samtools view -q 20 | samtools depth -a -`; sequences without reads
(also: only MAPQ < 20 reads) are omitted by both. Region output is still
limited to the sequences seen in the stream. Tests: contract test of the
region helper (mutations caught), original window scripts on full vs region
input, fan-out failure handling, effect-size scripts incl. `Nonmap.loc`.
M9 on synthetic data: every CRAQ file identical to unpatched CRAQ 1.10 (up to
CRAQ's own hash-order rows); only the depth tables are gone. 30 M-line table:
long-read stream readers 18.7 s together vs 100.8 s for the originals in
sequence. Real-data validation: compare a JI1006 run with 0.1.2.

### Concurrent mapping and segmented CRAQ passes (2026-10-06, 0.1.4)

Mapping: `map_reads` threads come from `asmqc.plan.map_threads` (largest
remainder of cores by input bytes per read type), so HiFi and Illumina map
at the same time; in 0.1.1 they ran one after the other (4,086 s + 2,912 s at
a mean load of 46 and 51 of 64 cores). BAM records are identical with any
thread count (checked on synthetic data).

CRAQ: the clip scan (`caculate_breakpoint_depth.pl`) and indel scan
(`caculate_clipDI_cov.pl`) keep no state across sequences and print in Perl
hash order, except that the indel scan prints all deletions before all
insertions (a later step keys on position, so an insertion wins over a
deletion at the same site). The patched drivers run the BAM filter, the
scans and the depth stream per segment of whole sequences (contiguous in
header order, so concatenated parts follow table order), CRAQ's `-t`
segments at a time. Region access uses `samtools view -M -L segment.bed` on
the indexed BAM (CRAQ links the input BAM without its index; the drivers
resolve the link). Joins: `cat` for table-ordered outputs and the clip scan,
`asmqc_merge.pl` for the indel scan (D then I), the effective-size counts and
`search_dep0.pl` (string-sorted sequences). `tests/unit/test_craq_patch.py`
compares segmented with genome-wide results for 1, 2, 3 and 8 segments on a
BAM with clipped and indel-bearing read clusters. M9 and M11 on synthetic
data are identical to 0.1.3, and so are CRAQ's filtered BAMs (records).
