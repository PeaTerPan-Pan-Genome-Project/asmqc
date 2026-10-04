# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versions
follow semantic versioning; results are comparable only within one
`MAJOR.MINOR` version.

## [Unreleased]

### Changed

- Run time with long reads: CRAQ (M9) gets at most 8 threads and half of
  `--mem-gb`; the other modules use the rest and run while CRAQ runs. With
  short reads, CRAQ's long- and short-read passes run concurrently. No
  result changes.
- `m09_craq/out_final.Report` lists sequences in assembly order (was Perl
  hash order, different between runs).

### Added

- Per-rule benchmarks (wall time, peak memory, CPU load) in
  `logs/benchmarks/`.

## [0.1.0rc6] - 2026-10-04

Sixth release candidate. Not validated against SPEC §11.2; results are not
for comparison.

### Fixed

- Report: the row "Errors within 100 bp of a BUSCO exon" (M11) was missing
  from the template in 0.1.0rc5.

## [0.1.0rc5] - 2026-10-03

Fifth release candidate. Not validated against SPEC §11.2; results are not
for comparison.

### Changed

- M4: a chromosome is T2T only if both arms are capped and it has no gap
  (no N-run ≥ 10 bp). `telomeres.tsv` gains `chromosome_gaps` and `t2t`.
- Report: every value has a plain label, unit and definition (tooltip and
  table), with its `qc_summary.tsv` column name; module explanations
  rewritten, in detail for M8, M9 and M11. Combined report column headers
  carry the same labels and definitions.

### Added

- M6 writes `busco_cds.bed.gz`, the coding exons of the Complete
  single-copy BUSCOs.
- Report: at-a-glance values, per-chromosome table, rDNA karyoplot and
  largest arrays (M5), CSE list (M9), homopolymer errors in BUSCO coding
  exons with frameshift count, errors near exons and affected genes (M11),
  flag-code meanings, glossary and output file guide.
- README: measured wall times from the first full runs.

## [0.1.0rc4] - 2026-10-02

Fourth release candidate. Not validated against SPEC §11.2; results are not
for comparison.

### Fixed

- Read mapping failed on pea chromosomes longer than 536,870,912 bp (2^29):
  BAMs were indexed as BAI. Every BAM and VCF index is now CSI. CRAQ (M9)
  gets `.bai` links to the CSI index and a samtools shim that makes its
  internal indexing CSI too.

## [0.1.0rc3] - 2026-10-02

Third release candidate. Not validated against SPEC §11.2; results are not
for comparison.

### Added

- M6 synteny with Caméor v2 (report only): a dotplot from shared Complete
  single-copy BUSCOs, `synteny.tsv` (best-matching Caméor chromosome and
  orientation per chromosome) and `synteny_points.tsv`; small multiples in the
  combined report. The anchors are a committed table; the Caméor sequence is
  not in the image.

## [0.1.0rc2] - 2026-10-02

Second release candidate, after the first runs on a real assembly (Caméor
v2). Not validated against SPEC §11.2; results are not for comparison.

### Changed

- M5 rDNA: a BLAST hit counts as a copy only if it covers at least 50 % of
  the library subunit; a cluster is an array with at least 3 (45S) or 10 (5S)
  copies. Smaller clusters are counted in the new columns
  `m05_rdna45s_fragments` and `m05_rdna5s_fragments`. `rdna_arrays.tsv` gains
  a `class` column.

### Fixed

- M6 failed on assemblies with genes: BUSCO 6.1.0 writes 10 columns in
  `full_table.tsv` for found BUSCOs.

## [0.1.0rc1] - 2026-10-01

Release candidate to test the release path. Not validated on real
assemblies (SPEC §11.2); results are not for comparison.

### Added

- Modules M1, M2, M4, M5, M6, M7, M8, M9 and M11 as specified in
  `docs/design/SPEC.md`.
- `asmqc run`, `asmqc version`, `asmqc aggregate` and `asmqc test`.
- Self-contained `report.html` per assembly and a combined report.
- Singularity definition with environments from explicit lock files and
  md5-verified reference data.
