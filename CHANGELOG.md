# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versions
follow semantic versioning; results are comparable only within one
`MAJOR.MINOR` version.

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
