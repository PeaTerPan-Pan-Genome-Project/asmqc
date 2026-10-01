# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versions
follow semantic versioning; results are comparable only within one
`MAJOR.MINOR` version.

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
