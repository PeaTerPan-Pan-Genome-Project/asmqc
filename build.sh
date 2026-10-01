#!/usr/bin/env bash
# Build asmqc_<version>.sif (SPEC §10). Writes VERSION.json (version, git
# commit, lock-file sha256), then runs apptainer build; the image's %test
# runs `asmqc test`. Extra arguments go to apptainer build (e.g. --fakeroot).
#
# Usage: ./build.sh [apptainer build options]
set -euo pipefail
cd "$(dirname "$0")"

version=$(python3 -c 'import tomllib; print(tomllib.load(open("pyproject.toml", "rb"))["project"]["version"])')
commit=$(git rev-parse HEAD)
if [[ -n $(git status --porcelain --untracked-files=no) ]]; then
    commit="${commit}-dirty"
fi

python3 - "$version" "$commit" <<'EOF'
import hashlib, json, pathlib, sys
version, commit = sys.argv[1:3]
locks = {p.stem: hashlib.sha256(p.read_bytes()).hexdigest()
         for p in sorted(pathlib.Path("envs").glob("*.lock"))}
pathlib.Path("VERSION.json").write_text(json.dumps(
    {"version": version, "git_commit": commit, "lockfile_sha256": locks}, indent=1) + "\n")
EOF

apptainer build "$@" "asmqc_${version}.sif" Singularity
sha256sum "asmqc_${version}.sif" | tee "asmqc_${version}.sif.sha256"
