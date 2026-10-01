Bootstrap: docker
From: mambaorg/micromamba:2.9.0

# asmqc image (SPEC §5.1, §10). Build with ./build.sh, which writes
# VERSION.json first. Every environment comes from its committed explicit
# lock; reference data are downloaded and md5-verified here, and the build
# fails on any mismatch.

%labels
    org.opencontainers.image.title asmqc
    org.opencontainers.image.description Comparable QC of chromosome-level pea genome assemblies
    org.opencontainers.image.licenses GPL-3.0-only
    org.opencontainers.image.source https://github.com/PeaTerPan-Pan-Genome-Project/asmqc

%files
    VERSION.json /opt/asmqc/VERSION.json
    LICENSE /opt/asmqc/LICENSE
    envs /opt/asmqc/envs
    src /opt/asmqc/src
    workflow /opt/asmqc/workflow
    templates /opt/asmqc/templates
    refs/refs.tsv /opt/asmqc/refs/refs.tsv
    refs/fetch_refs.py /opt/asmqc/refs/fetch_refs.py
    tests/make_testdata.py /opt/asmqc/tests/make_testdata.py

%post
    set -eu
    for env in core quast busco merqury craq; do
        micromamba create -y -q -p /opt/envs/$env -f /opt/asmqc/envs/$env.lock
    done
    micromamba clean -a -y -q

    /opt/envs/core/bin/python /opt/asmqc/refs/fetch_refs.py /opt/asmqc/refs --drop-archives

    mkdir -p /opt/asmqc/bin
    printf '#!/bin/sh\nexec /opt/envs/core/bin/python -m asmqc.cli "$@"\n' > /opt/asmqc/bin/asmqc
    chmod 755 /opt/asmqc/bin/asmqc
    find /opt/asmqc -name __pycache__ -prune -exec rm -rf {} +
    chmod -R a+rX /opt/asmqc /opt/envs

%environment
    export ASMQC_HOME=/opt/asmqc
    export ASMQC_ENV_ROOT=/opt/envs
    export ASMQC_REFS=/opt/asmqc/refs
    export PYTHONPATH=/opt/asmqc/src
    export PYTHONNOUSERSITE=1
    export PATH=/opt/asmqc/bin:$PATH

%runscript
    exec /opt/asmqc/bin/asmqc "$@"

%test
    /opt/asmqc/bin/asmqc test
