# Module 9: read-back structural validation with CRAQ (SPEC §8.9).

M09_LONG = PLAN["m09"]["read_type"]
M09_NGS = PLAN["m09"]["short_reads"]


rule m09_craq:
    input:
        fa=W / "prep" / "asm.fa",
        sms=W / "map" / f"{M09_LONG}.bam",
        ngs=[W / "map" / "illumina.bam"] if M09_NGS else [],
    output:
        # CRAQ refuses an existing output/ dir, which Snakemake would create for
        # a declared output inside it; a sentinel outside it is declared instead
        touch(W / "m09" / "craq.done"),
    params:
        pre=env("craq", "m09"),
        ngs="-ngs inputs/illumina.bam" if M09_NGS else "",
        shim=Path(workflow.basedir) / "bin" / "craq_shim",
        patch=Path(workflow.basedir) / "craq_patch" / "src",
        parallel=1 if M09_NGS else 0,
    threads: M09_THREADS
    resources:
        mem_mb=BUDGET["craq_mem_mb"],
    # longest job: start it as soon as its BAMs exist, before side modules
    priority: 50
    log:
        log("m09_craq"),
    benchmark:
        bench("m09", "craq")
    shell:
        # CRAQ insists on <bam>.bai and indexes its filtered BAMs with plain
        # "samtools index" (BAI, fails beyond 2^29 bp). Inputs are linked into
        # its run dir with .bai -> .csi links (CRAQ only checks the file
        # exists), and a samtools shim on PATH turns "index" into "index -c".
        # With short reads, the bash shim runs CRAQ's long- and short-read
        # passes concurrently (see workflow/bin/craq_shim/bash).
        # CRAQ runs from a copy of its bin/ and src/ with the streaming
        # replacements in workflow/craq_patch/src (same output, a fraction of
        # the memory and time); the env itself is not modified.
        "{params.pre}( mkdir -p {W}/m09/craq/inputs && cd {W}/m09/craq && rm -rf output lr.status sw"
        " && c=$(dirname $(readlink -f $(command -v craq))) && mkdir sw"
        " && cp -rL $c sw/bin && cp -rL $c/../src sw/src && cp {params.patch}/*.pl sw/src/"
        " && for b in {input.sms} {input.ngs}; do n=$(basename $b);"
        " ln -sf $b inputs/$n && ln -sf $b.csi inputs/$n.csi && ln -sf $b.csi inputs/$n.bai; done"
        " && ASMQC_CRAQ_PARALLEL={params.parallel} ASMQC_CRAQ_LR_STATUS={W}/m09/craq/lr.status"
        " PATH={params.shim}:$PATH perl {W}/m09/craq/sw/bin/craq -g {input.fa} -sms inputs/$(basename {input.sms})"
        " {params.ngs} -t {threads} ) > {log} 2>&1"


rule m09:
    input:
        report=rules.m09_craq.output,
        bam=W / "map" / f"{M09_LONG}.bam",
        fai=W / "prep" / "asm.fa.fai",
        agp=[W / "prep" / "asm.agp"] if config["agp"] else [],
    output:
        W / "m09" / "summary.json",
    params:
        pre=env("core", "m09_summary"),
        agp=f"--agp {W}/prep/asm.agp" if config["agp"] else "",
        ngs="yes" if M09_NGS else "no",
        out=R / "m09_craq",
    log:
        log("m09"),
    benchmark:
        bench("m09", "m09")
    shell:
        "{params.pre}{PY} -m asmqc.m09_craq --craq-dir {W}/m09/craq --bam {input.bam}"
        " --fai {input.fai} {params.agp} --long-read-type {M09_LONG}"
        " --short-reads-used {params.ngs} --outdir {params.out} --work {W}/m09 > {log} 2>&1"
