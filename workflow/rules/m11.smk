# Module 11: homopolymer and short-STR errors (SPEC §8.11).

M11_ARGS = (
    f"--bam {W}/map/{PLAN['m11']['read_type']}.bam --fa {W}/prep/asm.fa"
    f" --fai {W}/prep/asm.fa.fai --read-type {PLAN['m11']['read_type']}"
    f" --reads-used-in-assembly {config['reads_used_in_assembly']} --work {W}/m11"
)


rule m11_callable:
    input:
        bam=W / "map" / f"{PLAN['m11']['read_type']}.bam",
        fai=W / "prep" / "asm.fa.fai",
    output:
        W / "m11" / "callable.bed",
    params:
        pre=env("core", "m11_callable"),
    threads: SIDE_THREADS
    log:
        log("m11_callable"),
    benchmark:
        bench("m11", "callable")
    shell:
        "{params.pre}{PY} -m asmqc.m11_homopolymer callable {M11_ARGS}"
        " --threads {threads} > {log} 2>&1"


rule m11_call:
    input:
        rules.m11_callable.output,
    output:
        W / "m11" / "calls.norm.vcf.gz",
    params:
        pre=env("core", "m11_call"),
    threads: SIDE_THREADS
    log:
        log("m11_call"),
    benchmark:
        bench("m11", "call")
    shell:
        "{params.pre}{PY} -m asmqc.m11_homopolymer call {M11_ARGS}"
        " --threads {threads} > {log} 2>&1"


rule m11:
    input:
        rules.m11_call.output,
    output:
        W / "m11" / "summary.json",
    params:
        pre=env("core", "m11"),
        out=R / "m11_homopolymer",
    log:
        log("m11"),
    benchmark:
        bench("m11", "m11")
    shell:
        "{params.pre}{PY} -m asmqc.m11_homopolymer classify {M11_ARGS}"
        " --outdir {params.out} > {log} 2>&1"
