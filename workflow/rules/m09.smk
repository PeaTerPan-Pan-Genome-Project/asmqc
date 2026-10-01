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
        ngs=f"-ngs {W}/map/illumina.bam" if M09_NGS else "",
    threads: workflow.cores
    log:
        log("m09_craq"),
    benchmark:
        bench("m09", "craq")
    shell:
        "{params.pre}( mkdir -p {W}/m09/craq && cd {W}/m09/craq && rm -rf output"
        " && craq -g {input.fa} -sms {input.sms} {params.ngs} -t {threads} ) > {log} 2>&1"


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
