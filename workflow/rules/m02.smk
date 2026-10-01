# Module 2: contiguity and anchoring (SPEC §8.2).


rule m02_quast:
    input:
        W / "prep" / "asm.fa",
    output:
        W / "m02" / "quast" / "report.tsv",
    params:
        pre=env("quast", "m02_quast"),
    threads: workflow.cores
    log:
        log("m02_quast"),
    benchmark:
        bench("m02", "quast")
    shell:
        "{params.pre}quast.py --large --min-contig 0 --split-scaffolds --no-icarus"
        " --no-plots --threads {threads} -o {W}/m02/quast {input} > {log} 2>&1"


rule m02:
    input:
        scan=rules.scan.output,
        quast=W / "m02" / "quast" / "report.tsv",
        agp=[W / "prep" / "asm.agp"] if config["agp"] else [],
    output:
        W / "m02" / "summary.json",
    params:
        pre=env("core", "m02"),
        agp=f"--agp {W}/prep/asm.agp" if config["agp"] else "",
        out=R / "m02_contiguity",
    log:
        log("m02"),
    benchmark:
        bench("m02", "m02")
    shell:
        "{params.pre}{PY} -m asmqc.m02_contiguity --scan {W}/scan"
        " --quast-report {input.quast} {params.agp} --outdir {params.out}"
        " --work {W}/m02 > {log} 2>&1"
