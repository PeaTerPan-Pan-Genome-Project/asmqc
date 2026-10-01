# Module 1: integrity, format and ENA rules (SPEC §8.1).

M01_AGP = (
    f"--agp {W}/prep/asm.agp --agp-original {config['agp']}" if config["agp"] else ""
)


rule m01:
    input:
        scan=rules.scan.output,
        md5=W / "prep" / "asm.md5",
        agp=[W / "prep" / "asm.agp"] if config["agp"] else [],
    output:
        W / "m01" / "summary.json",
    params:
        pre=env("core", "m01"),
        agp=M01_AGP,
        name=Path(config["assembly"]).name,
        out=R / "m01_integrity",
    log:
        log("m01"),
    benchmark:
        bench("m01", "m01")
    shell:
        "{params.pre}{PY} -m asmqc.m01_integrity --scan {W}/scan --md5 {input.md5}"
        " --assembly-name {params.name:q} {params.agp} --outdir {params.out}"
        " --work {W}/m01 > {log} 2>&1"
