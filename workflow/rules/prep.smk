# Shared stage: decompressed, renamed working copy of the assembly (and AGP).

PREP_OUT = {
    "fa": W / "prep" / "asm.fa",
    "fai": W / "prep" / "asm.fa.fai",
    "md5": W / "prep" / "asm.md5",
}
if config["agp"]:
    PREP_OUT["agp"] = W / "prep" / "asm.agp"


rule prep:
    input:
        fa=config["assembly"],
        map=W / "chromosome_map.json",
    output:
        **PREP_OUT,
    params:
        agp=f"--agp {config['agp']}" if config["agp"] else "",
        pre=env("core", "prep"),
    log:
        log("prep"),
    benchmark:
        bench("prep", "prep")
    shell:
        "{params.pre}"
        "( {PY} -m asmqc.prep --assembly {input.fa} {params.agp} --map {input.map}"
        " --outdir {W}/prep && samtools faidx {output.fa} ) > {log} 2>&1"
