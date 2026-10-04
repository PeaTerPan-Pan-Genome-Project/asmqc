# Module 6: gene-space completeness, BUSCO (SPEC §8.6).

BUSCO_DL = Path(config["refs"]) / "busco_downloads"
LINEAGE = "fabales_odb12.2"


rule m06_busco:
    input:
        W / "prep" / "asm.fa",
    output:
        W / "m06" / "busco" / f"short_summary.specific.{LINEAGE}.busco.json",
        W / "m06" / "busco" / f"run_{LINEAGE}" / "full_table.tsv",
    params:
        pre=env("busco", "m06_busco"),
        check=env("core", "m06_check"),
    threads: SIDE_THREADS
    log:
        log("m06_busco"),
    benchmark:
        bench("m06", "busco")
    shell:
        # the lineage is asserted before and after the run (SPEC §8.6)
        "( {params.check}{PY} -m asmqc.m06_busco check-lineage"
        " --lineage-dir {BUSCO_DL}/lineages/{LINEAGE}"
        " && {params.pre}cd {W}/m06 && busco --in {input} --mode genome"
        " --lineage_dataset {LINEAGE} --offline --opt-out-run-stats"
        " --download_path {BUSCO_DL} --cpu {threads} --out busco --out_path {W}/m06 -f"
        " ) > {log} 2>&1"


rule m06:
    input:
        busco=rules.m06_busco.output,
        fai=W / "prep" / "asm.fa.fai",
    output:
        W / "m06" / "summary.json",
    params:
        pre=env("core", "m06"),
        out=R / "m06_busco",
    log:
        log("m06"),
    benchmark:
        bench("m06", "m06")
    shell:
        "{params.pre}{PY} -m asmqc.m06_busco summarise --busco-dir {W}/m06/busco"
        " --fai {input.fai} --outdir {params.out} --work {W}/m06 > {log} 2>&1"
