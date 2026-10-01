# Module 8: QV, k-mer completeness, spectra-cn with Merqury (SPEC §8.8).


rule m08_merqury:
    input:
        db=W / "meryl" / "reads.meryl",
        fa=W / "prep" / "asm.fa",
    output:
        W / "m08" / "merqury" / f"{config['label']}.qv",
    params:
        pre=env("merqury", "m08"),
    threads: workflow.cores
    log:
        log("m08_merqury"),
    benchmark:
        bench("m08", "merqury")
    shell:
        "{params.pre}( mkdir -p {W}/m08/merqury && cd {W}/m08/merqury"
        " && export OMP_NUM_THREADS={threads}"
        " && merqury.sh {input.db} {input.fa} {config[label]} ) > {log} 2>&1"


rule m08:
    input:
        rules.m08_merqury.output,
    output:
        W / "m08" / "summary.json",
    params:
        pre=env("core", "m08_summary"),
        out=R / "m08_merqury",
    log:
        log("m08"),
    benchmark:
        bench("m08", "m08")
    shell:
        "{params.pre}{PY} -m asmqc.m08_merqury --merqury-dir {W}/m08/merqury"
        " --prefix {config[label]} --read-type {PLAN[m08][read_type]}"
        " --reads-used-in-assembly {config[reads_used_in_assembly]}"
        " --outdir {params.out} --work {W}/m08 > {log} 2>&1"
