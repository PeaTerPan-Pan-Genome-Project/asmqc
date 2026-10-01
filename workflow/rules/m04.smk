# Module 4: telomeres (SPEC §8.4).


rule m04_tidk:
    input:
        W / "prep" / "asm.fa",
    output:
        W / "m04" / "tidk" / "asm_telomeric_repeat_windows.tsv",
    params:
        pre=env("core", "m04_tidk"),
    log:
        log("m04_tidk"),
    benchmark:
        bench("m04", "tidk")
    shell:
        "{params.pre}tidk search --string TTTAGGG --window 10000 --output asm"
        " --dir {W}/m04/tidk {input} > {log} 2>&1"


rule m04:
    input:
        windows=rules.m04_tidk.output,
        fai=W / "prep" / "asm.fa.fai",
    output:
        W / "m04" / "summary.json",
    params:
        pre=env("core", "m04"),
        out=R / "m04_telomeres",
    log:
        log("m04"),
    benchmark:
        bench("m04", "m04")
    shell:
        "{params.pre}{PY} -m asmqc.m04_telomeres --windows {input.windows} --fai {input.fai}"
        " --outdir {params.out} --work {W}/m04 > {log} 2>&1"
