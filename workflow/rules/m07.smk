# Module 7: redundancy of unplaced scaffolds (SPEC §8.7).


rule m07_align:
    input:
        fa=W / "prep" / "asm.fa",
        fai=W / "prep" / "asm.fa.fai",
    output:
        paf=W / "m07" / "redundancy.paf",
        lst=W / "m07" / "unplaced_ge1kb.txt",
    params:
        pre=env("core", "m07_align"),
    threads: workflow.cores
    log:
        log("m07_align"),
    benchmark:
        bench("m07", "align")
    shell:
        "{params.pre}( {PY} -m asmqc.m07_redundancy lists --fai {input.fai} --work {W}/m07"
        " && touch {output.paf}"
        " && if [ -s {output.lst} ]; then"
        " samtools faidx {input.fa} -r {W}/m07/chroms.txt > {W}/m07/chroms.fa"
        " && samtools faidx {input.fa} -r {output.lst} > {W}/m07/unplaced_ge1kb.fa"
        " && mm2plus -x asm5 -c --secondary=yes -N 5 -t {threads}"
        " {W}/m07/chroms.fa {W}/m07/unplaced_ge1kb.fa > {output.paf}; fi"
        " ) > {log} 2>&1"


rule m07:
    input:
        paf=rules.m07_align.output.paf,
        fai=W / "prep" / "asm.fa.fai",
    output:
        W / "m07" / "summary.json",
    params:
        pre=env("core", "m07"),
        out=R / "m07_redundancy",
    log:
        log("m07"),
    benchmark:
        bench("m07", "m07")
    shell:
        "{params.pre}{PY} -m asmqc.m07_redundancy classify --fai {input.fai}"
        " --paf {input.paf} --outdir {params.out} --work {W}/m07 > {log} 2>&1"
