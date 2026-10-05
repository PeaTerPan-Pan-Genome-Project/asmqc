# Module 5: organelle and rDNA inventory (SPEC §8.5).

REFS = Path(config["refs"])


rule m05_organelle_ref:
    input:
        REFS / "plastid_NC_014057.1.fa",
        REFS / "mito_PP555264.1.fa",
    output:
        W / "m05" / "organelles.fa",
    shell:
        "cat {input} > {output}"


rule m05_organelles:
    input:
        ref=rules.m05_organelle_ref.output,
        fa=W / "prep" / "asm.fa",
    output:
        W / "m05" / "organelles.paf",
    params:
        pre=env("core", "m05_organelles"),
    threads: SIDE_THREADS
    log:
        log("m05_organelles"),
    benchmark:
        bench("m05", "organelles")
    shell:
        "{params.pre}mm2plus -x asm20 -c -t {threads} {input.ref} {input.fa}"
        " > {output} 2> {log}"


rule m05_rdna:
    input:
        lib=REFS / "rdna_library.fasta",
        fa=W / "prep" / "asm.fa",
    output:
        W / "m05" / "rdna_blast.tsv",
    params:
        pre=env("core", "m05_rdna"),
        db=W / "m05" / "blastdb" / "asm",
    threads: SIDE_THREADS
    log:
        log("m05_rdna"),
    benchmark:
        bench("m05", "rdna")
    shell:
        "{params.pre}( makeblastdb -dbtype nucl -in {input.fa} -out {params.db}"
        " && blastn -task blastn -query {input.lib} -db {params.db} -evalue 1e-10"
        " -outfmt '6 std qlen slen' -max_target_seqs 1000000 -max_hsps 1000000"
        " -num_threads {threads} -out {output} ) > {log} 2>&1"


rule m05:
    input:
        paf=rules.m05_organelles.output,
        blast=rules.m05_rdna.output,
        fai=W / "prep" / "asm.fa.fai",
    output:
        W / "m05" / "summary.json",
    params:
        pre=env("core", "m05"),
        out=R / "m05_organelle_rdna",
    log:
        log("m05"),
    benchmark:
        bench("m05", "m05")
    shell:
        "{params.pre}{PY} -m asmqc.m05_organelle_rdna --paf {input.paf}"
        " --blast {input.blast} --fai {input.fai} --outdir {params.out}"
        " --work {W}/m05 > {log} 2>&1"
