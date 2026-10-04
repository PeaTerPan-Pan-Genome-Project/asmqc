# Shared read stages (SPEC §5.2): one mm2-plus mapping per read type, and the
# read k-mer database for M8.


def read_units(rt):
    """Mapping units: R1,R2 pairs for Illumina, single files otherwise."""
    if rt == "illumina":
        il = config["reads"]["illumina"]
        return [f"{a},{b}" for a, b in zip(il[0::2], il[1::2])]
    return config["reads"]["hifi" if rt == "hifi" else "ont"]


rule map_reads:
    input:
        fa=W / "prep" / "asm.fa",
        fai=W / "prep" / "asm.fa.fai",
        reads=lambda w: [f for u in read_units(w.rt) for f in u.split(",")],
    output:
        bam=W / "map" / "{rt}.bam",
        csi=W / "map" / "{rt}.bam.csi",  # CSI: BAI cannot hold positions > 2^29
    wildcard_constraints:
        rt="illumina|hifi|ont_r9|ont_r10",
    params:
        pre=lambda w: env("core", f"map_{w.rt}"),
        units=lambda w: " ".join(f"'{u}'" for u in read_units(w.rt)),
    threads: workflow.cores
    log:
        log("map_{rt}"),
    benchmark:
        bench("map", "{rt}")
    shell:
        "{params.pre}{PY} -m asmqc.mapping --fa {input.fa} --read-type {wildcards.rt}"
        " --label {config[label]} --threads {threads} --mem-mb {config[mem_mb]}"
        " --out {output.bam} {params.units} > {log} 2>&1"


M08_RT = PLAN["m08"]["read_type"]
M08_FILES = (config["reads"]["illumina"] if M08_RT == "illumina"
             else config["reads"]["hifi"] if M08_RT == "hifi" else [])


def meryl_cmd(wildcards, threads):
    """meryl count per read file (k = 21, fixed), then union-sum."""
    mem = max(1, BUDGET["side_mem_mb"] // 1024 - 2)
    parts = [f"{W}/meryl/part{i}.meryl" for i in range(len(M08_FILES))]
    count = [f"meryl count k=21 memory={mem} threads={threads} output {p} {f}"
             for p, f in zip(parts, M08_FILES)]
    return " && ".join([f"mkdir -p {W}/meryl", *count,
                        f"meryl union-sum output {W}/meryl/reads.meryl {' '.join(parts)}",
                        f"rm -rf {' '.join(parts)}"])


rule meryl_reads:
    input:
        M08_FILES,
    output:
        directory(W / "meryl" / "reads.meryl"),
    params:
        pre=env("merqury", "meryl"),
        cmd=meryl_cmd,
    threads: SIDE_THREADS
    resources:
        mem_mb=BUDGET["side_mem_mb"],
    log:
        log("meryl_reads"),
    benchmark:
        bench("meryl", "reads")
    shell:
        "{params.pre}( {params.cmd} ) > {log} 2>&1"
