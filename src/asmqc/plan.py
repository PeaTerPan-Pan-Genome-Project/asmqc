"""Module switching and read-type selection (SPEC §3, §4.3).

The plan is decided once, before Snakemake runs, and written into the
Snakemake config. Modules may still report `skipped_no_input` themselves when
the condition depends on the assembly content (M7: no unplaced >= 1 kb).
"""

from dataclasses import asdict, dataclass

from asmqc.schema import MODULES
from asmqc.validate import RunOptions

ASSEMBLY_ONLY = ("m01", "m02", "m04", "m05", "m06", "m07")


@dataclass
class ModulePlan:
    module: str
    status: str  # run | skipped_no_input | skipped_by_user
    read_type: str | None = None  # illumina | hifi | ont_r9 | ont_r10
    short_reads: bool = False  # M9: Illumina passed as -ngs
    reason: str | None = None


def normalise_modules(selection: list[str] | None) -> list[str] | None:
    """'1,4,11' style selection -> ['m01', 'm04', 'm11']; raises ValueError."""
    if selection is None:
        return None
    out = []
    for item in selection:
        m = f"m{int(item.lstrip('m')):02d}" if item.lstrip("m").isdigit() else item
        if m not in MODULES:
            raise ValueError(f"--modules: unknown module {item!r}")
        out.append(m)
    return out


def make_plan(opts: RunOptions) -> dict[str, ModulePlan]:
    selected = set(normalise_modules(opts.modules) or MODULES)
    ont_type = f"ont_{opts.ont_chemistry}" if opts.ont else None
    short = "illumina" if opts.illumina else "hifi" if opts.hifi else None
    long = "hifi" if opts.hifi else ont_type

    plan: dict[str, ModulePlan] = {}
    for m in MODULES:
        if m not in selected:
            plan[m] = ModulePlan(m, "skipped_by_user", reason="not in --modules")
        elif m in ASSEMBLY_ONLY:
            plan[m] = ModulePlan(m, "run")
        elif m in ("m08", "m11"):
            plan[m] = (ModulePlan(m, "run", read_type=short) if short else
                       ModulePlan(m, "skipped_no_input", reason="no Illumina or HiFi reads"))
        elif m == "m09":
            plan[m] = (ModulePlan(m, "run", read_type=long, short_reads=bool(opts.illumina))
                       if long else
                       ModulePlan(m, "skipped_no_input", reason="no HiFi or ONT reads"))
    return plan


def read_types_needed(plan: dict[str, ModulePlan]) -> list[str]:
    """Read types to map, in a fixed order (shared map_<readtype> stage)."""
    needed = {p.read_type for p in plan.values() if p.status == "run" and p.read_type}
    if plan["m09"].status == "run" and plan["m09"].short_reads:
        needed.add("illumina")
    return [t for t in ("illumina", "hifi", "ont_r9", "ont_r10") if t in needed]


def describe(plan: dict[str, ModulePlan]) -> str:
    lines = []
    for p in plan.values():
        extra = []
        if p.read_type:
            extra.append(f"reads: {p.read_type}")
        if p.short_reads:
            extra.append("+ illumina (-ngs)")
        if p.reason:
            extra.append(p.reason)
        lines.append(f"  {p.module}  {p.status:<17} {'; '.join(extra)}".rstrip())
    return "\n".join(lines)


def as_dict(plan: dict[str, ModulePlan]) -> dict[str, dict]:
    return {m: asdict(p) for m, p in plan.items()}


def budget(cores: int, mem_mb: int, craq: bool) -> dict[str, int]:
    """Thread and memory shares for CRAQ (M9) and the other rules.

    CRAQ runs for hours but uses -t only for samtools view, so it gets a
    small fixed thread share and the other rules the rest; they then run
    while CRAQ runs instead of after it. Memory is split in half because
    CRAQ's depth scripts and meryl are the two large consumers. Without M9
    everything gets the whole machine. Read mapping always uses all cores
    (CRAQ waits for its BAMs anyway).
    """
    if not craq:
        return {"craq_threads": 0, "side_threads": cores, "craq_mem_mb": 0,
                "side_mem_mb": mem_mb}
    craq_threads = max(1, min(8, cores // 2))
    return {"craq_threads": craq_threads, "side_threads": max(1, cores - craq_threads),
            "craq_mem_mb": mem_mb // 2, "side_mem_mb": mem_mb - mem_mb // 2}
