"""Fixed `qc_summary.tsv` schema (SPEC §7.2).

The column order here is the contract. Changing it is a minor or major
version bump (SPEC principle E). tests/unit/test_schema.py compares it with
SPEC §7.2.

Value kinds and their formatting:
  str    as is
  enum   one of the listed values
  int    integer
  pct    0-100, 2 decimals
  rate   4 significant digits
  dec2   2 decimals (QV, coverage, AQI: not percentages, counts or rates)
Missing values are written as NA.
"""

from dataclasses import dataclass

NA = "NA"
MODULES = ["m01", "m02", "m04", "m05", "m06", "m07", "m08", "m09", "m11"]
STATUSES = ("ok", "skipped_no_input", "skipped_by_user", "failed")
YES_NO = ("yes", "no")


@dataclass(frozen=True)
class Column:
    name: str
    kind: str
    values: tuple[str, ...] = ()

    @property
    def module(self) -> str | None:
        prefix = self.name[:3]
        return prefix if prefix in MODULES and not self.name.endswith("_status") else None


def _c(kind: str, *names: str, values: tuple[str, ...] = ()) -> list[Column]:
    return [Column(n, kind, values) for n in names]


COLUMNS: list[Column] = [
    # identity and run
    *_c("str", "label", "asmqc_version", "assembly_md5", "run_date"),
    *_c("enum", "ena_rules", values=("PASS", "FAIL")),
    *_c("int", "n_flags_ena_blocking", "n_flags_warning", "n_flags_info"),
    *_c("enum", *(f"{m}_status" for m in MODULES), values=STATUSES),
    # module 1
    *_c("int", "m01_n_seq", "m01_total_bp", "m01_n_lt20bp", "m01_n_lt200bp",
        "m01_n_terminal_n", "m01_n_duplicate_names", "m01_n_invalid_chars",
        "m01_n_iupac", "m01_n_gaps_ge10", "m01_gap_bp_ge10", "m01_n_gaps_ge100",
        "m01_n_seq_gt50pct_n"),
    *_c("pct", "m01_softmask_pct"),
    *_c("enum", "m01_agp_consistent", values=YES_NO),
    *_c("int", "m01_n_round_lengths"),
    *_c("rate", "m01_round_lengths_expected"),
    *_c("int", "m01_n_round_agp_cuts", "m01_n_agp_cuts"),
    # module 2
    *_c("int", "m02_scaffold_n50", "m02_scaffold_l50", "m02_scaffold_n90", "m02_longest_bp"),
    *_c("enum", "m02_contig_method", values=("agp", "nsplit10")),
    *_c("int", "m02_contig_n50", "m02_contig_l50", "m02_contig_n90", "m02_n_contigs",
        "m02_contig_n50_nsplit10", "m02_chrom_bp"),
    *_c("pct", "m02_anchored_pct"),
    *_c("int", "m02_n_unplaced", "m02_unplaced_bp", "m02_n_gaps", "m02_gap_bp"),
    # module 4
    *_c("int", "m04_capped_arms", "m04_t2t_chromosomes", "m04_wrong_orientation_arms",
        "m04_interstitial_arrays", "m04_unplaced_with_telomere"),
    # module 5
    *_c("int", "m05_plastid_scaffolds_n", "m05_plastid_scaffolds_bp",
        "m05_mito_scaffolds_n", "m05_mito_scaffolds_bp",
        "m05_chrom_plastid_like_bp", "m05_chrom_mito_like_bp"),
    *_c("str", "m05_rdna45s_loci"),
    *_c("int", "m05_rdna45s_copies"),
    *_c("str", "m05_rdna5s_loci"),
    *_c("int", "m05_rdna5s_copies", "m05_rdna_only_scaffolds_n", "m05_rdna_only_scaffolds_bp"),
    # module 6
    *_c("pct", "m06_complete_pct", "m06_single_pct", "m06_duplicated_pct",
        "m06_fragmented_pct", "m06_missing_pct"),
    *_c("int", "m06_n_markers"),
    *_c("pct", "m06_internal_stop_pct"),
    *_c("int", "m06_complete_on_unplaced", "m06_dup_both_on_chrom", "m06_dup_any_on_unplaced"),
    *_c("str", "m06_lineage"),
    # module 7
    *_c("int", "m07_n_duplicate", "m07_duplicate_bp", "m07_n_partial_overlap",
        "m07_partial_overlap_bp", "m07_n_repeat_like", "m07_repeat_like_bp",
        "m07_n_unique", "m07_unique_bp", "m07_n_short", "m07_short_bp",
        "m07_total_minus_duplicate_bp"),
    # module 8
    *_c("enum", "m08_read_type", values=("illumina", "hifi")),
    *_c("enum", "m08_reads_independent", values=YES_NO),
    *_c("int", "m08_k"),
    *_c("dec2", "m08_kmer_coverage"),
    *_c("enum", "m08_low_coverage", values=YES_NO),
    *_c("dec2", "m08_qv"),
    *_c("rate", "m08_error_rate"),
    *_c("pct", "m08_completeness_pct"),
    # module 9
    *_c("enum", "m09_long_read_type", values=("hifi", "ont_r9", "ont_r10")),
    *_c("enum", "m09_short_reads_used", values=YES_NO),
    *_c("dec2", "m09_coverage"),
    *_c("enum", "m09_low_coverage", values=YES_NO),
    *_c("dec2", "m09_aqi", "m09_r_aqi", "m09_s_aqi"),
    *_c("int", "m09_n_cre", "m09_n_cse", "m09_cse_near_agp_junction", "m09_cse_inside_contig"),
    # module 11
    *_c("enum", "m11_read_type", values=("illumina", "hifi")),
    *_c("enum", "m11_reads_independent", values=YES_NO),
    *_c("int", "m11_callable_bp", "m11_hp_errors"),
    *_c("rate", "m11_hp_errors_per_mb", "m11_hp_errors_per_10k_runs"),
    *_c("int", "m11_dinuc_errors"),
    *_c("rate", "m11_dinuc_errors_per_mb", "m11_other_indels_per_mb", "m11_snv_per_mb"),
    *_c("pct", "m11_hp_pct_of_errors"),
    *_c("rate", "m11_hp_ins_del_ratio"),
    *_c("pct", "m11_hp_at_pct"),
    *_c("int", "m11_het_calls"),
]

BY_NAME: dict[str, Column] = {c.name: c for c in COLUMNS}
HEADER: list[str] = [c.name for c in COLUMNS]


def module_columns(module: str) -> list[str]:
    """Metric columns a module must report in its summary.json."""
    return [c.name for c in COLUMNS if c.module == module]


def format_value(column: str, value) -> str:
    """Format one value for qc_summary.tsv; raises ValueError on a bad value."""
    if value is None or value == NA:
        return NA
    col = BY_NAME[column]
    match col.kind:
        case "str":
            s = str(value)
            if "\t" in s or "\n" in s:
                raise ValueError(f"{column}: tab or newline in value")
            return s
        case "enum":
            if value not in col.values:
                raise ValueError(f"{column}: {value!r} not in {col.values}")
            return value
        case "int":
            if isinstance(value, bool) or int(value) != value:
                raise ValueError(f"{column}: {value!r} is not an integer")
            return str(int(value))
        case "pct":
            if not 0 <= value <= 100:
                raise ValueError(f"{column}: {value!r} outside 0-100")
            return f"{value:.2f}"
        case "dec2":
            return f"{value:.2f}"
        case "rate":
            return f"{value:.4g}"
    raise AssertionError(col.kind)
