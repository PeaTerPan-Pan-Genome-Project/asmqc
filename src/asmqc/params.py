"""Every fixed parameter named in SPEC §8, in one place.

Modules read their values from here; the manifest records the whole table
under `parameters` (SPEC §7.3). Changing a value is a minor or major version
bump (SPEC principle E).
"""

PARAMS: dict[str, dict] = {
    "m01": {
        "min_len_ena": 20,
        "min_len_warning": 200,
        "n_fraction_warning": 0.5,
        "gap_min_bp": 10,
        "gap_long_bp": 100,
        "round_unit_bp": 1000,
    },
    "m02": {"nsplit_min_gap": 10},
    "m04": {
        "motif": "TTTAGGG",
        "window": 10000,
        "min_repeats": 25,
        "terminal_bp": 50000,
    },
    "m05": {
        "organelle_min_identity": 0.95,
        "organelle_scaffold_min_cov": 0.80,
        "rdna_evalue": 1e-10,
        "rdna_array_link_bp": 20000,
        "rdna_only_min_frac": 0.80,
    },
    "m06": {
        "lineage": "fabales_odb12.2",
        "lineage_date": "2026-05-13",
        "n_markers": 7702,
    },
    "m07": {
        "min_len_bp": 1000,
        "min_identity": 0.99,
        "min_mapq": 20,
        "duplicate_min_cov": 0.90,
        "partial_min_cov": 0.50,
        "repeat_min_cov": 0.50,
        "max_secondary": 5,
    },
    "m08": {"k": 21, "low_coverage": 20},
    "m09": {"low_coverage": 20, "junction_window_bp": 10000},
    "m11": {
        "min_mapq": 20,
        "min_baseq": 20,
        "max_depth": 500,
        "callable_depth_low": 0.5,
        "callable_depth_high": 2.0,
        "depth_sample_bp": 1000,
        "min_qual": 30,
        "hp_min_len": 4,
        "str2_min_copies": 2,
    },
}
