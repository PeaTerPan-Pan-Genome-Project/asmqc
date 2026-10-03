"""Plain-language documentation of every reported value, for the HTML reports.

COLUMNS maps each qc_summary.tsv column to (label, unit, description). The
description is shown as a tooltip and must say what the value is and how it
is computed. tests/unit/test_docs.py fails if a column is missing.

Units drive the display format: "bp" adds a readable size, "%" a percent
sign, "count" thousands separators; "" leaves the value as written.
"""

COLUMNS: dict[str, tuple[str, str, str]] = {
    # identity and run
    "label": ("Assembly label", "", "The --label given to asmqc; names every output."),
    "asmqc_version": ("asmqc version", "",
                      "Image version. Results are comparable only within one MAJOR.MINOR "
                      "version."),
    "assembly_md5": ("Assembly md5", "",
                     "MD5 checksum of the decompressed assembly FASTA as given, before any "
                     "renaming. Identifies exactly which file was checked."),
    "run_date": ("Run date", "", "Date the run started (UTC)."),
    "ena_rules": ("ENA submission rules", "",
                  "PASS if no sequence breaks an ENA rule (length ≥ 20 bp, no leading or "
                  "trailing N, unique names, IUPAC characters only); FAIL otherwise. "
                  "Only these rules decide it."),
    "n_flags_ena_blocking": ("ENA-blocking flags", "count",
                             "Number of problems that ENA would reject. Any value above 0 "
                             "makes the ENA result FAIL."),
    "n_flags_warning": ("Warnings", "count",
                        "Number of findings worth checking that do not block submission."),
    "n_flags_info": ("Information flags", "count",
                     "Number of observations reported for completeness (e.g. IUPAC codes, "
                     "interstitial telomeric repeats)."),
    **{f"{m}_status": (f"{m.upper()} status", "",
                       "ok: finished; skipped_no_input: the inputs it needs were not given; "
                       "skipped_by_user: excluded with --modules; failed: see logs/.")
       for m in ("m01", "m02", "m04", "m05", "m06", "m07", "m08", "m09", "m11")},
    # M1
    "m01_n_seq": ("Number of sequences", "count",
                  "All FASTA records: chr1–chr7 plus unplaced scaffolds."),
    "m01_total_bp": ("Assembly size", "bp", "Sum of all sequence lengths, N included."),
    "m01_n_lt20bp": ("Sequences < 20 bp", "count",
                     "ENA rejects sequences shorter than 20 bp."),
    "m01_n_lt200bp": ("Sequences < 200 bp", "count",
                      "Very short sequences; allowed by ENA but rarely useful."),
    "m01_n_terminal_n": ("Sequences with terminal N", "count",
                         "Sequences that start or end with N. ENA rejects them; trim the Ns."),
    "m01_n_duplicate_names": ("Duplicate names", "count",
                              "Names used by more than one sequence. ENA rejects them."),
    "m01_n_invalid_chars": ("Invalid characters", "count",
                            "Characters outside the IUPAC nucleotide alphabet. ENA rejects "
                            "them."),
    "m01_n_iupac": ("IUPAC ambiguity codes", "count",
                    "Bases other than A, C, G, T and N (e.g. R, Y). Allowed, but often a sign "
                    "of unresolved consensus."),
    "m01_n_gaps_ge10": ("Gaps (N-runs ≥ 10 bp)", "count",
                        "Runs of at least 10 N. Each is a gap between two contigs."),
    "m01_gap_bp_ge10": ("Bases in gaps", "bp", "Total length of the N-runs of ≥ 10 bp."),
    "m01_n_gaps_ge100": ("Gaps ≥ 100 bp", "count",
                         "N-runs of at least 100 bp; scaffolders usually insert 100 N."),
    "m01_n_seq_gt50pct_n": ("Sequences > 50 % N", "count",
                            "Sequences that are mostly gap."),
    "m01_softmask_pct": ("Soft-masked bases", "%",
                         "Share of bases in lower case. Lower case usually marks repeats; 0 "
                         "means the FASTA is not soft-masked."),
    "m01_agp_consistent": ("AGP matches FASTA", "",
                           "yes: every FASTA sequence is built exactly as the AGP says. no: "
                           "names, lengths or coordinates differ (see flags). NA: no AGP."),
    "m01_n_round_lengths": ("Unplaced sequences with round lengths", "count",
                            "Unplaced sequences whose length is a multiple of 1,000 bp. Many "
                            "more than expected suggests breaks placed on a fixed grid during "
                            "manual curation."),
    "m01_round_lengths_expected": ("Round lengths expected by chance", "",
                                   "Number of unplaced sequences ÷ 1,000: the count expected "
                                   "if lengths were random."),
    "m01_n_round_agp_cuts": ("Round AGP cut coordinates", "count",
                             "AGP cut coordinates that are multiples of 1,000. By chance about "
                             "0.1 % are; many more point to breaks placed on a grid."),
    "m01_n_agp_cuts": ("AGP cut coordinates", "count",
                       "All component start and end coordinates examined for roundness."),
    # M2
    "m02_scaffold_n50": ("Scaffold N50", "bp",
                         "Half of the assembly is in sequences at least this long. For a "
                         "chromosome-level assembly it is about one chromosome."),
    "m02_scaffold_l50": ("Scaffold L50", "count",
                         "Number of sequences that together make half of the assembly."),
    "m02_scaffold_n90": ("Scaffold N90", "bp",
                         "90 % of the assembly is in sequences at least this long."),
    "m02_longest_bp": ("Longest sequence", "bp", "Length of the longest sequence."),
    "m02_contig_method": ("Contig definition", "",
                          "agp: contigs are the AGP components. nsplit10: no consistent AGP, "
                          "contigs are the pieces between N-runs of ≥ 10 bp."),
    "m02_contig_n50": ("Contig N50", "bp",
                       "Half of the assembly is in contigs at least this long. The main "
                       "measure of continuity for chromosome-level assemblies."),
    "m02_contig_l50": ("Contig L50", "count",
                       "Number of contigs that together make half of the assembly."),
    "m02_contig_n90": ("Contig N90", "bp", "90 % of the assembly is in contigs at least this "
                                           "long."),
    "m02_n_contigs": ("Number of contigs", "count", "Contigs under the contig definition."),
    "m02_contig_n50_nsplit10": ("Contig N50 (split at N-runs)", "bp",
                                "Contig N50 with contigs split at every N-run of ≥ 10 bp, "
                                "whatever the AGP says. Lower than the AGP value when AGP "
                                "components themselves contain N-runs."),
    "m02_chrom_bp": ("Bases in chr1–chr7", "bp", "Total length of the seven chromosomes."),
    "m02_anchored_pct": ("Anchored to chromosomes", "%",
                         "Share of the assembly in chr1–chr7."),
    "m02_n_unplaced": ("Unplaced sequences", "count",
                       "Sequences not assigned to chr1–chr7."),
    "m02_unplaced_bp": ("Bases in unplaced sequences", "bp",
                        "Total length of the unplaced sequences."),
    "m02_n_gaps": ("Gaps", "count",
                   "AGP gap lines (contig method agp) or N-runs of ≥ 10 bp (nsplit10)."),
    "m02_gap_bp": ("Bases in gaps", "bp", "Total length of those gaps."),
    # M4
    "m04_capped_arms": ("Capped chromosome arms (of 14)", "count",
                        "Arms with a telomeric repeat array (TTTAGGG) within 50 kb of the end, "
                        "on the expected strand."),
    "m04_t2t_chromosomes": ("Telomere-to-telomere chromosomes (of 7)", "count",
                            "Chromosomes with both arms capped and no gap (no N-run of "
                            "≥ 10 bp)."),
    "m04_wrong_orientation_arms": ("Arms with wrong-strand telomere", "count",
                                   "Arms whose terminal telomeric array runs on the "
                                   "unexpected strand: often an inverted end piece."),
    "m04_interstitial_arrays": ("Interstitial telomeric arrays", "count",
                                "Telomeric repeat arrays more than 50 kb from both ends. Can "
                                "be genuine in pea, or mark a misjoin."),
    "m04_unplaced_with_telomere": ("Unplaced sequences with a telomere", "count",
                                   "Unplaced sequences with a telomeric array within 50 kb of "
                                   "an end: chromosome ends that exist but were not anchored."),
    # M5
    "m05_plastid_scaffolds_n": ("Plastid scaffolds", "count",
                                "Unplaced sequences ≥ 80 % covered by the pea plastid genome "
                                "at ≥ 95 % identity."),
    "m05_plastid_scaffolds_bp": ("Bases in plastid scaffolds", "bp", ""),
    "m05_mito_scaffolds_n": ("Mitochondrial scaffolds", "count",
                             "Unplaced sequences ≥ 80 % covered by the pea mitochondrial "
                             "genome at ≥ 95 % identity."),
    "m05_mito_scaffolds_bp": ("Bases in mitochondrial scaffolds", "bp", ""),
    "m05_chrom_plastid_like_bp": ("Plastid-like bases in chromosomes", "bp",
                                  "Chromosome bases aligned to the plastid genome. Nuclear "
                                  "insertions of plastid DNA are normal biology."),
    "m05_chrom_mito_like_bp": ("Mitochondria-like bases in chromosomes", "bp",
                               "Chromosome bases aligned to the mitochondrial genome."),
    "m05_rdna45s_loci": ("45S rDNA array locations", "",
                         "Chromosomes carrying a 45S (18S-5.8S-25S) array of ≥ 3 copies; "
                         "'unplaced' if arrays lie on unplaced sequences. The large 45S arrays "
                         "are the nucleolus organiser regions (NORs)."),
    "m05_rdna45s_copies": ("45S rDNA copies in arrays", "count",
                           "18S copies in the 45S arrays. Measures how much of the arrays "
                           "was assembled, not the copy number in the plant (thousands)."),
    "m05_rdna45s_fragments": ("45S fragments", "count",
                              "Clusters with fewer than 3 copies: dispersed rDNA pieces."),
    "m05_rdna5s_loci": ("5S rDNA array locations", "",
                        "Chromosomes carrying a 5S array of ≥ 10 copies."),
    "m05_rdna5s_copies": ("5S rDNA copies in arrays", "count", "5S copies in the 5S arrays."),
    "m05_rdna5s_fragments": ("5S fragments", "count",
                             "Clusters with fewer than 10 copies: mostly dispersed 5S-like "
                             "sequence."),
    "m05_rdna_only_scaffolds_n": ("rDNA-only scaffolds", "count",
                                  "Unplaced sequences ≥ 80 % covered by one rDNA array: "
                                  "pieces of arrays that could not be anchored."),
    "m05_rdna_only_scaffolds_bp": ("Bases in rDNA-only scaffolds", "bp", ""),
    # M6
    "m06_complete_pct": ("BUSCO complete", "%",
                         "Share of the 7,702 conserved single-copy Fabales genes found complete "
                         "(single + duplicated)."),
    "m06_single_pct": ("BUSCO complete, single copy", "%",
                       "Found complete exactly once: the expected state."),
    "m06_duplicated_pct": ("BUSCO complete, duplicated", "%",
                           "Found complete more than once. In an inbred line mostly assembly "
                           "duplication; some genes are genuinely duplicated."),
    "m06_fragmented_pct": ("BUSCO fragmented", "%", "Found only partially."),
    "m06_missing_pct": ("BUSCO missing", "%", "Not found."),
    "m06_n_markers": ("BUSCO genes searched", "count",
                      "Number of genes in fabales_odb12.2 (always 7,702)."),
    "m06_internal_stop_pct": ("Complete BUSCOs with an internal stop codon", "%",
                              "Complete genes whose predicted protein contains a stop codon: "
                              "usually a frameshifting indel in the assembly (see M11)."),
    "m06_complete_on_unplaced": ("Complete BUSCOs on unplaced sequences", "count",
                                 "Complete genes with a copy on an unplaced sequence."),
    "m06_dup_both_on_chrom": ("Duplicated BUSCOs, all copies on chromosomes", "count",
                              "Duplicated genes whose copies all lie on chr1–chr7."),
    "m06_dup_any_on_unplaced": ("Duplicated BUSCOs with an unplaced copy", "count",
                                "Duplicated genes with a copy on an unplaced sequence: likely "
                                "false duplication by redundant scaffolds (see M7)."),
    "m06_lineage": ("BUSCO lineage", "", "Gene set and its date."),
    # M7
    "m07_n_duplicate": ("Duplicate scaffolds", "count",
                        "Unplaced scaffolds (≥ 1 kb) whose ≥ 90 % aligns in one piece to a "
                        "chromosome at ≥ 99 % identity: redundant copies of chromosome "
                        "sequence."),
    "m07_duplicate_bp": ("Bases in duplicate scaffolds", "bp", ""),
    "m07_n_partial_overlap": ("Partially overlapping scaffolds", "count",
                              "As duplicate, but the single alignment covers 50–90 %."),
    "m07_partial_overlap_bp": ("Bases in partially overlapping scaffolds", "bp", ""),
    "m07_n_repeat_like": ("Repeat-like scaffolds", "count",
                          "No single alignment covers ≥ 50 %, but many short alignments do: "
                          "scaffolds made of repeats present on the chromosomes."),
    "m07_repeat_like_bp": ("Bases in repeat-like scaffolds", "bp", ""),
    "m07_n_unique": ("Unique scaffolds", "count",
                     "Unplaced scaffolds with no substantial match on the chromosomes."),
    "m07_unique_bp": ("Bases in unique scaffolds", "bp", ""),
    "m07_n_short": ("Short scaffolds (< 1 kb)", "count", "Too short to classify."),
    "m07_short_bp": ("Bases in short scaffolds", "bp", ""),
    "m07_total_minus_duplicate_bp": ("Assembly size without duplicates", "bp",
                                     "Assembly size minus the duplicate scaffolds."),
    # M8
    "m08_read_type": ("Reads used for QV", "",
                      "illumina if given, else hifi. ONT reads are never used here."),
    "m08_reads_independent": ("Reads independent of the assembly", "",
                              "yes: the reads were not used to build the assembly, so the QV "
                              "measures accuracy. no: a self-consistency QV, higher than the "
                              "true accuracy."),
    "m08_k": ("k-mer size", "", "Length of the k-mers compared (21, fixed)."),
    "m08_kmer_coverage": ("k-mer coverage", "×",
                          "How many times a single-copy 21-mer occurs in the reads: the main "
                          "peak of the read k-mer histogram. Below 20× the QV is less reliable."),
    "m08_low_coverage": ("Low k-mer coverage", "", "yes if the k-mer coverage is below 20×."),
    "m08_qv": ("Consensus quality (QV)", "QV",
               "QV = −10·log10(error rate), from assembly k-mers missing in the reads. QV 30 = "
               "1 error per 1 kb, 40 = 1 per 10 kb, 50 = 1 per 100 kb."),
    "m08_error_rate": ("Estimated base error rate", "",
                       "Errors per base, estimated from assembly-only k-mers."),
    "m08_completeness_pct": ("k-mer completeness", "%",
                             "Share of reliable read k-mers found in the assembly. Missing "
                             "k-mers mean missing or collapsed sequence."),
    # M9
    "m09_long_read_type": ("Long reads used", "", "hifi if given, else ONT (r9 or r10)."),
    "m09_short_reads_used": ("Short reads also used", "",
                             "yes if Illumina reads were given to CRAQ as well."),
    "m09_coverage": ("Long-read depth", "×", "Median long-read depth over chr1–chr7."),
    "m09_low_coverage": ("Low long-read depth", "", "yes if the depth is below 20×."),
    "m09_aqi": ("Overall AQI", "",
                "Not reported by CRAQ 1.10; left NA until defined (SPEC §13.8)."),
    "m09_r_aqi": ("Regional AQI (R-AQI)", "",
                  "100·exp(−0.1 × regional errors per Mb). 100 = no small-scale errors; CRAQ "
                  "rates > 90 reference quality, 80–90 high, 60–80 draft, < 60 low."),
    "m09_s_aqi": ("Structural AQI (S-AQI)", "",
                  "100·exp(−0.1 × structural errors per Mb), same scale as R-AQI."),
    "m09_n_cre": ("Regional errors (CRE)", "count",
                  "Clip-based regional errors: places where many reads are clipped or broken, "
                  "pointing to small-scale misassembly (local mis-joins, collapsed or expanded "
                  "repeats)."),
    "m09_n_cse": ("Structural errors (CSE)", "count",
                  "Clip-based structural errors: breakpoints where long reads split to "
                  "distant places, pointing to misjoins, inversions or translocations."),
    "m09_cse_near_agp_junction": ("Structural errors at scaffold joins", "count",
                                  "CSEs within 10 kb of an AGP gap or contig join: suggests a "
                                  "scaffolding error."),
    "m09_cse_inside_contig": ("Structural errors inside contigs", "count",
                              "CSEs farther from any join: suggests an assembler error."),
    # M11
    "m11_read_type": ("Reads used", "", "illumina if given, else hifi."),
    "m11_reads_independent": ("Reads independent of the assembly", "",
                              "As for QV: errors are under-counted when the reads built the "
                              "assembly."),
    "m11_callable_bp": ("Callable region", "bp",
                        "Bases with read depth between half and twice the median: the part "
                        "of the assembly where errors can be called reliably."),
    "m11_hp_errors": ("Homopolymer errors", "count",
                      "Homozygous indels (QUAL ≥ 30) in runs of one base, e.g. AAAAAAA, "
                      "of ≥ 4 bp in the assembly or the reads: the reads consistently show a "
                      "different run length from the assembly."),
    "m11_hp_errors_per_mb": ("Homopolymer errors per Mb", "",
                             "Homopolymer errors per million callable bases."),
    "m11_hp_errors_per_10k_runs": ("Homopolymer errors per 10,000 runs", "",
                                   "Errors per 10,000 homopolymer runs (≥ 4 bp) in the "
                                   "callable region: how often a run is wrong."),
    "m11_dinuc_errors": ("Dinucleotide-repeat errors", "count",
                         "Homozygous indels of whole repeat units in dinucleotide repeats "
                         "(e.g. ATATAT)."),
    "m11_dinuc_errors_per_mb": ("Dinucleotide-repeat errors per Mb", "", ""),
    "m11_other_indels_per_mb": ("Other indels per Mb", "",
                                "Homozygous indels outside homopolymers and dinucleotide "
                                "repeats."),
    "m11_snv_per_mb": ("Substitutions per Mb", "", "Homozygous single-base differences."),
    "m11_hp_pct_of_errors": ("Homopolymer share of all errors", "%",
                             "Homopolymer errors as a share of all homozygous differences."),
    "m11_hp_ins_del_ratio": ("Insertions per deletion (homopolymers)", "",
                             "Above 1: the reads have longer runs than the assembly, i.e. "
                             "assembled runs are too short, typical of ONT-based assemblies."),
    "m11_hp_at_pct": ("A/T homopolymer errors", "%",
                      "Share of homopolymer errors in A or T runs."),
    "m11_het_calls": ("Heterozygous calls (not errors)", "count",
                      "Heterozygous differences. In an inbred line mostly reads mapped to the "
                      "wrong repeat copy; never counted as errors."),
}

# Explanatory paragraphs per module, shown above its numbers (SPEC §7.4).
MODULE_TEXT: dict[str, list[str]] = {
    "m01": [
        "Checks the FASTA file itself: names, characters, gaps, and the rules ENA applies "
        "at submission.",
        "ENA rejects sequences shorter than 20 bp, leading or trailing N, duplicate names "
        "and characters outside the IUPAC alphabet. Any of these makes the ENA result FAIL; "
        "everything else is a warning or information.",
        "Round numbers: manual Hi-C curation sometimes breaks scaffolds on a fixed grid "
        "(e.g. every 1 kb) instead of at the true junction. That leaves unplaced pieces whose "
        "lengths, and AGP cut coordinates, are multiples of 1,000. By chance about 0.1 % of "
        "coordinates are; a much higher count suggests grid-placed breaks. Round numbers "
        "never affect the ENA result.",
    ],
    "m02": [
        "How continuous the assembly is, and how much of it is in chr1–chr7.",
        "In a chromosome-level assembly the scaffold N50 is about one chromosome and says "
        "little. The contig N50, the length of the gap-free pieces, is the discriminating "
        "number.",
        "Contigs are the AGP components when an AGP is given and matches the FASTA, "
        "otherwise the pieces between N-runs of at least 10 bp. The two contig N50 values "
        "differ when AGP components themselves contain N-runs.",
    ],
    "m04": [
        "Telomeres are arrays of the repeat TTTAGGG at the chromosome ends (CCCTAAA on the "
        "opposite strand at the start). tidk counts the repeat in 10-kb windows; windows "
        "with at least 25 copies form telomeric bands.",
        "An arm is capped when a band lies within 50 kb of the end on the expected strand. "
        "A chromosome is telomere-to-telomere (T2T) when both arms are capped and it has no "
        "gap (no N-run of ≥ 10 bp).",
        "A band on the unexpected strand (wrong orientation) usually means the end piece is "
        "inverted. Interstitial arrays lie more than 50 kb from both ends; pea carries some "
        "genuine ones, but they can also mark a misjoin. Unplaced sequences with a terminal "
        "band are chromosome ends that exist but were not anchored.",
    ],
    "m05": [
        "Report only; nothing here passes or fails.",
        "Organelles: sequences aligned to the pea plastid and mitochondrial genomes at "
        "≥ 95 % identity. An unplaced sequence ≥ 80 % covered is organelle DNA assembled "
        "as a scaffold. Organelle-like blocks inside chromosomes are normal nuclear "
        "insertions; a very large block can indicate a misjoin.",
        "rDNA: copies of the ribosomal RNA genes found with BLAST (a hit must cover ≥ 50 % "
        "of the subunit). Copies within 20 kb form a cluster; a cluster is an array with "
        "≥ 3 copies (45S, counted by 18S) or ≥ 10 copies (5S), otherwise a fragment. The "
        "large 45S arrays form the nucleolus organiser regions (NORs). Assembled copy numbers "
        "measure how much of an array was captured; the plant has thousands of copies.",
    ],
    "m06": [
        "BUSCO searches 7,702 genes that occur once in nearly every Fabales genome "
        "(fabales_odb12.2). Their recovery measures how complete the gene space is.",
        "Complete percentages saturate near 100 in good assemblies; the duplicated share "
        "and the placement of duplicated genes carry more signal. A duplicated gene with a "
        "copy on an unplaced sequence is a likely false duplication (compare M7).",
        "Internal stop codons in complete genes are usually frameshifting indels in the "
        "assembly; M11 counts such errors directly.",
        "Synteny (report only): each point is a gene found once in both this assembly and "
        "Caméor v2, plotted at its two positions. Collinear chromosomes give diagonals; "
        "segments elsewhere are translocations, orange (opposite orientation) segments "
        "inversions. The table gives, per chromosome, the Caméor chromosome holding most of "
        "its genes and the overall orientation.",
    ],
    "m07": [
        "A measurement, never a purge. Every unplaced scaffold of ≥ 1 kb is aligned to "
        "chr1–chr7.",
        "duplicate: one alignment covers ≥ 90 % of it at ≥ 99 % identity (MAPQ ≥ 20), i.e. "
        "a redundant copy of chromosome sequence. partial overlap: the same, covering "
        "50–90 %. repeat-like: no single alignment covers half, but many short ones do. "
        "Coverage is not summed for the duplicate class because pea is about 85 % repeats: "
        "summed coverage would call every repeat scaffold a duplicate.",
    ],
    "m08": [
        "Merqury compares the 21-mers (all 21-bp substrings) of the assembly with those of "
        "accurate reads (Illumina if given, else HiFi). A 21-mer present in the assembly "
        "but never in the reads most likely contains an assembly error.",
        "QV (consensus quality) puts that error rate on the Phred scale: QV = "
        "−10·log10(error rate). QV 30 means one error per 1,000 bases, QV 40 one per "
        "10,000, QV 50 one per 100,000; every 10 points is a tenfold lower error rate.",
        "k-mer completeness is the share of reliable read 21-mers (seen often enough not to "
        "be sequencing errors) that are present in the assembly. Values below about 95 % "
        "point to missing sequence, often collapsed repeats.",
        "k-mer coverage is the position of the main peak of the read 21-mer histogram, "
        "about the read depth. Below 20× real 21-mers are missed more often, which lowers "
        "completeness and makes the QV less reliable.",
        "If the reads were used to build or polish the assembly, errors they share with it "
        "are invisible: the QV then measures self-consistency and is higher than the true "
        "accuracy. QVs are comparable only between assemblies with the same declaration.",
        "Spectra-cn plot: x = how often a 21-mer occurs in the reads, y = how many distinct "
        "21-mers do so; colour = how often it occurs in the assembly. The main peak should be "
        "red (once). Black at the main peak is read sequence missing from the assembly; blue "
        "(twice) at the main peak is false duplication. The black spike at very low counts "
        "is sequencing errors in the reads and is expected.",
    ],
    "m09": [
        "CRAQ maps the reads back to the assembly and looks for places where many reads "
        "are clipped (aligned only in part) or split. Where the assembly is right, reads "
        "align end to end.",
        "CRE, clip-based regional error: a short region where reads are clipped, i.e. a "
        "small-scale error such as a local misjoin, an indel cluster, or a collapsed or "
        "expanded repeat. CSE, clip-based structural error: a breakpoint where long reads "
        "split and continue elsewhere, i.e. a misjoin, inversion or translocation.",
        "Sites where only about half of the reads are clipped are heterozygous and reported "
        "separately (CRH, CSH); they are not errors.",
        "R-AQI and S-AQI summarise CREs and CSEs: AQI = 100·exp(−0.1 × errors per Mb). 100 "
        "means none; CRAQ rates > 90 as reference quality, 80–90 high, 60–80 draft, < 60 "
        "low.",
        "A CSE within 10 kb of an AGP gap or contig join points to a scaffolding (Hi-C) "
        "error; one inside a contig to an assembler error. CRAQ also flags the very ends of "
        "sequences and gaps, where reads are clipped by necessity. Long-read depth below 20× "
        "makes both error types less reliable.",
    ],
    "m11": [
        "Assemblies built from Oxford Nanopore reads often get the length of homopolymer "
        "runs (one base repeated, e.g. AAAAAAAA) wrong; usually the assembled run is too "
        "short. Module 11 finds such errors with accurate reads (Illumina if given, else "
        "HiFi).",
        "Reads are mapped and variants called where the read depth is typical (the callable "
        "region: half to twice the median depth). A homozygous call, where all reads "
        "disagree with the assembly, is an assembly error in an inbred line. Heterozygous "
        "calls mostly come from reads placed on the wrong repeat copy and are not counted.",
        "Each error is classified: homopolymer (the inserted or deleted bases are one base "
        "and the run is ≥ 4 bp), dinucleotide repeat (whole units of a 2-bp repeat), other "
        "indel, or substitution.",
        "Rates are per million callable bases and per 10,000 homopolymer runs. An "
        "insertion/deletion ratio well above 1 means the reads have longer runs than the "
        "assembly: the assembled runs are too short. A high A/T share is typical, since long "
        "runs are mostly A or T.",
        "In genes, an indel whose length is not a multiple of 3 shifts the reading frame "
        "and usually truncates the protein. The table below counts errors inside the coding "
        "exons of the conserved BUSCO genes (M6); such errors appear in M6 as internal stop "
        "codons.",
    ],
}

FLAGS: dict[str, str] = {
    "seq_lt_20bp": "Sequence shorter than 20 bp. ENA rejects it; remove it.",
    "terminal_n": "Sequence starts or ends with N. ENA rejects it; trim the Ns.",
    "duplicate_name": "Two or more sequences share a name. ENA rejects it; rename.",
    "invalid_char": "Character outside the IUPAC nucleotide alphabet (or a bad name).",
    "empty_sequence": "Record with no sequence.",
    "seq_lt_200bp": "Sequence shorter than 200 bp: allowed, rarely useful.",
    "n_fraction_gt_50pct": "More than half of the sequence is N.",
    "agp_mismatch": "The AGP does not describe the FASTA exactly (names, lengths or "
                    "coordinates). Contigs are then taken from N-runs instead.",
    "crlf_line_endings": "Windows line endings (CRLF) in the FASTA.",
    "iupac_present": "IUPAC ambiguity codes in the sequence.",
    "inconsistent_line_width": "FASTA lines of varying width; harmless for most tools.",
    "round_length": "Unplaced sequence whose length is a multiple of 1,000 bp.",
    "round_agp_cut": "AGP cut at a coordinate that is a multiple of 1,000.",
    "organelle_scaffold": "Unplaced sequence that is plastid or mitochondrial DNA. ENA accepts "
                          "organelle sequence only when declared as such.",
    "wrong_orientation_telomere": "Terminal telomeric array on the unexpected strand.",
    "interstitial_telomere": "Telomeric array more than 50 kb from both chromosome ends.",
    "low_coverage": "Read coverage below 20×: results of this module are less reliable.",
    "reads_not_independent": "The reads were used to build the assembly; the module measures "
                             "self-consistency, not accuracy.",
}

GLOSSARY: list[tuple[str, str]] = sorted([
    ("AGP", "File describing how contigs and gaps are joined into scaffolds and "
            "chromosomes."),
    ("AQI", "Assembly Quality Index of CRAQ: 100·exp(−0.1 × errors per Mb); 100 means no "
            "errors found."),
    ("BUSCO", "Benchmarking Universal Single-Copy Orthologs: conserved genes expected once "
              "in every Fabales genome; their recovery measures gene-space completeness."),
    ("Callable region", "Bases with typical read depth, where variants can be called "
                        "reliably."),
    ("Contig", "Continuous sequence without gaps."),
    ("CRE / CSE", "CRAQ regional (small-scale) and structural (breakpoint) errors, from "
                  "clipped read alignments."),
    ("Homopolymer (HP)", "Run of one base, e.g. AAAAAAA. Long runs are hard for ONT and older "
                         "sequencing, so their length is often wrong in assemblies."),
    ("Dinucleotide repeat (STR2)", "Tandem repeat of a 2-bp unit, e.g. ATATATAT."),
    ("Homozygous / heterozygous call", "A difference seen in all reads (homozygous: an error "
                                       "in an inbred line) or in about half of them "
                                       "(heterozygous)."),
    ("k-mer", "Every substring of length k (here 21) of reads or assembly."),
    ("N50 / L50", "N50: the length such that half of the assembly is in pieces at least that "
                  "long. L50: how many pieces that takes."),
    ("NOR", "Nucleolus organiser region: chromosome locus of the 45S rDNA arrays."),
    ("QV", "Phred-scaled consensus quality: QV = −10·log10(error rate)."),
    ("rDNA (45S, 5S)", "Ribosomal RNA genes in long tandem arrays; 45S units contain 18S, "
                       "5.8S and 25S."),
    ("Spectra-cn plot", "Histogram of read k-mer counts coloured by how often each k-mer "
                        "occurs in the assembly. Black (read-only) k-mers at the main peak "
                        "are sequence missing from the assembly; k-mers at the main peak "
                        "found twice are false duplication."),
    ("T2T", "Telomere-to-telomere: both arms capped and no gap."),
    ("Telomere", "Chromosome end made of TTTAGGG repeats in pea; CCCTAAA on the opposite "
                 "strand at the start of a chromosome."),
], key=lambda g: g[0].lower())

FILES: list[tuple[str, str]] = [
    ("qc_summary.tsv", "One row of all values in this report, fixed column order."),
    ("flags.tsv", "Every flag with sequence, position and message."),
    ("run_manifest.json", "Versions, inputs, parameters, reference data, host, run times."),
    ("report.html", "This report."),
    ("m01_integrity/sequences.tsv", "Per-sequence length, composition, N ends, line width."),
    ("m01_integrity/gaps.tsv", "Every N-run of ≥ 10 bp."),
    ("m01_integrity/checksums.md5", "MD5 of the FASTA and AGP."),
    ("m02_contiguity/per_chromosome.tsv", "Per chromosome: length, contigs, gaps, contig N50."),
    ("m02_contiguity/quast_report.tsv", "QUAST statistics."),
    ("m04_telomeres/telomeres.tsv", "Status of every chromosome arm."),
    ("m04_telomeres/interstitial.tsv", "Interstitial telomeric arrays."),
    ("m04_telomeres/tidk_windows.tsv", "Raw telomeric repeat counts per 10-kb window."),
    ("m04_telomeres/karyoplot.png", "Telomere plot."),
    ("m05_organelle_rdna/rdna_arrays.tsv", "Every rDNA cluster: position, family, copies, "
                                           "array or fragment."),
    ("m05_organelle_rdna/organelle_scaffolds.tsv", "Unplaced sequences with organelle "
                                                   "similarity."),
    ("m05_organelle_rdna/organelle_on_chromosomes.tsv", "Organelle-like bases per "
                                                        "chromosome."),
    ("m06_busco/full_table.tsv", "BUSCO result per gene."),
    ("m06_busco/busco_cds.bed.gz", "Coding exons of the complete BUSCO genes."),
    ("m06_busco/synteny.tsv", "Best-matching Caméor v2 chromosome per chromosome."),
    ("m07_redundancy/redundancy.tsv", "Class of every unplaced scaffold."),
    ("m08_merqury/", "Merqury QV, completeness, spectra plots, assembly-only k-mers."),
    ("m09_craq/CRE.bed, CSE.bed", "Positions of CRAQ regional and structural errors."),
    ("m11_homopolymer/errors.bed.gz", "Every homopolymer and dinucleotide-repeat error."),
    ("m11_homopolymer/hom_calls.vcf.gz", "All homozygous calls used."),
    ("logs/", "Tool logs (local paths replaced by placeholders)."),
]


def fmt(column: str, value: str) -> str:
    """Display form of a qc_summary.tsv value."""
    if value in ("NA", "", None):
        return "not available"
    unit = COLUMNS.get(column, ("", "", ""))[1]
    try:
        if unit == "bp":
            n = int(value)
            return f"{n:,} bp" + (f" ({human_bp(n)})" if n >= 10_000 else "")
        if unit == "count":
            return f"{int(value):,}"
        if unit == "%":
            return f"{value} %"
        if unit in ("×", "QV"):
            return f"{value} {unit}" if unit == "×" else value
    except ValueError:
        return value
    return value


def human_bp(n: float) -> str:
    for unit, div in (("Gb", 1e9), ("Mb", 1e6), ("kb", 1e3)):
        if n >= div:
            return f"{n / div:.3g} {unit}"
    return f"{n:.0f} bp"
