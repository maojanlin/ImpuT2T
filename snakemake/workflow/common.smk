# Shared configuration and path helpers (included by master and chromosome workflows).

import sys

OUT = config["output_dir"]

SMALL_PAN = config["small_pan"]
QUERIES = config["queries"]
QUERY_NAMES = [item[0] for item in QUERIES]
SMALL_PAN_NAMES = [item[0] for item in SMALL_PAN]
THREADS = config["threads"]
RAGTAG_THREADS = config["ragtag_threads"]
SAMPLE = config["sample"]
SMALL_PAN_PREFIXES = config["small_pan_prefixes"]
TARGET_CHROMOSOMES = config["target_chromosomes"]
QUERY_TO_CHROMOSOMES = {QUERY_NAMES[i]: TARGET_CHROMOSOMES[i] for i in range(len(QUERY_NAMES))}

SMALL_PAN_PREFIXES_DICT = {item[0]: item[1] for item in SMALL_PAN_PREFIXES}
SMALL_PAN_DICT = {item[0]: item[1] for item in SMALL_PAN}
QUERIES_DICT = {item[0]: item[1] for item in QUERIES}

PANGENOME_NAMES = config["pangenome_names"]
PANGENOME_PATH = config["pangenome_path"]

FLAT_THRESHOLD = config.get("heavy_jobs_flat_threshold", 6000)

# Child chromosome workflow passes scope via --config; master sets globals in snakefile.
CHILD_MODE = (
    config.get("run_tag") is not None
    and config.get("query_name") is not None
    and config.get("chromosome") is not None
)

if CHILD_MODE:
    RUN_TAG = config["run_tag"]
    SAMPLE_LIST_FILE = config["sample_list_file"]
    QUERY_CHROMOSOME_PAIRS = [(config["query_name"], config["chromosome"])]
else:
    QUERY_CHROMOSOME_PAIRS = [
        (q, c) for q in QUERY_NAMES for c in QUERY_TO_CHROMOSOMES[q]
    ]


def patch_filtered_fasta(sample, query_name, chromosome, run_tag):
    return (
        f"{OUT}/impuT2T_patch/{run_tag}/"
        f"{sample}.{query_name}.{chromosome}.aggregate.patch.filtered.fasta"
    )


def patch_outputs():
    return [
        patch_filtered_fasta(SAMPLE, q, c, RUN_TAG)
        for q, c in QUERY_CHROMOSOME_PAIRS
    ]


def final_genome_fasta():
    return f"{OUT}/final_genome/{SAMPLE}.aggregate.patch.filtered.{RUN_TAG}.fasta"


# Only computed on master (after RUN_TAG is defined in snakefile).
if not CHILD_MODE:
    TOTAL_HEAVY_JOBS = 2 * len(PANGENOME_NAMES) * len(QUERY_CHROMOSOME_PAIRS)
    WORKFLOW_MODE = "flat" if TOTAL_HEAVY_JOBS < FLAT_THRESHOLD else "nested"
    JOBS_PER_CHR = 2 * len(PANGENOME_NAMES) + 2
    CHR_THREADS = min(THREADS, max(1, JOBS_PER_CHR))
    print(
        f"[impuT2T] mode={WORKFLOW_MODE} "
        f"heavy_jobs={TOTAL_HEAVY_JOBS} (threshold={FLAT_THRESHOLD}) "
        f"pangenome_donors={len(PANGENOME_NAMES)} "
        f"chr_threads={CHR_THREADS}",
        file=sys.stderr,
    )
