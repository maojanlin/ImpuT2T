# Shared configuration and path helpers (included by master and chromosome workflows).

import os
import sys

OUT = config["output_dir"]

# Paths in included .smk files resolve relative to workflow/, not the repo root.
WORKFLOW_ROOT = workflow.basedir
if os.path.basename(WORKFLOW_ROOT.rstrip("/")) == "workflow":
    WORKFLOW_ROOT = os.path.dirname(WORKFLOW_ROOT.rstrip("/"))
SCRIPTS = os.path.join(WORKFLOW_ROOT, "scripts")
RAGTAG_CONDA_ENV = os.path.join(WORKFLOW_ROOT, "envs", "ragtag.yaml")
CHROMOSOME_SMK = os.path.join(WORKFLOW_ROOT, "workflow", "chromosome.smk")

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
ALL_PANGENOME_NAMES = PANGENOME_NAMES
PANGENOME_PATH = config["pangenome_path"].rstrip("/")

FLAT_THRESHOLD = config.get("heavy_jobs_flat_threshold", 6000)

# Child chromosome workflow passes scope via --config; master sets globals in snakefile.
CHILD_MODE = (
    config.get("run_tag") is not None
    and config.get("query_name") is not None
    and config.get("chromosome") is not None
)


def pangenome_fasta_path(chromosome, name):
    return os.path.join(PANGENOME_PATH, chromosome, f"{name}.fasta")


def pangenome_names_for(chromosome):
    """Donors that have a per-chromosome split FASTA (chrX/chrY/chrM are often partial)."""
    return [
        name
        for name in ALL_PANGENOME_NAMES
        if os.path.isfile(pangenome_fasta_path(chromosome, name))
    ]


# run_tag / sample_list_file live in the YAML so switching panels is --configfile only.
# Child jobs still pass them via --config (overrides YAML for that invocation).
if "run_tag" not in config:
    raise KeyError(
        "config must set run_tag (e.g. in config.yaml or config.POP-10.yaml)"
    )
if "sample_list_file" not in config:
    raise KeyError(
        "config must set sample_list_file (path relative to snakemake/)"
    )
RUN_TAG = config["run_tag"]
SAMPLE_LIST_FILE = config["sample_list_file"]

if CHILD_MODE:
    QUERY_CHROMOSOME_PAIRS = [(config["query_name"], config["chromosome"])]
    PANGENOME_NAMES = pangenome_names_for(config["chromosome"])
    skipped = len(ALL_PANGENOME_NAMES) - len(PANGENOME_NAMES)
    if skipped:
        print(
            f"[impuT2T] {config['chromosome']}: "
            f"{len(PANGENOME_NAMES)}/{len(ALL_PANGENOME_NAMES)} pangenome donors "
            f"(skipped {skipped} without split FASTA)",
            file=sys.stderr,
        )
else:
    QUERY_CHROMOSOME_PAIRS = [
        (q, c) for q in QUERY_NAMES for c in QUERY_TO_CHROMOSOMES[q]
    ]
    PANGENOME_NAMES_BY_CHROMOSOME = {
        chromosome: pangenome_names_for(chromosome)
        for chromosome in {c for _, c in QUERY_CHROMOSOME_PAIRS}
    }


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


# Only computed on master (after RUN_TAG / SAMPLE_LIST_FILE from config).
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
