_Updated: 2026-07-09_

# ImpuT2T

This is a Snakemake workflow for ImpuT2T, a pipeline for pangenome-guided patching of diploid genome assemblies using the HPRC v2 panel.

First, the workflow assigns input contigs to chromosomes. This is done by running RagTag scaffolding on a high-quality subset of pangenome references (“small_pan”) and reaching consensus assignments.

Next, for each chromosome, the patching module performs several steps:
- Aligns the chromosome-partitioned contigs to all pangenome donor references.
- Selects which connections and donor haplotypes to use.
- Assembles new, patched scaffolds from these selections.
- Applies quality filters to the results.

This modular approach makes each step clear and allows for quality control at key stages.

All workflow files live under `snakemake/`. Run Snakemake from that directory. Reference databases are downloaded into `database/` at the **repository root** (not tracked by git).

## Quick start

```bash
git clone https://github.com/maojanlin/ImpuT2T.git
cd ImpuT2T
# Download reference data into database/ (see below)
# Prefer the smoke test first (see Test run), then a full-panel run with config.yaml
cd snakemake
snakemake --dry-run --configfile config.POP-10.yaml
```

## Dependencies

### Required on `PATH`

| Tool | Recommended version | Used for |
|------|---------------------|----------|
| [Snakemake](https://snakemake.github.io/) | ≥ 7 (8.x) | Workflow engine (`--configfile` on CLI) |
| [RagTag](https://github.com/malonge/RagTag) | 2.0.1 | Chromosome scaffolding (`ragtag.py`; see below) |
| [minimap2](https://github.com/lh3/minimap2) | ≥ 2.24 | Pangenome alignment |
| [samtools](https://www.htslib.org/) | ≥ 1.17 | FASTA indexing |
| Python 3 | ≥ 3.9 | Most rules call `python3 scripts/...` directly |

### Python packages

Install on the machine or environment where Snakemake runs (not managed by a workflow-wide conda env yet):

```bash
pip install numpy scikit-learn pysam
```

`ragtag_utilities` is bundled under `snakemake/scripts/ragtag_utilities/`.

### RagTag

Only the **RagTag scaffolding** step (`run_ragtag`) calls `ragtag.py`. We recommend installing RagTag yourself so it is on your `PATH`:

```bash
conda install -c bioconda ragtag=2.0.1
```

Alternatively, you can rely on the workflow `snakemake/envs/ragtag.yaml` (RagTag 2.0.1). Pass `--use-conda` and Snakemake will create that env for the scaffolding rule only:

```bash
snakemake -j 48 --use-conda --configfile config.yaml
```

## Reference data

After cloning the repository, download reference files into `database/` at the repo root. The subfolders `high_quality_set/`, `chr_agc/`, and `chr_split/` are included in git as empty placeholders. Paths in `snakemake/config.yaml` use `../database/`.

```
ImpuT2T/                   # git clone root
  database/                # you download references here
    human472.agc
    high_quality_set/      # full-genome FASTAs → small_pan
    chr_agc/               # per-chromosome AGC archives (download)
    chr_split/             # per-donor per-chr FASTAs → pangenome_path (agc extract)
      chr1/
        CHM13.fasta
        ...
  snakemake/
    config.yaml
```

### 1. Small pangenome FASTAs (RagTag scaffolding)

Full-genome references listed in `small_pan` in `config.yaml` (CHM13, GRCh38, CN1, YAO, KSA001, HG002 by default).

These assemblies are in the [human472 AGC archive](https://zenodo.org/records/14854401) on Zenodo. Download the archive (~3.4 GB) and extract only the genomes you list in `small_pan`:

```bash
# From the repository root (high_quality_set/ already exists after git clone)

# Install a recent AGC (human472.agc requires bioconda 3.2.x or newer)
conda install -c bioconda agc

# Download
wget -O database/human472.agc \
  https://zenodo.org/records/14854401/files/human472.agc?download=1

OUTDIR=database/high_quality_set

# List sample names inside the archive (internal IDs, no .fa suffix)
agc listset database/human472.agc

# Extract the default small_pan set (matches config.yaml)
for id in \
  100000_CHM13.pri \
  110000_GRCh38.pri \
  120001_CN1.pat \
  120002_CN1.mat \
  120003_YAO.pat \
  120004_YAO.mat \
  120005_KSA001.pat \
  120006_KSA001.mat \
  200001_HG002.pat \
  200002_HG002.mat
do
  agc getset database/human472.agc "$id" > "${OUTDIR}/${id}.fa"
done
```

`config.yaml` points `small_pan` at `../database/high_quality_set/*.fa`.

### 2. Chromosome-partitioned HPRC v2 pangenome (patching)

Per-chromosome donor FASTAs for minimap2 alignment and local patching. The workflow expects one FASTA per donor under `database/chr_split/<chromosome>/` (default `pangenome_path`: `../database/chr_split`):

```
database/chr_split/
  chr1/
    CHM13.fasta
    GRCh38.fasta
    HG002.1.fasta
    ...
  chr2/
    ...
```

Donor IDs should match `pangenome_names` in the config and `SAMPLE_LIST_FILE` in `snakefile`. Not every donor has every chromosome (chrX, chrY, chrM may be partial); missing files are skipped automatically.

### Download (chromosome-partitioned AGC archives)

Twenty-five per-chromosome AGC archives (chr1–chr22, chrX, chrY, chrM) are on the [Langmead Lab Index zone](https://benlangmead.github.io/aws-indexes/) (AWS Open Data; free HTTPS/S3 access):

| Chromosome | URL |
|------------|-----|
| chr1–chr22, chrX, chrY, chrM | `https://genome-idx.s3.amazonaws.com/agc/hprc2/pangenome_<chr>.agc` |

Example for chr1: [pangenome_chr1.agc](https://genome-idx.s3.amazonaws.com/agc/hprc2/pangenome_chr1.agc)

**1. Download archives** into `database/chr_agc/` (from the repository root):

```bash
BASE=https://genome-idx.s3.amazonaws.com/agc/hprc2
for chr in $(seq 1 22 | sed 's/^/chr/') chrX chrY chrM; do
  wget -c -O "database/chr_agc/pangenome_${chr}.agc" \
    "${BASE}/pangenome_${chr}.agc"
done
```

<!--
Or with AWS CLI:

```bash
for chr in $(seq 1 22 | sed 's/^/chr/') chrX chrY chrM; do
  aws s3 cp --no-sign-request \
    "s3://genome-idx/agc/hprc2/pangenome_${chr}.agc" \
    "database/chr_agc/pangenome_${chr}.agc"
done
```
-->

**2. Extract donors** into `database/chr_split/` with [AGC](https://github.com/refresh-bio/agc) (`conda install -c bioconda agc`). Sample names inside each archive match `pangenome_names` (e.g. `CHM13`, `HG002.1`):

```bash
# Default: donors from snakemake/subsample_lists/sample_id.POP-10.log
database/extract_chr_pangenome.sh

# Or pass a donor list (one ID per line, same as SAMPLE_LIST_FILE)
database/extract_chr_pangenome.sh snakemake/subsample_lists/sample_id.POP-10.log
```

<!--
Manual single-donor example:

```bash
agc getset database/chr_agc/pangenome_chr1.agc CHM13 > database/chr_split/chr1/CHM13.fasta
```
-->

## Configuration

`run_tag` and `sample_list_file` live in the **config YAML** (not the snakefile). Switching between a smoke test and a full-panel run is only a different `--configfile` — no snakefile edits.

| Config | Purpose |
|--------|---------|
| `config.yaml` | Full HPRC panel; set your own `queries` / `sample` |
| `config.POP-10.yaml` | Smoke test: POP-10 donors + shipped CN1 contig subsets |

```bash
cd snakemake
cp config.yaml config.my_run.yaml   # optional: personalize a full-panel run
```


### Keys to change

**Required for a production run:** `queries`, `target_chromosomes`, and `sample` (bold below). For the smoke test, `config.POP-10.yaml` already sets these.

| Key | Description |
|-----|-------------|
| **`queries`** | `[haplotype_name, path/to/query.asm.fa]` — replace `path/to/query.hap*.fa` |
| `small_pan` | `[name, path/to/reference.fa]` — default: `../database/high_quality_set/` |
| `small_pan_prefixes` | AGP ID prefixes matching `small_pan` entries |
| **`target_chromosomes`** | Lists of chromosomes for each haplotype (autosomes ± X/Y/M). For example, a maternal haplotype should exclude chrY from its list. |
| `pangenome_path` | Root of chromosome-split donors — default: `../database/chr_split` |
| `pangenome_names` | Donor IDs (should match `sample_list_file` and files under `pangenome_path`) |
| **`sample`** | Sample label used in output paths |
| `run_tag` | Tag for patch / final output paths (e.g. `run0`, `test-POP-10`) |
| `sample_list_file` | Donor list under `subsample_lists/` (one ID per line) |
| `output_dir` | Directory for all run artifacts (under `snakemake/`) |
| `threads` | Total cores for Snakemake |
| `ragtag_threads` | Threads per RagTag job (max 3) |
| `heavy_jobs_flat_threshold` | Auto flat vs nested mode (default `6000`) |

When changing the donor panel, update **`pangenome_names`** and **`sample_list_file`** together (and use a new `run_tag` / `output_dir` so outputs do not collide).

## Test run (POP-10)

Use this before a full-panel run. It exercises the whole pipeline on two chromosomes and 10 donors (flat mode).

**Shipped query data** (`testdata/`, ~6.4 Mb total):

| File | Contigs |
|------|---------|
| `CN1.hap1.chr20.contigs.fa` | 4 contigs (partial chr20), ~1.5 Mb |
| `CN1.hap2.chr22.contigs.fa` | 3 contigs (partial chr22), ~5.0 Mb |

**1. Extract only POP-10 donors** (after downloading `database/chr_agc/`):

```bash
# From the repository root
database/extract_chr_pangenome.sh snakemake/subsample_lists/sample_id.POP-10.log
```

Also ensure `database/high_quality_set/` is populated (RagTag references; see above).

**2. Dry run and execute** (no snakefile changes):

```bash
cd snakemake
snakemake --dry-run --configfile config.POP-10.yaml
snakemake --cores 48 --keep-going --rerun-incomplete --latency-wait 120 \
  --configfile config.POP-10.yaml
```

**3. Check results** under `snakemake/output_test/final_genome/` (`CN1-test.aggregate.patch.filtered.test-POP-10.*`).

When the smoke test succeeds, switch to a full-panel config (edit `queries` / `sample` in `config.yaml` or a copy), extract the full donor list into `chr_split/`, and run with `--configfile config.yaml`.

## Running

```bash
cd snakemake

# Dry run (full panel — after editing queries in config.yaml)
snakemake --dry-run --configfile config.yaml

# Full run (RagTag already on PATH)
snakemake --cores 48 --keep-going --rerun-incomplete --latency-wait 120 \
  --configfile config.yaml
```

For large pangenomes, the workflow runs in **nested** mode (one child Snakemake per chromosome). To limit concurrent chromosome jobs, add `--resources nested_workflow=2` to either command above.

**Mode selection:** if `2 × n_donors × n_chromosome_pairs < heavy_jobs_flat_threshold`, **flat** mode runs everything in one process; otherwise **nested** mode is used. Snakemake prints the chosen mode at startup.

## Outputs

### Main results

These three files under `{output_dir}/final_genome/` are the primary deliverables:

| Path | Description |
|------|-------------|
| `{sample}.aggregate.patch.filtered.{RUN_TAG}.fasta` | Patched diploid genome (all chromosomes merged) |
| `{sample}.aggregate.patch.filtered.{RUN_TAG}.agp` | AGP for the merged genome |
| `{sample}.aggregate.patch.filtered.{RUN_TAG}.edge.log` | Junction edge log: inter-contig connections, junction scores, logistic-regression probabilities, and used/unused flags |

### Intermediate files

| Path | Description |
|------|-------------|
| `{output_dir}/impuT2T_patch/{RUN_TAG}/*.aggregate.patch.filtered.fasta` | Per-chromosome filtered patched sequences |
| `{output_dir}/impuT2T_patch/{RUN_TAG}/*.aggregate.edge.log` | Per-chromosome edge log (pre-filter) |
| `{output_dir}/impuT2T_patch/{RUN_TAG}/*.aggregate.patch.filtered.edge.log` | Per-chromosome edge log after post-patch quality filter |
| `{output_dir}/ragtag_output/`, `consensus_scaffold/`, `contigs_assignment/` | Prepare-stage outputs |

Logs and benchmarks are written under `{output_dir}/logs/` and `{output_dir}/benchmarks/`.

## Workflow overview

1. **Prepare** — RagTag scaffold against `small_pan`; consensus chromosome assignment; per-chromosome contig FASTAs.
2. **Pangenome** — minimap2 alignment per donor; local connection inference; aggregate layout, haplotype selection, tiebreak; post-patch quality filter.
3. **Final** — concatenate per-chromosome filtered FASTAs into one genome.
