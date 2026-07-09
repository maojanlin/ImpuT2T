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
git clone <repository-url> impuT2T_workflow
cd impuT2T_workflow
# Download reference data into database/ (see below)
cd snakemake
# Edit config.yaml (query FASTA paths) and run
snakemake -n --configfile config.yaml
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

After cloning the repository, download reference files into `database/` at the repo root. The subfolders `high_quality_set/`, `chr_agc/`, and `chr_split/` are included in git as empty placeholders; downloaded files are gitignored. Paths in `snakemake/config.yaml` use `../database/`.

```
ImpuT2T/          # git clone root
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

`config.yaml` points `small_pan` at `../database/high_quality_set/*.fa`. Remove CN1 or other entries from both the config and the extract loop if you do not need them.

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

Donor IDs must match `pangenome_names` in the config and `SAMPLE_LIST_FILE` in `snakefile`. Not every donor has every chromosome (chrX, chrY, chrM may be partial); missing files are skipped automatically.

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

Or with AWS CLI:

```bash
for chr in $(seq 1 22 | sed 's/^/chr/') chrX chrY chrM; do
  aws s3 cp --no-sign-request \
    "s3://genome-idx/agc/hprc2/pangenome_${chr}.agc" \
    "database/chr_agc/pangenome_${chr}.agc"
done
```

**2. Extract donors** into `database/chr_split/` with [AGC](https://github.com/refresh-bio/agc) (`conda install -c bioconda agc`). Sample names inside each archive match `pangenome_names` (e.g. `CHM13`, `HG002.1`):

```bash
# Default: donors from snakemake/subsample_lists/sample_id.POP-10.log
database/extract_chr_pangenome.sh

# Or pass a donor list (one ID per line, same as SAMPLE_LIST_FILE)
database/extract_chr_pangenome.sh snakemake/subsample_lists/sample_id.POP-200.log
```

Manual single-donor example:

```bash
agc getset database/chr_agc/pangenome_chr1.agc CHM13 > database/chr_split/chr1/CHM13.fasta
```

## Configuration

The sample config is `snakemake/config.yaml`. Edit it in place, or copy it for a new run:

```bash
cd snakemake
cp config.yaml config.my_run.yaml
```

Larger pre-built configs (full HPRC panel, 200-donor subsamples, etc.) are also in `snakemake/`: `config.HG002.yaml`, `config.CN1.200.yaml`, `config.CN1.50.yaml`, `config.PAN027.yaml`, …

### Keys to change

| Key | Description |
|-----|-------------|
| `queries` | `[haplotype_name, path/to/query.asm.fa]` — replace `path/to/query.hap*.fa` |
| `small_pan` | `[name, path/to/reference.fa]` — default: `../database/high_quality_set/` |
| `small_pan_prefixes` | AGP ID prefixes matching `small_pan` entries |
| `target_chromosomes` | Per-haplotype chromosome lists (autosomes ± X/Y/M) |
| `pangenome_path` | Root of chromosome-split donors — default: `../database/chr_split` |
| `pangenome_names` | Donor IDs (should match `SAMPLE_LIST_FILE` and files under `pangenome_path`) |
| `sample` | Sample label used in output paths |
| `output_dir` | Directory for all run artifacts (under `snakemake/`) |
| `threads` | Total cores for Snakemake |
| `ragtag_threads` | Threads per RagTag job (max 3) |
| `heavy_jobs_flat_threshold` | Auto flat vs nested mode (default `6000`) |

### Edit `snakemake/snakefile`

These values are not in the YAML; defaults in the shipped `snakefile` match `config.yaml`:

```python
configfile: "config.yaml"
RUN_TAG = "run0"
SAMPLE_LIST_FILE = "subsample_lists/sample_id.log"
```

`SAMPLE_LIST_FILE` is one donor ID per line and must agree with `pangenome_names`. Pre-made lists are in `snakemake/subsample_lists/` (`POP-5`, `POP-10`, `POP-50`, `POP-200`, `EUR-*`, `AFR-*`, …). For a larger panel, switch both `pangenome_names` and `SAMPLE_LIST_FILE` together (e.g. use `config.CN1.200.yaml` + `sample_id.POP-200.log` as a reference).

## Running

```bash
cd snakemake

# Dry run
snakemake --dry-run --configfile config.yaml

# Full run (RagTag already on PATH)
snakemake  --cores 48 --keep-going --latency-wait 120 --configfile config.yaml
```

For large pangenomes, the workflow runs in **nested** mode (one child Snakemake per chromosome). To limit concurrent chromosome jobs, add `--resources nested_workflow=2` to either command above.

**Mode selection:** if `2 × n_donors × n_chromosome_pairs < heavy_jobs_flat_threshold`, **flat** mode runs everything in one process; otherwise **nested** mode is used. Snakemake prints the chosen mode at startup.

## Outputs

| Path | Description |
|------|-------------|
| `{output_dir}/impuT2T_patch/{RUN_TAG}/*.aggregate.patch.filtered.fasta` | Per-chromosome filtered patched sequences |
| `{output_dir}/final_genome/{sample}.aggregate.patch.filtered.{RUN_TAG}.fasta` | Merged genome across chromosomes |
| `{output_dir}/final_genome/{sample}.aggregate.patch.filtered.{RUN_TAG}.agp` | Merged AGP |
| `{output_dir}/ragtag_output/`, `consensus_scaffold/`, `contigs_assignment/` | Intermediate prepare-stage outputs |

Logs and benchmarks are written under `{output_dir}/logs/` and `{output_dir}/benchmarks/`.

## Workflow overview

1. **Prepare** — RagTag scaffold against `small_pan`; consensus chromosome assignment; per-chromosome contig FASTAs.
2. **Pangenome** — minimap2 alignment per donor; local connection inference; aggregate layout, haplotype selection, tiebreak; post-patch quality filter.
3. **Final** — concatenate per-chromosome filtered FASTAs into one genome.
