_Updated: 2026-07-22_

# ImpuT2T

Snakemake workflow for **pangenome-guided patching** of diploid assemblies using the HPRC v2 panel.

RagTag assigns contigs to chromosomes; ImpuT2T then aligns each chromosome to many pangenome donors, chooses junctions and haplotypes, assembles patched scaffolds, and applies a post-hoc quality filter.

Detailed docs: [Wiki](https://github.com/maojanlin/ImpuT2T/wiki) (or `docs/` in this repo).

## Dependencies

Tools expected on `PATH`:

- **Snakemake** ≥ 7
- **RagTag** 2.0.1
- **minimap2**
- **samtools**
- **Python 3** with `numpy`, `scikit-learn`, `pysam`
- **[agc](https://github.com/refresh-bio/agc)** (needed for database setup)

Install notes and environment tips are in the [Wiki → Dependencies](https://github.com/maojanlin/ImpuT2T/wiki/Dependencies).

## Quick start

### 1. Download reference data

```bash
git clone https://github.com/maojanlin/ImpuT2T.git
cd ImpuT2T

# Smoke-test refs only (POP-10, chr20 + chr22) — for the quick test below:
database/setup_database.sh --test

# Or full panel (all chromosomes + donors) — for real runs:
database/setup_database.sh

```

### 2. Set path parameters

Edit `snakemake/config.yaml` (or use the shipped `config.POP-10.yaml` for the smoke test). The main path entries:

| Parameter | Meaning |
|-----------|---------|
| **`queries`** | Paths to your diploid haplotype assemblies (typically two FASTAs: maternal and paternal). Each file is scaffolded with RagTag, then patched. |
| **`target_chromosomes`** | Which chromosomes to process (e.g. `["chr1", …, "chr22", "chrX"]`). Smoke test uses `["chr20", "chr22"]`. |
| **`sample`** | Sample / run label used in output filenames under `output/` (e.g. `CN1`, `CN1-test`). |

`config.POP-10.yaml` already points at `testdata/` and `database/` — no edits needed for the smoke test. More options: [Wiki → Configuration](https://github.com/maojanlin/ImpuT2T/wiki/Configuration).

### 3. Run (smoke test)

```bash
cd snakemake
snakemake --cores 48 --keep-going --rerun-incomplete --latency-wait 120 \
  --configfile config.POP-10.yaml
```

### Check results

Primary outputs under `snakemake/output_test/final_genome/`:

| File | Content |
|------|---------|
| `CN1-test.aggregate.patch.filtered.test-POP-10.fasta` | Patched sequences |
| `CN1-test.aggregate.patch.filtered.test-POP-10.agp` | AGP layout |
| `CN1-test.aggregate.patch.filtered.test-POP-10.edge.log` | Junction scores and probabilities |

When the smoke test succeeds, run `database/setup_database.sh` (full panel if you used `--test` earlier), set **`queries`** / **`target_chromosomes`** / **`sample`** in `config.yaml`, and:

```bash
cd snakemake
snakemake --cores 48 --keep-going --rerun-incomplete --latency-wait 120 \
  --configfile config.yaml
```
