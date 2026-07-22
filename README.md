_Updated: 2026-07-22_

# ImpuT2T

Snakemake workflow for **pangenome-guided patching** of diploid assemblies using the HPRC v2 panel.

RagTag assigns contigs to chromosomes; ImpuT2T then aligns each chromosome to many pangenome donors, chooses junctions and haplotypes, assembles patched scaffolds, and applies a post-hoc quality filter.

Detailed docs (dependencies, full-panel config, nested mode, outputs): see the [Wiki](https://github.com/maojanlin/ImpuT2T/wiki) (or `docs/` in this repo).

## Quick start

**Prerequisites on `PATH`:** Snakemake ≥ 7, RagTag 2.0.1, minimap2, samtools, Python 3 (`numpy`, `scikit-learn`, `pysam`), and [agc](https://github.com/refresh-bio/agc) (`conda install -c bioconda agc`). See the Wiki for install details.

### 1. Download reference data

```bash
git clone https://github.com/maojanlin/ImpuT2T.git
cd ImpuT2T
database/setup_database.sh          # smoke-test refs (POP-10, chr20 + chr22)
# database/setup_database.sh --full # later: all chromosomes + full donor panel
```

### 2. Set path parameters

For the **shipped smoke test**, nothing to edit — `snakemake/config.POP-10.yaml` already points at `testdata/` and `database/`.

For **your own assembly**, copy and edit `snakemake/config.yaml` (at least **`queries`**, **`target_chromosomes`**, and **`sample`**). Details: [Wiki → Configuration](https://github.com/maojanlin/ImpuT2T/wiki/Configuration).

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

When the smoke test succeeds, run `database/setup_database.sh --full`, set your query paths in `config.yaml`, and:

```bash
cd snakemake
snakemake --cores 48 --keep-going --rerun-incomplete --latency-wait 120 \
  --configfile config.yaml
```
