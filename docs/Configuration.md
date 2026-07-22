# Configuration

`run_tag` and `sample_list_file` live in the **config YAML** (not the snakefile). Switching panels is only `--configfile`.

| Config | Purpose |
|--------|---------|
| `config.POP-10.yaml` | Smoke test: POP-10 donors + shipped CN1 contig subsets |
| `config.yaml` | Full HPRC panel; set your own `queries` / `sample` |

```bash
cd snakemake
cp config.yaml config.my_run.yaml
```

## Keys

**Required for a production run:** `queries`, `target_chromosomes`, and `sample`. The smoke-test config already sets these.

| Key | Description |
|-----|-------------|
| **`queries`** | `[haplotype_name, path/to/query.asm.fa]` |
| `small_pan` | `[name, path/to/reference.fa]` — default `../database/high_quality_set/` |
| `small_pan_prefixes` | AGP ID prefixes matching `small_pan` |
| **`target_chromosomes`** | Per-haplotype chromosome lists (e.g. maternal without chrY) |
| `pangenome_path` | Default `../database/chr_split` |
| `pangenome_names` | Must match `sample_list_file` and files under `pangenome_path` |
| **`sample`** | Label in output paths |
| `run_tag` | Tag for patch / final outputs |
| `sample_list_file` | Donor list under `subsample_lists/` |
| `output_dir` | Run artifacts under `snakemake/` |
| `threads` / `ragtag_threads` | CPU allocation |
| `heavy_jobs_flat_threshold` | Flat vs nested mode (default `6000`) |

When changing the donor panel, update **`pangenome_names`** and **`sample_list_file`** together (and use a new `run_tag` / `output_dir`).
