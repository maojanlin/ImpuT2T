# Running

## Smoke test (recommended first)

```bash
database/setup_database.sh --test
cd snakemake
snakemake --cores 48 --keep-going --rerun-incomplete --latency-wait 120 \
  --configfile config.POP-10.yaml
```

Shipped queries: `testdata/CN1.hap1.chr20.contigs.fa` and `testdata/CN1.hap2.chr22.contigs.fa` (~6.4 Mb).

## Full panel

```bash
database/setup_database.sh
# edit snakemake/config.yaml (queries, sample, …)
cd snakemake
snakemake --cores 48 --keep-going --rerun-incomplete --latency-wait 120 \
  --configfile config.yaml
```

## Options

| Flag | When to use |
|------|-------------|
| `--use-conda` | RagTag not on `PATH`; uses `envs/ragtag.yaml` for scaffolding only |
| `--resources nested_workflow=2` | Limit concurrent per-chromosome child Snakemake jobs |
| `--dry-run` / `-n` | Print the DAG without running |

**Mode selection:** if `2 × n_donors × n_chromosome_pairs < heavy_jobs_flat_threshold`, **flat** mode runs in one process; otherwise **nested** mode (one child Snakemake per chromosome). Snakemake prints the chosen mode at startup.
