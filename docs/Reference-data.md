# Reference data

All reference files go under `database/` at the repository root (gitignored downloads; empty placeholders are tracked).

```
ImpuT2T/
  database/
    human472.agc
    high_quality_set/      # full genomes → small_pan (RagTag)
    chr_agc/               # per-chromosome AGC archives
    chr_split/             # per-donor per-chr FASTAs → pangenome_path
  snakemake/
```

## One-command setup

```bash
# From repository root
database/setup_database.sh              # smoke: POP-10, chr20 + chr22 only
database/setup_database.sh --full       # all chromosomes + full donor panel
```

The script:

1. Downloads [human472.agc](https://zenodo.org/records/14854401) and extracts the default `small_pan` genomes into `high_quality_set/`.
2. Downloads HPRC v2 chromosome AGC archives from the [Index zone](https://benlangmead.github.io/aws-indexes/) (`https://genome-idx.s3.amazonaws.com/agc/hprc2/pangenome_<chr>.agc`).
3. Calls `database/extract_chr_pangenome.sh` to write `chr_split/<chr>/<donor>.fasta`.

Requires `agc` and `wget` (or `curl`) on `PATH`.

## Manual / partial steps

Extract only (archives already downloaded):

```bash
database/extract_chr_pangenome.sh snakemake/subsample_lists/sample_id.POP-10.log
database/extract_chr_pangenome.sh snakemake/subsample_lists/sample_id.log
```

Missing chromosome AGC files are skipped. Donor IDs in the list must match sample names inside each archive (e.g. `CHM13`, `HG002.1`).
