# Dependencies

## Required on `PATH`

| Tool | Recommended version | Used for |
|------|---------------------|----------|
| [Snakemake](https://snakemake.github.io/) | ≥ 7 (8.x) | Workflow engine |
| [RagTag](https://github.com/malonge/RagTag) | 2.0.1 | Chromosome scaffolding (`ragtag.py`) |
| [minimap2](https://github.com/lh3/minimap2) | ≥ 2.24 | Pangenome alignment |
| [samtools](https://www.htslib.org/) | ≥ 1.17 | FASTA indexing |
| [agc](https://github.com/refresh-bio/agc) | ≥ 3.2 (bioconda) | Extracting AGC archives |
| Python 3 | ≥ 3.9 | Most rules call `python3 scripts/...` |

## Python packages

```bash
pip install numpy scikit-learn pysam
```

`ragtag_utilities` is bundled under `snakemake/scripts/ragtag_utilities/`.

## RagTag

Recommended: install RagTag yourself so it is on `PATH`:

```bash
conda install -c bioconda ragtag=2.0.1
```

**Fallback:** the workflow ships `snakemake/envs/ragtag.yaml`. Pass `--use-conda` and Snakemake builds that env for the scaffolding rule only:

```bash
snakemake --cores 48 --use-conda --configfile config.POP-10.yaml
```

`minimap2`, `samtools`, and the Python packages are still required on `PATH` either way; `--use-conda` does not install those.
