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

## Targeted tiebreaking

The default tiebreak settings are:

```yaml
tiebreak_backend: targeted
tiebreak_cache_mode: reuse
tiebreak_dump_info: "off"
```

`targeted` extracts only the terminal sequences needed by ambiguous edges selected by the greedy path, grouped by their tied donors. It retains full alignment context for connection selection and uses the same 25 kb map-length parameter and scoring rules as the legacy implementation. Full PAF/reference reads and reference indexing can still dominate a small batch.

Cache policies are `reuse` (reuse validated results and compute misses), `refresh` (recompute requested results), and `only` (no alignment; fail if any requested result is missing, stale, or unavailable). Completed no-connection results are cached; execution failures are not. Existing legacy `.info` files are not imported as trusted cache records.

Set `tiebreak_dump_info: "on"` to export targeted four-column donor `.info` files plus `evidence.tsv` and `info.metadata.json`. This exports cached results too and never adds alignment work. These tables contain only the current requested edges; the metadata identifies the current donor set. Missing evidence is recorded in the audit table, not represented by a fabricated score. Previously exported files are not deleted when export is later disabled.

The equivalent aggregation CLI options are `--tiebreak_backend targeted`, `--tiebreak_cache_mode reuse`, and `--tiebreak_dump_info off`. `--tiebreak_dump_info` alone enables export. The existing `--tiebreak_no_reuse` flag is an alias for refresh and conflicts with cache-only mode. `--tiebreak_jobs` and `--tiebreak_threads` still control donor parallelism and threads per aligner.

Use `tiebreak_backend: legacy` to run the donor-wide implementation. Legacy mode always writes its required `.info` files; export-off only controls optional targeted exports. Strict cache-only mode requires the targeted backend.

Results are stored under `<tiebreak_dir>/targeted/cache.sqlite3`, with optional tables under `targeted/info/`. Cache identity includes content digests for the PAF, FASTAs, existing query index, scorer/dependencies, minimap2 binary/version, and parameters. File digests are reused while file identity, size, modification time and change time remain unchanged. Inputs must remain unchanged during scoring; detected changes abort the batch without caching it. Keep minimap2 on PATH even for cache-only mode so its identity can be validated. Temporary terminal FASTA/PAF files are cleaned up after each donor batch. `tiebreak_resolutions.tsv` remains the final-choice audit log.

Settings propagate to nested chromosome jobs as well as flat workflows. Use a fresh output directory for validation; updating the original repository does not update existing copied experiment workflows.
