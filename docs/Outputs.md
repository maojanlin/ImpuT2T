# Outputs

## Main results

Under `{output_dir}/final_genome/`:

| File | Content |
|------|---------|
| `{sample}.aggregate.patch.filtered.{RUN_TAG}.fasta` | Patched diploid genome |
| `{sample}.aggregate.patch.filtered.{RUN_TAG}.agp` | AGP for the merged genome |
| `{sample}.aggregate.patch.filtered.{RUN_TAG}.edge.log` | Junction edge log (scores, probabilities, used/unused) |

Smoke-test example: `snakemake/output_test/final_genome/CN1-test.aggregate.patch.filtered.test-POP-10.*`.

## Intermediate files

| Path | Description |
|------|-------------|
| `{output_dir}/impuT2T_patch/{RUN_TAG}/*.aggregate.patch.filtered.fasta` | Per-chromosome filtered patches |
| `{output_dir}/impuT2T_patch/{RUN_TAG}/*.aggregate.edge.log` | Per-chromosome edge log (pre-filter) |
| `{output_dir}/impuT2T_patch/{RUN_TAG}/*.aggregate.patch.filtered.edge.log` | Edge log after post-patch QC |
| `{output_dir}/ragtag_output/`, `consensus_scaffold/`, `contigs_assignment/` | Prepare stage |

Logs and benchmarks: `{output_dir}/logs/`, `{output_dir}/benchmarks/`.

## Workflow overview

1. **Prepare** — RagTag against `small_pan`; consensus chromosome assignment.
2. **Pangenome** — minimap2; local connections; aggregate layout; haplotype selection; post-patch filter.
3. **Final** — merge per-chromosome products into one genome.
