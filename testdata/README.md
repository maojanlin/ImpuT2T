# Smoke-test query assemblies

CN1 contig subsets used by `snakemake/config.POP-10.yaml`. Contigs are original
assembly IDs corresponding to selected `*_impuT2T_*` blocks from a prior
chromosome assignment (not the renamed headers).

| File | Contigs (impuT2T IDs) | Approx. size |
|------|------------------------|--------------|
| `CN1.hap1.chr20.contigs.fa` | hap1 chr20: IDs **3, 4, 5, 6** (`h1tg000334l`, `h1tg000115l`, `h1tg000201l`, `h1tg000105l`) | ~1.5 Mb |
| `CN1.hap2.chr22.contigs.fa` | hap2 chr22: IDs **6, 7, 9** (`h2tg000160l`, `h2tg000126l`, `h2tg000069l`) | ~5.0 Mb |

Total ~6.4 Mb — small enough to keep in git. Enough to exercise RagTag
assignment, multi-donor alignment, patching across gaps, and final merge.
