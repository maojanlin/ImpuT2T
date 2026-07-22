# Reference database

Empty folders `high_quality_set/`, `chr_agc/`, and `chr_split/` are tracked as placeholders.
Downloads are gitignored.

**Recommended:** from the repository root, run:

```bash
database/setup_database.sh              # smoke test (POP-10, chr20 + chr22)
database/setup_database.sh --full       # full panel
```

This downloads Zenodo `human472.agc` and Index-zone `pangenome_chr*.agc`, then calls
`extract_chr_pangenome.sh`. See the [Wiki → Reference data](https://github.com/maojanlin/ImpuT2T/wiki/Reference-data)
or `docs/Reference-data.md`.
