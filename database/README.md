# Reference database (local download)

The subfolders below are tracked in git as empty placeholders. After cloning, download and extract reference data here (see the [main README](../README.md)).

```
database/
  human472.agc              # Zenodo AGC archive → extract to high_quality_set/
  high_quality_set/         # full-genome FASTAs for RagTag (small_pan)
  chr_agc/                  # Index-zone per-chromosome archives (pangenome_chr*.agc)
  chr_split/                # per-donor per-chromosome FASTAs for patching (pangenome_path)
    chr1/
      CHM13.fasta
      ...
```

Index zone AGC URLs: `https://genome-idx.s3.amazonaws.com/agc/hprc2/pangenome_<chr>.agc`  
(chr1–chr22, chrX, chrY, chrM — 25 files)

Extract chr_split FASTAs: `database/extract_chr_pangenome.sh [donor_list]`

Paths in `snakemake/config.yaml` point here via `../database/`.
