# Wiki home

These pages are meant to be copied into the GitHub Wiki
(`https://github.com/maojanlin/ImpuT2T/wiki`). The repo README stays a short
1–2–3 quick start; details live here.

Suggested Wiki sidebar / pages:

1. [Home](Home.md) (this page)
2. [Dependencies](Dependencies.md)
3. [Reference-data](Reference-data.md)
4. [Configuration](Configuration.md)
5. [Running](Running.md)
6. [Outputs](Outputs.md)

## Create the Wiki on GitHub

1. Open the repository on GitHub → **Wiki** tab → **Create the first page**.
2. Paste content from `docs/*.md` into matching Wiki pages (filenames become page titles).
3. Optional: clone the wiki as a separate git repo:

```bash
git clone https://github.com/maojanlin/ImpuT2T.wiki.git
# edit .md files, then commit and push
```

## Quick path (from README)

```bash
git clone https://github.com/maojanlin/ImpuT2T.git
cd ImpuT2T
database/setup_database.sh
cd snakemake
snakemake --cores 48 --keep-going --rerun-incomplete --latency-wait 120 \
  --configfile config.POP-10.yaml
```
