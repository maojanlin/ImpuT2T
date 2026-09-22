# Targeted tiebreak validation

Run with the workflow's Python dependencies (pysam, numpy, sklearn, networkx) and minimap2 on PATH:

```bash
python -m unittest discover -s tests -v
```

To compare with an independent, unmodified workflow checkout as well as the legacy backend:

```bash
IMPUT2T_BASELINE_SCRIPT=/path/to/original/snakemake/scripts/from_paf_to_multi_connections.py \
  python -m unittest discover -s tests -v
```

The optional baseline also supplies the sibling original aggregate script for end-to-end AGP/edge-log parity. Use the same minimap2 binary/version for both implementations. Tests cover reverse and short terminals, split/duplicate alignments, request membership, donor batching, warm cache, parameter/content invalidation, cache-only misses, execution-failure fallback, atomic cache rollback, optional exports, and end-to-end decisions. All generated data and outputs use temporary directories.
