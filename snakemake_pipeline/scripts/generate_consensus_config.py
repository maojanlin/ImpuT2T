#!/usr/bin/env python3
"""
Snakemake script to generate JSON config for consensus_scaffold.py
This script is called by Snakemake and has access to snakemake object
"""
import json
import os

# Get inputs and params from Snakemake
agp_files = snakemake.input.agp_files
fai_file = str(snakemake.input.fai_file)  # FAI file for this specific query
config_file = str(snakemake.output.config_file)
agp_names = snakemake.params.agp_names
agp_prefixes = snakemake.params.agp_prefixes

# Create config dictionary
config = {
    "agp_files": agp_files,
    "agp_names": agp_names,
    "agp_prefixes": agp_prefixes,
    "fai_file": fai_file
}

# Create output directory
os.makedirs(os.path.dirname(config_file), exist_ok=True)

# Write JSON config file
with open(config_file, 'w') as f:
    json.dump(config, f, indent=2)

print(f"Generated config file: {config_file}")
print(f"  - AGP files: {len(agp_files)}")
print(f"  - FAI file: {fai_file}")
print(f"  - Query: {snakemake.wildcards.query_name}")

