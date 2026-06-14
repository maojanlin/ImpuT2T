# Per-chromosome child workflow (invoked by workflow/nested.smk).
# Requires --configfile plus --config run_tag, sample_list_file, query_name, chromosome.

include: "common.smk"
include: "pangenome.smk"


rule all:
    input:
        patch_outputs(),
