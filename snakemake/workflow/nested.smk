# Nested meta-workflow: one Snakemake child invocation per (sample, query, chromosome).

rule run_chromosome:
    input:
        assigned_dir=f"{OUT}/contigs_assignment/{{sample}}_{{query_name}}_assigned",
    output:
        done=f"{OUT}/checkpoints/{{sample}}_{{query_name}}_{{chromosome}}.done",
    params:
        target=lambda wildcards: patch_filtered_fasta(
            wildcards.sample,
            wildcards.query_name,
            wildcards.chromosome,
            RUN_TAG,
        ),
        configfile=ACTIVE_CONFIGFILE,
    benchmark:
        f"{OUT}/benchmarks/{{sample}}_{{query_name}}_{{chromosome}}_run_chromosome.txt"
    log:
        f"{OUT}/logs/{{sample}}_{{query_name}}_{{chromosome}}_run_chromosome.txt"
    threads:
        CHR_THREADS
    shell:
        """
        mkdir -p $(dirname {output.done})
        snakemake -s workflow/chromosome.smk \
            --configfile {params.configfile} \
            --cores {threads} \
            --rerun-incomplete \
            --keep-going \
            --latency-wait 60 \
            --config run_tag={RUN_TAG} \
                     sample_list_file={SAMPLE_LIST_FILE} \
                     query_name={wildcards.query_name} \
                     chromosome={wildcards.chromosome} \
            -- {params.target}
        touch {output.done}
        """
