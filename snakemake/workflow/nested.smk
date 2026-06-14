# Nested meta-workflow: one Snakemake child invocation per (sample, query, chromosome).

rule run_chromosome:
    input:
        assigned_dir=f"{OUT}/contigs_assignment/{{sample}}_{{query_name}}_assigned",
    output:
        # Must declare patch outputs here so the master DAG can wire
        # aggregate_final_genome → run_chromosome (impuT2T_patch lives in chromosome.smk).
        patch_fasta=f"{OUT}/impuT2T_patch/{RUN_TAG}/{{sample}}.{{query_name}}.{{chromosome}}.aggregate.patch.filtered.fasta",
        done=f"{OUT}/checkpoints/{{sample}}_{{query_name}}_{{chromosome}}.done",
    params:
        target=lambda wildcards: patch_filtered_fasta(
            wildcards.sample,
            wildcards.query_name,
            wildcards.chromosome,
            RUN_TAG,
        ),
        configfile=os.path.join(WORKFLOW_ROOT, ACTIVE_CONFIGFILE),
        chromosome_smk=CHROMOSOME_SMK,
        workflow_root=WORKFLOW_ROOT,
    benchmark:
        f"{OUT}/benchmarks/{{sample}}_{{query_name}}_{{chromosome}}_run_chromosome.txt"
    log:
        f"{OUT}/logs/{{sample}}_{{query_name}}_{{chromosome}}_run_chromosome.txt"
    threads:
        CHR_THREADS
    # Limit concurrent nested Snakemake invocations if .snakemake metadata races occur:
    #   snakemake --resources nested_workflow=2 ...
    resources:
        nested_workflow=1
    shell:
        """
        mkdir -p $(dirname {output.done})
        cd "{params.workflow_root}" || exit 1
        snakemake -s "{params.chromosome_smk}" \
            --configfile "{params.configfile}" \
            --cores {threads} \
            --nolock \
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
