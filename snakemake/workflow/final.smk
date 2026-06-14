# Final genome aggregation across all chromosomes.

rule aggregate_final_genome:
    input:
        filtered_fasta=lambda wildcards: [
            f"{OUT}/impuT2T_patch/{wildcards.run_tag}/{SAMPLE}.{q}.{c}.aggregate.patch.filtered.fasta"
            for q, c in QUERY_CHROMOSOME_PAIRS
        ],
    output:
        final_genome=f"{OUT}/final_genome/{{sample}}.aggregate.patch.filtered.{{run_tag}}.fasta",
        final_agp=f"{OUT}/final_genome/{{sample}}.aggregate.patch.filtered.{{run_tag}}.agp",
        final_edge_log=f"{OUT}/final_genome/{{sample}}.aggregate.patch.filtered.{{run_tag}}.edge.log",
    params:
        paths=lambda wildcards: " ".join(
            [
                f"{OUT}/impuT2T_patch/{wildcards.run_tag}/{SAMPLE}.{q}.{c}.aggregate.patch.filtered.fasta"
                for q, c in QUERY_CHROMOSOME_PAIRS
            ]
        ),
        names=lambda wildcards: " ".join(
            [f"{q}_{c}" for q, c in QUERY_CHROMOSOME_PAIRS]
        ),
        agp_paths=lambda wildcards: " ".join(
            [
                f"{OUT}/impuT2T_patch/{wildcards.run_tag}/{SAMPLE}.{q}.{c}.aggregate.patch.filtered.agp"
                for q, c in QUERY_CHROMOSOME_PAIRS
            ]
        ),
        edge_paths=lambda wildcards: " ".join(
            [
                f"{OUT}/impuT2T_patch/{wildcards.run_tag}/{SAMPLE}.{q}.{c}.aggregate.edge.log"
                for q, c in QUERY_CHROMOSOME_PAIRS
            ]
        ),
    benchmark:
        f"{OUT}/benchmarks/{{sample}}_aggregate_final_genome_{{run_tag}}.txt"
    log:
        f"{OUT}/logs/{{sample}}_aggregate_final_genome_{{run_tag}}.txt"
    threads:
        1
    shell:
        """
        mkdir -p $(dirname {output.final_genome})
        python3 scripts/collect_chromosomes_into_one.py \
            -i {params.paths} \
            -n {params.names} \
            -o {output.final_genome} \
            -a {params.agp_paths} \
            -g {output.final_agp}

        cat {params.edge_paths} > {output.final_edge_log}
        """
