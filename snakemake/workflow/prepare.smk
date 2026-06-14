# Genome preparation: ragtag scaffolding through contig assignment.

rule run_ragtag:
    input:
        query=lambda wildcards: QUERIES_DICT[wildcards.query_name],
        small_pan=lambda wildcards: SMALL_PAN_DICT[wildcards.pan_name],
    output:
        ragtag_paf=f"{OUT}/ragtag_output/{{sample}}_{{query_name}}_{{pan_name}}/ragtag.scaffold.asm.paf",
        ragtag_agp=f"{OUT}/ragtag_output/{{sample}}_{{query_name}}_{{pan_name}}/ragtag.scaffold.agp",
    benchmark:
        f"{OUT}/benchmarks/{{sample}}_{{query_name}}_{{pan_name}}.run_ragtag.txt"
    log:
        f"{OUT}/logs/{{sample}}_{{query_name}}_{{pan_name}}.run_ragtag.txt"
    conda:
        RAGTAG_CONDA_ENV
    threads:
        RAGTAG_THREADS
    shell:
        """
        ragtag.py scaffold \\
            {input.small_pan} \\
            {input.query} \\
            -o $(dirname {output.ragtag_paf}) \\
            -t {threads}
        """


rule index_query:
    input:
        query=lambda wildcards: QUERIES_DICT[wildcards.query_name],
    output:
        index=f"{OUT}/indexes/{{sample}}_{{query_name}}.fai",
    benchmark:
        f"{OUT}/benchmarks/{{sample}}_{{query_name}}_index_query.txt"
    log:
        f"{OUT}/logs/{{sample}}_{{query_name}}_index_query.txt"
    shell:
        """
        samtools faidx {input.query} -o {output.index}
        """


rule generate_consensus_config:
    input:
        agp_files=lambda wildcards: expand(
            f"{OUT}/ragtag_output/{{sample}}_{{query_name}}_{{pan_name}}/ragtag.scaffold.agp",
            sample=[SAMPLE],
            query_name=[wildcards.query_name],
            pan_name=SMALL_PAN_NAMES,
        ),
        fai_file=f"{OUT}/indexes/{{sample}}_{{query_name}}.fai",
    output:
        config_file=f"{OUT}/consensus_scaffold/{{sample}}_{{query_name}}_config.json",
    benchmark:
        f"{OUT}/benchmarks/{{sample}}_{{query_name}}_generate_consensus_config.txt"
    log:
        f"{OUT}/logs/{{sample}}_{{query_name}}_generate_consensus_config.txt"
    params:
        agp_names=lambda wildcards: [
            f"{SAMPLE}.{wildcards.query_name}.{p}" for p in SMALL_PAN_NAMES
        ],
        agp_prefixes=lambda wildcards: [
            SMALL_PAN_PREFIXES_DICT[p] for p in SMALL_PAN_NAMES
        ],
    script:
        os.path.join(SCRIPTS, "generate_consensus_config.py")


rule consensus_scaffold:
    input:
        config_file=f"{OUT}/consensus_scaffold/{{sample}}_{{query_name}}_config.json",
    output:
        consensus=f"{OUT}/consensus_scaffold/{{sample}}_{{query_name}}_consensus.csv",
    params:
        extra_flags=lambda wildcards: " ".join(
            []
            + (["--X"] if "chrX" in QUERY_TO_CHROMOSOMES[wildcards.query_name] else [])
            + (["--Y"] if "chrY" in QUERY_TO_CHROMOSOMES[wildcards.query_name] else [])
            + (["--M"] if "chrM" in QUERY_TO_CHROMOSOMES[wildcards.query_name] else [])
        ),
    benchmark:
        f"{OUT}/benchmarks/{{sample}}_{{query_name}}_consensus_scaffold.txt"
    log:
        f"{OUT}/logs/{{sample}}_{{query_name}}_consensus_scaffold.txt"
    threads:
        THREADS
    shell:
        """
        python3 scripts/consensus_scaffold.py --config {input.config_file} -o {output.consensus} {params.extra_flags}
        """


rule assign_contigs_to_chromosomes:
    input:
        consensus=f"{OUT}/consensus_scaffold/{{sample}}_{{query_name}}_consensus.csv",
        fasta=lambda wildcards: QUERIES_DICT[wildcards.query_name],
    output:
        assigned_dir=directory(f"{OUT}/contigs_assignment/{{sample}}_{{query_name}}_assigned"),
    params:
        output_prefix=f"{OUT}/contigs_assignment/{{sample}}_{{query_name}}_assigned",
    benchmark:
        f"{OUT}/benchmarks/{{sample}}_{{query_name}}_assign_contigs_to_chromosomes.txt"
    log:
        f"{OUT}/logs/{{sample}}_{{query_name}}_assign_contigs_to_chromosomes.txt"
    threads:
        1
    shell:
        """
        python3 scripts/assign_aggrgate_contigs.py \
            --contigs_file {input.fasta} \
            --contigs_csv {input.consensus} \
            --output_path {params.output_prefix}
        """
