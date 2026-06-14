# Pangenome alignment, local patching, and per-chromosome aggregate patch.

rule alignment_pangenome:
    input:
        pangenome=lambda wildcards: pangenome_fasta_path(
            wildcards.chromosome, wildcards.pangenome_name
        ),
        assigned_dir=f"{OUT}/contigs_assignment/{{sample}}_{{query_name}}_assigned",
    output:
        alignment=f"{OUT}/pangenome_alignment_{{sample}}_{{query_name}}/{{chromosome}}/{{sample}}_{{query_name}}_to_{{pangenome_name}}_{{chromosome}}.paf",
    benchmark:
        f"{OUT}/benchmarks/{{sample}}_{{query_name}}_to_{{pangenome_name}}_{{chromosome}}_alignment_pangenome.txt"
    log:
        f"{OUT}/logs/{{sample}}_{{query_name}}_to_{{pangenome_name}}_{{chromosome}}_alignment_pangenome.txt"
    threads:
        1
    shell:
        """
        mkdir -p $(dirname {output.alignment})
        minimap2 -x asm5 -t {threads} {input.pangenome} {input.assigned_dir}/{wildcards.chromosome}.assigned.fa -o {output.alignment}
        """


rule impuT2T_local:
    input:
        assigned_dir=f"{OUT}/contigs_assignment/{{sample}}_{{query_name}}_assigned",
        alignment=f"{OUT}/pangenome_alignment_{{sample}}_{{query_name}}/{{chromosome}}/{{sample}}_{{query_name}}_to_{{pangenome_name}}_{{chromosome}}.paf",
        pangenome=lambda wildcards: pangenome_fasta_path(
            wildcards.chromosome, wildcards.pangenome_name
        ),
    output:
        local_patched=f"{OUT}/local_alignment_{{sample}}_{{query_name}}/{{chromosome}}/{{sample}}.{{query_name}}.{{pangenome_name}}.agp",
    params:
        output_prefix=f"{OUT}/local_alignment_{{sample}}_{{query_name}}/{{chromosome}}/{{sample}}.{{query_name}}.{{pangenome_name}}",
    benchmark:
        f"{OUT}/benchmarks/{{sample}}_{{query_name}}_{{chromosome}}_{{pangenome_name}}_local.txt"
    log:
        f"{OUT}/logs/{{sample}}_{{query_name}}_{{chromosome}}_{{pangenome_name}}_local.txt"
    threads:
        1
    shell:
        """
        mkdir -p $(dirname {output.local_patched})
        python3 scripts/from_paf_to_multi_connections.py \
            -fl {input.alignment} \
            -q {input.assigned_dir}/{wildcards.chromosome}.assigned.fa \
            -o {params.output_prefix} \
            -r {input.pangenome} \
            --info_only
        """


rule tiebreak_manifest:
    """
    One TSV per (sample, query, chromosome) for aggregate_paths --tiebreak_manifest:
    donor_id \\t pangenome PAF \\t assigned query FA \\t pangenome chr FASTA
    """
    input:
        sample_list=SAMPLE_LIST_FILE,
        assigned_fa=f"{OUT}/contigs_assignment/{{sample}}_{{query_name}}_assigned/{{chromosome}}.assigned.fa",
    output:
        f"{OUT}/local_alignment_{{sample}}_{{query_name}}/{{chromosome}}/tiebreak_manifest.tsv",
    params:
        pg_path=lambda wildcards: PANGENOME_PATH,
        out=OUT,
    benchmark:
        f"{OUT}/benchmarks/{{sample}}_{{query_name}}_{{chromosome}}_tiebreak_manifest.txt"
    log:
        f"{OUT}/logs/{{sample}}_{{query_name}}_{{chromosome}}_tiebreak_manifest.txt"
    threads:
        1
    shell:
        """
        mkdir -p "$(dirname "{output}")"
        : > "{output}"
        while IFS= read -r line || [ -n "$line" ]; do
            donor="$(echo "$line" | sed 's/\\r$//' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
            [ -z "$donor" ] && continue
            case "$donor" in '#'*) continue ;; esac
            ref_fa="{params.pg_path}/{wildcards.chromosome}/$donor.fasta"
            [ -f "$ref_fa" ] || continue
            printf '%s\\t{params.out}/pangenome_alignment_{wildcards.sample}_{wildcards.query_name}/{wildcards.chromosome}/{wildcards.sample}_{wildcards.query_name}_to_%s_{wildcards.chromosome}.paf\\t{params.out}/contigs_assignment/{wildcards.sample}_{wildcards.query_name}_assigned/{wildcards.chromosome}.assigned.fa\\t%s\\n' \
                "$donor" "$donor" "$ref_fa" >> "{output}"
        done < "{input.sample_list}"
        """


rule impuT2T_patch:
    input:
        assigned_dir=f"{OUT}/contigs_assignment/{{sample}}_{{query_name}}_assigned",
        local_patched=lambda wildcards: expand(
            f"{OUT}/local_alignment_{{sample}}_{{query_name}}/{{chromosome}}/{{sample}}.{{query_name}}.{{pangenome_name}}.agp",
            sample=[wildcards.sample],
            query_name=[wildcards.query_name],
            chromosome=[wildcards.chromosome],
            pangenome_name=pangenome_names_for(wildcards.chromosome),
        ),
        tiebreak_manifest=f"{OUT}/local_alignment_{{sample}}_{{query_name}}/{{chromosome}}/tiebreak_manifest.tsv",
    output:
        patch_fasta=f"{OUT}/impuT2T_patch/{{run_tag}}/{{sample}}.{{query_name}}.{{chromosome}}.aggregate.patch.fasta",
        patch_agp=f"{OUT}/impuT2T_patch/{{run_tag}}/{{sample}}.{{query_name}}.{{chromosome}}.aggregate.agp",
        filtered_fasta=f"{OUT}/impuT2T_patch/{{run_tag}}/{{sample}}.{{query_name}}.{{chromosome}}.aggregate.patch.filtered.fasta",
        releveant_seq=temp(
            f"{OUT}/impuT2T_patch/{{run_tag}}/{{sample}}.{{query_name}}.{{chromosome}}.aggregate.relevant_seq.fasta"
        ),
        releveant_seq_fai=temp(
            f"{OUT}/impuT2T_patch/{{run_tag}}/{{sample}}.{{query_name}}.{{chromosome}}.aggregate.relevant_seq.fasta.fai"
        ),
    params:
        output_prefix=lambda wildcards: (
            f"{OUT}/impuT2T_patch/{wildcards.run_tag}/"
            f"{wildcards.sample}.{wildcards.query_name}.{wildcards.chromosome}.aggregate"
        ),
        local_prefix=lambda wildcards: (
            f"{OUT}/local_alignment_{wildcards.sample}_{wildcards.query_name}/"
            f"{wildcards.chromosome}/{wildcards.sample}.{wildcards.query_name}."
        ),
        extra_flags=lambda wildcards: " ".join(
            []
            + (["--sample_num 350"] if wildcards.chromosome == "chrX" else [])
            + (["--sample_num 119"] if wildcards.chromosome == "chrY" else [])
        ),
    benchmark:
        f"{OUT}/benchmarks/{{sample}}_{{query_name}}_{{chromosome}}_{{run_tag}}_patch.txt"
    log:
        f"{OUT}/logs/{{sample}}_{{query_name}}_{{chromosome}}_{{run_tag}}_patch.txt"
    threads:
        THREADS
    shell:
        """
        mkdir -p $(dirname {output.patch_fasta})
        python3 scripts/aggregate_paths_info_semi-greedy_greedy-tiebreak_0517.py \
            -p {params.local_prefix} \
            -l {SAMPLE_LIST_FILE} \
            -o {params.output_prefix} \
            --contig_fasta {input.assigned_dir}/{wildcards.chromosome}.assigned.fa \
            --ref_path {PANGENOME_PATH}/{wildcards.chromosome} \
            {params.extra_flags} \
            --tiebreak \
            --tiebreak_manifest {input.tiebreak_manifest} \
            --tiebreak_jobs {threads} \
            --tiebreak_threads 1 \
            --timeout 100 --retry_timeout 100  --no_optional_tie_skips  \
            --second_pass --weight_ratio 0.5 >> {log} 2>&1

        python3 scripts/extract_fasta_patch_ratio.0402.py \
            -i {output.patch_agp} \
            -f {output.patch_fasta} \
            --edge-info {params.output_prefix}.edge.log \
            --disable-endpoint-segment-flank-filter \
            -o {params.output_prefix}.patch >> {log} 2>&1
        """
