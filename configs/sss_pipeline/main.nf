#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * ==============================================================================
 * SSS (Symmetrical Sequence Similarity score) Linear/Log Phylogenetic Pipeline
 * ==============================================================================
 * Inputs:
 * - params.input_dir: Original FASTA directory (e.g., input_fasta/<taxon>)
 * - params.homolog_dir: Homolog FASTA directory (e.g., results/<taxon>/homologs)
 * - params.lineage_file: Lineage TSV for TCS evaluation
 * - params.output_dir: Output base directory (e.g., results_sss/<taxon>/linear)
 * - params.model: Distance transformation ('linear' or 'log')
 * - params.mad_script: Path to mad.py
 * - params.tcs_script: Path to tcs.py
 */

process merge_sequences {
    label 'minimal'

    input:
        tuple val(id), path(ori_fa)
        path homolog_dir

    output:
        tuple val(id), path(ori_fa), path("${id}.amp.fa"), path("${id}.homologs.fa"), emit: amp_fasta

    script:
        """
        HOM_FILE="${homolog_dir}/${id}.homologs.fa"
        if [ -s "\$HOM_FILE" ]; then
            cp "\$HOM_FILE" "${id}.homologs.fa"
            python3 -m src.pwa_pipeline merge "${ori_fa}" "\$HOM_FILE" "${id}.amp.fa"
        else
            touch "${id}.homologs.fa"
            cp "${ori_fa}" "${id}.amp.fa"
        fi
        """
}

process compute_matrix {
    label 'short'
    publishDir "${params.output_dir}/matrix", mode: 'copy'

    input:
        tuple val(id), path(ori_fa), path(amp_fa), path(hom_fa)

    output:
        tuple val(id), path(ori_fa), path(amp_fa), path(hom_fa), path("${id}.amp.phy"), emit: matrix_tuple

    script:
        def model = params.model ?: 'linear'
        def threads = task.cpus ?: 4
        """
        python3 -m src.sss_distance matrix "${amp_fa}" "${id}.amp.phy" \
            --model "${model}" \
            --threads ${threads} \
            --sensitivity 7.5 \
            --evalue 10.0
        """
}

process infer_nj_tree {
    label 'short'
    publishDir "${params.output_dir}/tree_amp", mode: 'copy'

    input:
        tuple val(id), path(ori_fa), path(amp_fa), path(hom_fa), path(amp_phy)

    output:
        tuple val(id), path(ori_fa), path(amp_fa), path(hom_fa), path("${id}.amp.nwk"), emit: tree_tuple

    script:
        """
        if command -v rapidnj &>/dev/null; then
            rapidnj "${amp_phy}" -i pd -x "${id}.amp.nwk"
        elif command -v panjep &>/dev/null; then
            panjep -t ${task.cpus ?: 4} "${amp_phy}" > "${id}.amp.nwk"
        else
            echo "Error: Neither rapidnj nor panjep found" >&2
            exit 1
        fi

        # Clean negative branches
        python3 -m src.sss_distance clean-tree "${id}.amp.nwk"
        """
}

process root_and_prune {
    label 'short'
    publishDir path: "${params.output_dir}/tree_rooted", mode: 'copy', pattern: '*.amp.rooted.nwk'
    publishDir path: "${params.output_dir}/tree_enr", mode: 'copy', pattern: '*.enr.nwk'

    input:
        tuple val(id), path(ori_fa), path(amp_fa), path(hom_fa), path(amp_nwk)
        path mad_script

    output:
        tuple val(id), path("${id}.amp.rooted.nwk"), path("${id}.enr.nwk"), emit: pruned_trees

    script:
        """
        # 1. MAD Rooting
        python3 "${mad_script}" "${amp_nwk}" -n
        head -n 1 "${amp_nwk}.rooted" > "${id}.amp.rooted.nwk"

        # 2. GoTree Pruning of Homolog Tips -> enr (enriched / pruned) tree
        if [ -s "${hom_fa}" ]; then
            awk '/^>/{print substr(\$1,2)}' "${hom_fa}" > tips.txt
            if command -v gotree &>/dev/null; then
                gotree prune -i "${id}.amp.rooted.nwk" -f tips.txt -o "${id}.enr.nwk"
            else
                cp "${id}.amp.rooted.nwk" "${id}.enr.nwk"
            fi
            rm -f tips.txt
        else
            cp "${id}.amp.rooted.nwk" "${id}.enr.nwk"
        fi

        python3 -m src.sss_distance clean-tree "${id}.enr.nwk"
        """
}

process compute_tcs {
    label 'minimal'

    input:
        tuple val(id), path(amp_rooted_nwk), path(enr_nwk)
        path lineage_file
        path tcs_script

    output:
        path "${id}.tcs.tsv", emit: tcs_row

    script:
        """
        TCS_AMP=\$(python3 "${tcs_script}" "${amp_rooted_nwk}" "${lineage_file}" 2>/dev/null | cut -f2 || echo "0.0")
        TCS_ENR=\$(python3 "${tcs_script}" "${enr_nwk}" "${lineage_file}" 2>/dev/null | cut -f2 || echo "0.0")
        printf '%s\\t%s\\t%s\\n' "${id}" "\${TCS_AMP}" "\${TCS_ENR}" > "${id}.tcs.tsv"
        """
}

workflow {
    def input_dir = params.input_dir
    def homolog_dir = file(params.homolog_dir)
    def mad_script = file(params.mad_script)
    def tcs_script = file(params.tcs_script)

    channel
        .fromPath("${input_dir}/*.{fa,fasta}", checkIfExists: true)
        .map { f ->
            def name = f.getSimpleName()
            def base = name.replaceFirst(/\.(fa|fasta)$/, '')
            tuple(base, f)
        }
        .set { fasta_ch }

    merge_sequences(fasta_ch, homolog_dir)
    compute_matrix(merge_sequences.out.amp_fasta)
    infer_nj_tree(compute_matrix.out.matrix_tuple)
    root_and_prune(infer_nj_tree.out.tree_tuple, mad_script)

    if( params.lineage_file && file(params.lineage_file).exists() ) {
        compute_tcs(root_and_prune.out.pruned_trees, file(params.lineage_file), tcs_script)

        compute_tcs.out.tcs_row
            .map { f -> f.text.trim() }
            .collectFile(
                name: 'tcs_report.tsv',
                seed: "family_id\ttcs_amp\ttcs_enr",
                newLine: true,
                sort: true,
                storeDir: "${params.output_dir}/report"
            )
    }
}
