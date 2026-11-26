#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * ==============================================================================
 * E-value x addfragments (keeplength on/off) x MaxAlign (on/off) Benchmark
 * ==============================================================================
 * Methods:
 *   1. msa_ml         : IQ-TREE 2 (-m JTT+I+G4 -B 1000 -keep-ident) + MAD + GoTree
 *   2. fastme_poisson : MSA Poisson matrix + FastME (-m N -n -s) + MAD + GoTree
 *   3. msa_nj         : MSA Poisson matrix + RapidNJ (-i pd) + MAD + GoTree
 *   4. gs_poisson     : MSA Poisson matrix + GS2 (-d -s -l) + GoTree
 */

process mafft_align_ori {
    label 'short'
    publishDir "${params.output_dir}/msa_ori", mode: 'copy'

    input:
        tuple val(id), path(fa)

    output:
        tuple val(id), path("${id}.msa.fa"), emit: ori_msa

    script:
        def threads = task.cpus ?: 4
        """
        set -euo pipefail
        mafft --anysymbol --quiet --auto --thread ${threads} --threadit 0 "${fa}" \
        | awk '/^>/{if(x) print x; print; x=""; next} {x=x\$0} END{print x}' \
        > "${id}.msa.fa"
        """
}

process mmseqs_search_eval {
    label 'mmseqs'
    publishDir "${params.output_dir}/homologs/e_${eval_val}", mode: 'copy',
        saveAs: { filename -> filename.endsWith('.homologs.fa') ? filename : null }

    input:
        tuple val(id), path(fa), val(eval_val)

    output:
        tuple val(id), val(eval_val), path("${id}.homologs.fa"), emit: homolog_seqs

    script:
        def db_target   = "${params.database_dir}/${params.database}"
        def tmp_root    = params.tmp_dir ?: './tmp'
        def max_depth   = (params.max_depth ?: 5).toString()
        def random_seed = (params.random_seed ?: 42).toString()
        def threads     = task.cpus ?: 8
        def tmp_abs     = new File(tmp_root.toString()).absolutePath
        """
        set -euo pipefail
        export AMPLIPHY_RANDOM_SEED="${random_seed}"

        # Create dummy taxonomy status (disabled) for fast search
        printf "disabled\tbenchmark\n" > _tax_status.tsv

        AMPLIPHY_RANDOM_SEED="${random_seed}" ampliphy_mmseqs_search.sh \
            --id "${id}" \
            --input "${fa}" \
            --db-target "${db_target}" \
            --tmp-root "${tmp_abs}/e_${eval_val}" \
            --taxonomy-status _tax_status.tsv \
            --mmseqs-options "-e ${eval_val} --split-memory-limit 25G" \
            --max-depth "${max_depth}" \
            --max-seqs "0" \
            --threads ${threads} \
            --db-load-mode "0"
        """
}

process build_variant_msa {
    label 'short'
    publishDir "${params.output_dir}/e_${eval_val}/${variant}/msa", mode: 'copy',
        pattern: "${id}.msa_final.fa",
        saveAs: { "${id}.fa" }

    input:
        tuple val(id), val(eval_val), val(variant), path(ori_msa), path(hom_fa)
        path maxalign_pl

    output:
        tuple val(id), val(eval_val), val(variant), path("${id}.msa_final.fa"), path(hom_fa), emit: variant_msa

    script:
        def threads        = task.cpus ?: 4
        def use_keeplength = variant.contains('keep_on')
        def use_maxalign   = variant.contains('maxa_on')
        def kl_flag        = use_keeplength ? '--keeplength' : ''
        """
        set -euo pipefail

        if [[ ! -s "${hom_fa}" ]]; then
            cp "${ori_msa}" _amp.fa
        else
            mafft --anysymbol --quiet ${kl_flag} --thread ${threads} --threadit 0 --auto \
                --addfragments "${hom_fa}" "${ori_msa}" \
            | awk '/^>/{if(x) print x; print; x=""; next} {x=x\$0} END{print x}' \
            > _amp.fa
        fi

        if [[ "${use_maxalign}" == "true" ]]; then
            awk '/^>/{print substr(\$1,2)}' "${ori_msa}" > _keepids.txt
            awk -v kf="_keepids.txt" '
                BEGIN { while ((getline l < kf) > 0) keep[l] = 1 }
                /^>/ {
                    sid = substr(\$1, 2)
                    if (sid in keep) { print ">+" substr(\$0, 2) }
                    else             { print }
                    next
                }
                { print }
            ' _amp.fa > _preprocessed.fa

            perl "${maxalign_pl}" -p _preprocessed.fa \
            | awk '/^>/{if(x) print x; print; x=""; next} {x=x\$0} END{print x}' \
            > "${id}.msa_final.fa"
        else
            cp _amp.fa "${id}.msa_final.fa"
        fi
        """
}

process compute_poisson_matrix {
    label 'minimal'
    publishDir "${params.output_dir}/e_${eval_val}/${variant}/matrix", mode: 'copy',
        pattern: "${id}.poisson.phy"

    input:
        tuple val(id), val(eval_val), val(variant), path(msa_fa), path(hom_fa)

    output:
        tuple val(id), val(eval_val), val(variant), path("${id}.poisson.phy"), path(hom_fa), emit: matrix_tuple

    script:
        """
        python3 -m src.msa_distance matrix "${msa_fa}" "${id}.poisson.phy"
        """
}

process infer_msa_ml {
    label 'tree'
    publishDir "${params.output_dir}/e_${eval_val}/${variant}/msa_ml/tree_amp", mode: 'copy',
        pattern: "${id}.amp.rooted.nwk", saveAs: { "${id}.amp.nwk" }
    publishDir "${params.output_dir}/e_${eval_val}/${variant}/msa_ml/tree_enr", mode: 'copy',
        pattern: "${id}.enr.nwk"

    input:
        tuple val(id), val(eval_val), val(variant), path(msa_fa), path(hom_fa)
        path mad_script

    output:
        tuple val(id), val(eval_val), val(variant), val('msa_ml'), path("${id}.amp.rooted.nwk"), path("${id}.enr.nwk"), emit: trees

    script:
        def threads = task.cpus ?: 4
        """
        set -euo pipefail
        NSEQ=\$(grep -c '^>' "${msa_fa}" || true)
        if [[ "\${NSEQ}" -lt 4 ]]; then
            iqtree -s "${msa_fa}" -pre "${id}.ml" -T ${threads} -m JTT+I+G4 -keep-ident
        else
            iqtree -s "${msa_fa}" -pre "${id}.ml" -T ${threads} -m JTT+I+G4 -B 1000 -keep-ident
        fi
        cp "${id}.ml.treefile" "${id}.amp.nwk"

        python3 "${mad_script}" "${id}.amp.nwk" -n
        head -n 1 "${id}.amp.nwk.rooted" > "${id}.amp.rooted.nwk"

        if [[ -s "${hom_fa}" ]]; then
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
        python3 -m src.msa_distance clean-tree "${id}.enr.nwk"
        """
}

process infer_fastme_poisson {
    label 'short'
    publishDir "${params.output_dir}/e_${eval_val}/${variant}/fastme_poisson/tree_amp", mode: 'copy',
        pattern: "${id}.amp.rooted.nwk", saveAs: { "${id}.amp.nwk" }
    publishDir "${params.output_dir}/e_${eval_val}/${variant}/fastme_poisson/tree_enr", mode: 'copy',
        pattern: "${id}.enr.nwk"

    input:
        tuple val(id), val(eval_val), val(variant), path(phy_file), path(hom_fa)
        path mad_script

    output:
        tuple val(id), val(eval_val), val(variant), val('fastme_poisson'), path("${id}.amp.rooted.nwk"), path("${id}.enr.nwk"), emit: trees

    script:
        """
        set -euo pipefail
        fastme -m N -n -s -i "${phy_file}" -o "${id}.amp.nwk"
        python3 -m src.msa_distance clean-tree "${id}.amp.nwk"

        python3 "${mad_script}" "${id}.amp.nwk" -n
        head -n 1 "${id}.amp.nwk.rooted" > "${id}.amp.rooted.nwk"

        if [[ -s "${hom_fa}" ]]; then
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
        python3 -m src.msa_distance clean-tree "${id}.enr.nwk"
        """
}

process infer_msa_nj {
    label 'short'
    publishDir "${params.output_dir}/e_${eval_val}/${variant}/msa_nj/tree_amp", mode: 'copy',
        pattern: "${id}.amp.rooted.nwk", saveAs: { "${id}.amp.nwk" }
    publishDir "${params.output_dir}/e_${eval_val}/${variant}/msa_nj/tree_enr", mode: 'copy',
        pattern: "${id}.enr.nwk"

    input:
        tuple val(id), val(eval_val), val(variant), path(phy_file), path(hom_fa)
        path mad_script

    output:
        tuple val(id), val(eval_val), val(variant), val('msa_nj'), path("${id}.amp.rooted.nwk"), path("${id}.enr.nwk"), emit: trees

    script:
        """
        set -euo pipefail
        rapidnj "${phy_file}" -i pd -x "${id}.amp.nwk"
        python3 -m src.msa_distance clean-tree "${id}.amp.nwk"

        python3 "${mad_script}" "${id}.amp.nwk" -n
        head -n 1 "${id}.amp.nwk.rooted" > "${id}.amp.rooted.nwk"

        if [[ -s "${hom_fa}" ]]; then
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
        python3 -m src.msa_distance clean-tree "${id}.enr.nwk"
        """
}

process infer_gs_poisson {
    label 'short'
    publishDir "${params.output_dir}/e_${eval_val}/${variant}/gs_poisson/tree_amp", mode: 'copy',
        pattern: "${id}.amp.nwk"
    publishDir "${params.output_dir}/e_${eval_val}/${variant}/gs_poisson/tree_enr", mode: 'copy',
        pattern: "${id}.enr.nwk"

    input:
        tuple val(id), val(eval_val), val(variant), path(phy_file), path(hom_fa)

    output:
        tuple val(id), val(eval_val), val(variant), val('gs_poisson'), path("${id}.amp.nwk"), path("${id}.enr.nwk"), emit: trees

    script:
        def gs2_bin = params.gs2_bin ?: 'gs2'
        """
        set -euo pipefail
        "${gs2_bin}" -d -s -l "${phy_file}" > "${id}.amp.nwk"
        python3 -m src.msa_distance clean-tree "${id}.amp.nwk"

        # GS first-split root is preserved directly (no MAD)
        if [[ -s "${hom_fa}" ]]; then
            awk '/^>/{print substr(\$1,2)}' "${hom_fa}" > tips.txt
            if command -v gotree &>/dev/null; then
                gotree prune -i "${id}.amp.nwk" -f tips.txt -o "${id}.enr.nwk"
            else
                cp "${id}.amp.nwk" "${id}.enr.nwk"
            fi
            rm -f tips.txt
        else
            cp "${id}.amp.nwk" "${id}.enr.nwk"
        fi
        python3 -m src.msa_distance clean-tree "${id}.enr.nwk"
        """
}

process compute_tcs {
    label 'minimal'

    input:
        tuple val(id), val(eval_val), val(variant), val(method), path(amp_tree), path(enr_tree)
        path lineage_file
        path tcs_script

    output:
        path "${id}.${eval_val}.${variant}.${method}.tcs.tsv", emit: tcs_row

    script:
        """
        TCS_AMP=\$(python3 "${tcs_script}" "${amp_tree}" "${lineage_file}" 2>/dev/null | cut -f2 || echo "0.0")
        TCS_ENR=\$(python3 "${tcs_script}" "${enr_tree}" "${lineage_file}" 2>/dev/null | cut -f2 || echo "0.0")
        TCSR2_AMP=\$(python3 -m src.tcsr2 "${amp_tree}" "${lineage_file}" 2>/dev/null || echo "nan")
        TCSR2_ENR=\$(python3 -m src.tcsr2 "${enr_tree}" "${lineage_file}" 2>/dev/null || echo "nan")
        printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "${id}" "${eval_val}" "${variant}" "${method}" "\${TCS_AMP}" "\${TCS_ENR}" "\${TCSR2_AMP}" "\${TCSR2_ENR}" \
            > "${id}.${eval_val}.${variant}.${method}.tcs.tsv"
        """
}

workflow {
    def eval_list    = (params.evalues ?: '1e-100,1e-10,1e-3,1e-1,1e-0').toString().split(',').collect { it.trim() }.findAll { it }
    def variant_list = (params.variants ?: 'maxa_on_keep_on,maxa_on_keep_off,maxa_off_keep_on,maxa_off_keep_off').toString().split(',').collect { it.trim() }.findAll { it }

    def mad_script      = file(params.mad_script)
    def tcs_script      = file(params.tcs_script)
    def maxalign_script = file(params.maxalign_script)

    channel
        .fromPath("${params.input_dir}/*.{fa,fasta}", checkIfExists: true)
        .map { f ->
            def base = f.getSimpleName().replaceFirst(/\.(fa|fasta)$/, '')
            tuple(base, f)
        }
        .set { fasta_ch }

    // 1. Original MSA once per HOG
    mafft_align_ori(fasta_ch)

    // 2. MMseqs2 search + seeded random sampling per (HOG, E-value)
    fasta_ch
        .combine(channel.fromList(eval_list))
        .set { search_inputs }

    mmseqs_search_eval(search_inputs)

    // 3. Build 4 MSA variants per (HOG, E-value)
    mafft_align_ori.out.ori_msa
        .combine(mmseqs_search_eval.out.homolog_seqs, by: 0)
        .combine(channel.fromList(variant_list))
        .map { id, ori_msa, eval_val, hom_fa, variant ->
            tuple(id, eval_val, variant, ori_msa, hom_fa)
        }
        .set { variant_inputs }

    build_variant_msa(variant_inputs, maxalign_script)

    // 4. Compute Poisson distance matrix once per (HOG, E-value, variant)
    compute_poisson_matrix(build_variant_msa.out.variant_msa)

    // 5. Infer trees with all 4 methods
    infer_msa_ml(build_variant_msa.out.variant_msa, mad_script)
    infer_fastme_poisson(compute_poisson_matrix.out.matrix_tuple, mad_script)
    infer_msa_nj(compute_poisson_matrix.out.matrix_tuple, mad_script)
    infer_gs_poisson(compute_poisson_matrix.out.matrix_tuple)

    // 6. Compute TCS summary
    if (params.lineage_file && file(params.lineage_file).exists()) {
        infer_msa_ml.out.trees
            .mix(
                infer_fastme_poisson.out.trees,
                infer_msa_nj.out.trees,
                infer_gs_poisson.out.trees
            )
            .set { all_trees }

        compute_tcs(all_trees, file(params.lineage_file), tcs_script)

        compute_tcs.out.tcs_row
            .map { f -> f.text.trim() }
            .collectFile(
                name: 'tcs_summary.tsv',
                seed: "family_id\tevalue\tvariant\tmethod\ttcs_amp\ttcs_enr\ttcsr2_amp\ttcsr2_enr",
                newLine: true,
                sort: true,
                storeDir: "${params.output_dir}/report"
            )
    }
}
