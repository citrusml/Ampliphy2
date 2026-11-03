nextflow.enable.dsl=2

process mafft_amplify {
    label params.minimal ? 'minimal' : 'short'
    publishDir "${params.output_dir}/msa", mode: 'copy'

    input:
        tuple val(id), path(msa), path(homologs)

    output:
        path "${id}.amp.fa"

    script:
        def mafft_preset       = params.mafft_preset.toLowerCase()
        def mafft_options      = params.mafft_options
        def threads            = params.minimal ? 1 : params.threads
        def mafft_amplify_mode = (params.mafft_amplify_mode ?: (params.mafft_mode ?: 'addfragments')).toString().toLowerCase()

        Map mafft_preset_flags = [
            'auto'  : '--auto',
            'fast'  : '--retree 2 --maxiterate 0',
            'linsi' : '--localpair --maxiterate 1000',
            'ginsi' : '--genafpair --maxiterate 1000',
            'einsi' : '--globalpair --maxiterate 1000'
        ]

        if( !mafft_preset_flags.containsKey(mafft_preset) ) {
            error "Invalid MAFFT preset: ${mafft_preset}. Valid options are: ${mafft_preset_flags.keySet().join(', ')}"
        }

        def mafft_flag = (mafft_preset_flags[mafft_preset] + ( mafft_options ? " ${mafft_options}" : "" )).trim()

        """
        set -euo pipefail

        if [[ ! -s "${homologs}" ]]; then
            echo "Warning: No non-identical homologs were retrieved for ${id}; using the original MSA without amplification." >&2
            cp "${msa}" "${id}.amp.fa"
        elif [[ "${mafft_amplify_mode}" == "denovo" || "${mafft_amplify_mode}" == "main" || "${mafft_amplify_mode}" == "scratch" || "${mafft_amplify_mode}" == "full" ]]; then
            # 一からMSAを再構築: 元配列と集めてきたホモログ配列を結合し、de novoアライメントを実行
            sed '/^>/!s/-//g' "${msa}" > combined.fa
            cat "${homologs}" >> combined.fa
            mafft --anysymbol --quiet --thread ${threads} --threadit 0 ${mafft_flag} combined.fa \
            | awk '/^>/{if(x) print x; print; x=""; next} {x=x\$0} END{print x}' \
            > "${id}.amp.fa"
        elif [[ "${mafft_amplify_mode}" == "add" ]]; then
            mafft --anysymbol --quiet --keeplength --thread ${threads} --threadit 0 ${mafft_flag} \
            --add "${homologs}" "${msa}" \
            | awk '/^>/{if(x) print x; print; x=""; next} {x=x\$0} END{print x}' \
            > "${id}.amp.fa"
        else
            # デフォルト: 既存MSAに断片として追加 (--addfragments)
            mafft --anysymbol --quiet --keeplength --thread ${threads} --threadit 0 ${mafft_flag} \
            --addfragments "${homologs}" "${msa}" \
            | awk '/^>/{if(x) print x; print; x=""; next} {x=x\$0} END{print x}' \
            > "${id}.amp.fa"
        fi
        """
}
