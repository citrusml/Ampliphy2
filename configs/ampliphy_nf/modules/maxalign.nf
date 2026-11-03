nextflow.enable.dsl=2

process maxalign {
    label params.minimal ? 'minimal' : 'short'
    publishDir "${params.output_dir}/msa", mode: 'copy'

    input:
        tuple val(id), path(amp_fa), path(ref_fa)
        path maxalign_pl

    output:
        path "${id}.maxaligned.fa"

    script:
        def options = params.maxalign_options ?: ''
        """
        set -euo pipefail

        # Extract original sequence IDs from ref_fa and mark them as must-keep (+prefix)
        awk '/^>/{print substr(\$1,2)}' "${ref_fa}" > _keepids.txt
        awk -v kf="_keepids.txt" '
            BEGIN { while ((getline l < kf) > 0) keep[l] = 1 }
            /^>/ {
                id = substr(\$1, 2)
                if (id in keep) { print ">+" substr(\$0, 2) }
                else            { print }
                next
            }
            { print }
        ' "${amp_fa}" > _preprocessed.fa

        perl "${maxalign_pl}" -p ${options} _preprocessed.fa \
        | awk '/^>/{if(x) print x; print; x=""; next} {x=x\$0} END{print x}' \
        > "${id}.maxaligned.fa"
        """
}
