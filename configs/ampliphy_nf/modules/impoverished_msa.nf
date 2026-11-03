nextflow.enable.dsl=2

process impoverished_msa {
    label params.minimal ? 'minimal' : 'short'
    publishDir "${params.output_dir}/msa_imp", mode: 'copy'

    input:
        tuple val(id), path(ori_msa), path(amp_msa)

    output:
        tuple val(id), path("${id}.imp.fa"), emit: imp_msa

    script:
        """
        ampliphy_impoverished_msa.py msa --msa "${amp_msa}" --ori "${ori_msa}" --out "${id}.imp.fa"
        """
}
