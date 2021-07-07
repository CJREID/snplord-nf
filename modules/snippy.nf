
process SNIPPY {
    tag "${id}"

    publishDir "${params.outdir}/snippy" , mode: 'copy'

    input:
    tuple val(id), path(reads)
    path (ref)
    
    output:
    tuple val(id), path("${id}"), emit: snps

script:
    """
    snippy --force \
        --cpus ${task.cpus}\
        --outdir ${id} \
        --ref ${ref} \
        --R1 ${reads[0]}\
        --R2 ${reads[1]}
    """
}

process SNIPPYCORE {
    publishDir "${params.outdir}/snippy_core" , mode: 'copy'
    publishDir "${params.outdir}" , mode: 'copy' ,
        saveAs: { filename -> 
            if (filename == "core.aln") "core_aln.fa"
            else if (filename == "core.vcf") "core.vcf"
            else if (filename == "clean.full.aln") "full_aln.fa"
        }
    
    input:
    path (snps)
    path (ref)

    output:
    path "core.*"
    path "core.aln", emit: core_aln
    path "clean.full.aln", emit: full_aln

    script:
    """
    snippy-core \
        --ref ${ref} \
        ${snps}

    snippy-clean_full_aln \
        core.full.aln > clean.full.aln
    """
}
