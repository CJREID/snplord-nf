
process SNIPPY {
    echo true

    tag { $pair_id }

    publishDir "${params.outdir}/${ref.baseName}/snippy/", mode: 'copy'

    input:
        tuple val(pair_id), file(reads)
        path (ref)
        
    output:
        tuple val(pair_id), path("${pair_id}"), emit: aln_folders
    
    shell:
        """
        snippy --force \
        --cpus ${task.cpus}\
        --outdir "${pair_id}" \
        --ref ${ref} \
        --R1 ${reads[0]}\
        --R2 ${reads[1]}
        """
}

process SNIPPYCORE {
    echo true

    publishDir "${params.outdir}/${ref.baseName}/snippy_core/" , mode: 'copy',
        saveAs: { filename -> 
            if (filename == "core.aln") "core_aln.fa"
            else if (filename == "core.vcf") "core.vcf"
            else if (filename == "clean.full.aln") "full_aln.fa"
        }

    input: 
        path (aln_folders)
        path (ref)

    output:
    path "core.*"
    path "core.aln", emit: core_aln
    path "clean.full.aln", emit: clean_full_aln

    shell:
    """
    snippy-core \
        --ref ${ref} \
        ${aln_folders[0]}
    snippy-clean_full_aln \
        core.full.aln > clean.full.aln
    """
}


