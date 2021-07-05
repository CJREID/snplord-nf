// Workflow to run Snippy, generate SNP counts and build a tree

/*

Workflow
Step 1. Take a directory of reads and one or more reference genomes

Step 2. Run Snippy

Step 3. Run Snippy-core

Step 4. Get SNP-counts from snippy-core alignment

Step 5. Build tree with IQTree based on snippy-core alignment

*/

#!/usr/bin/env nextflow

process SNIPPY {
    echo true

    tag { $pair_id }

    publishDir "${outdir}/snippy/", mode: 'copy'

    input:
        tuple val(pair_id), file(reads)
        path (ref)
        
    output:
        tuple val(pair_id), path("${pair_id}"), emit: aln_folders
    
    shell:
        """
        snippy --force \
        --cpus ${cpus}\
        --outdir ${pair_id} \
        --ref ${ref} \
        --R1 ${reads[0]}\
        --R2 ${reads[1]}
        """
}

process SNIPPYCORE {
    echo true

    publishDir "${outdir}/snippy_core" , mode: 'copy'
    
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
/*
process IQTREE {
    echo true
    
    publishDir "${outdir}/iqtree" , mode: 'copy'

    input:
    

    output:

    shell:
    """
    iqtree -s core -m MFP -bb 1000 -nt AUTO
    """
}
*/

// Print parameters to screen
params.read_dir = "/projects/AusGEM/Cam/Dairy_Effluent/reads/all_reads/test"
params.output_dir = "/projects/AusGEM/Cam/Dairy_Effluent/reads/all_reads/test/SNPlord_out"
params.ref = "$params.read_dir/GCF_000005845.2_ASM584v2_genomic.fa"
params.reads = "$params.read_dir/*.R{1,2}.fastq.gz"
params.cpu = '20'
params.alignment_tree = ["core", "full", "both"] // default 'both'

//Params for snippy process
outdir = params.output_dir
reads = params.reads
ref = params.ref
cpus = params.cpu

Channel.fromFilePairs(reads)
    .set {reads_ch}
Channel.value(ref)
    .set{ref_ch}

workflow{
    SNIPPY(reads_ch, ref_ch)
    SNIPPYCORE(SNIPPY.out.aln_folders.map{ it -> it[1] }.collect(), ref)
    
    /*
    if params.alignment_tree = "core"
        Channel.fromPath(core_aln)
        .set { iqtree_ch }

    else if params.alignment_tree = "full"
        Channel.fromPath(clean_full_aln)
        .set { iqtree_ch }
        
    else if params.alignment_tree = "both"
        Channel.fromPath(core_aln, clean_full_aln)
        .set { iqtree_ch }
    */

}
