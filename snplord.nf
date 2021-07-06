#!/usr/bin/env nextflow

// Workflow to run Snippy, generate SNP counts and build a tree

/*

Workflow
Step 1. Take a directory of reads and one or more reference genomes

Step 2. Run Snippy

Step 3. Run Snippy-core

Step 4. Get SNP-counts from snippy-core alignment

Step 5. Build tree with IQTree based on snippy-core alignment

*/

// Enable DSL2 for modular workflows
nextflow.enable.dsl=2

// Specify config and module files and processes
include { SNIPPY } from "./modules/snippy.nf"
include { SNIPPYCORE } from "./modules/snippy.nf"

// Define input channels
Channel.fromFilePairs(params.reads)
    .set {reads_ch}

// snplord workflow
workflow {
    ref_file = file(params.ref, type: 'file')
    SNIPPY(reads_ch, ref_file)
    SNIPPYCORE(SNIPPY.out.aln_folders.map{ it -> it[1] }.collect(), ref_file)
}

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