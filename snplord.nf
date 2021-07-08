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
include { SNP_DISTS } from "./modules/snp-dists.nf"
include { IQTREE } from "./modules/iqtree.nf"

// Define input channels
Channel.fromFilePairs(params.reads)
    .set {reads_ch}

/* Helper function for SNIPPYCORE
def list_to_string(list) {
    def stringList = list
    String result = {stringList.join(' ')}
    return result
}
*/

// snplord workflow
workflow {
    ref_file = file(params.ref, type: 'file')
    SNIPPY(reads_ch, ref_file)
    SNIPPYCORE(SNIPPY.out.snps.map{ it -> it[1] }.collect(), ref_file)

   if (!params.use_alignment || 'both') {
        SNP_DISTS(SNIPPYCORE.out.core_aln.concat( SNIPPYCORE.out.full_aln ))
        IQTREE(SNIPPYCORE.out.core_aln.concat( SNIPPYCORE.out.full_aln ), ref_file)
        }
   else if (params.use_alignment == 'core') {
       SNP_DISTS(SNIPPYCORE.out.core_aln)
        IQTREE(SNIPPYCORE.out.core_aln, ref_file)
        }

    else if (params.use_alignment == 'full') {
        SNP_DISTS(SNIPPYCORE.out.full_aln)
        IQTREE(SNIPPYCORE.out.full_aln, ref_file)
        }

    else {
        log.error "Alignment must be one of 'core', 'full', or 'both'"
        exit(1)
        }

}


 
