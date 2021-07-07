
process IQTREE {
    echo true
    
    publishDir "${params.outdir}/iqtree" , mode: 'copy'

    input:
    path (alignment)

    output:
    path ("*.treefile")
        
    shell:
    """
    iqtree -s ${alignment} -m MFP -bb 1000 -nt AUTO
    """
    }
    


