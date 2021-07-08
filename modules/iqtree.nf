
process IQTREE {
    echo true
    
    publishDir "${params.outdir}/iqtree" , mode: 'copy'

    input:
        path (alignment)
        path (ref)

    output:
        path ("*.treefile")
    
        
    shell:
    """
    iqtree -s ${alignment} -pre ${ref.baseName}.${alignment.baseName} -m MFP -bb 1000 -nt AUTO
    """
    }
    


