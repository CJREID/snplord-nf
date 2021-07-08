
process SNP_DISTS{
    echo true

    publishDir "${params.outdir}/snp-dists" , mode: 'copy'

    input:
        path (input)

    output:
        path ("*snp-dists.csv")

    script:
    """
    snp-dists -c -b ${input} > ${input}.snp-dists.csv
    """
}