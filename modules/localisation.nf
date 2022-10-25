/*
module file for localisation prediction.
Can be Psortb
*/


process psortb {
    publishDir "${params.outdir}/localisation", mode: 'copy'

    input:
    file(pan_genome_reference)

    output:
    path("results/*txt", emit: localisation)

    script:
    if (params.gram == "neg")
    """
    mkdir results
    psortb -i $pan_genome_reference -r results --negative
    """
    else
    """
    mkdir results
    psortb -i $pan_genome_reference -r results --positive
    """
}