/*
module file for qc components at start of pipeline.
These include:
    - K-mer distance sketch using MASH
    - Species indentification with Kraken2
*/

process kmer_distance {
    publishDir "${params.outdir}/mash", mode: 'copy', overwrite: 'true'

    label "multithreaded"

    input:
      file(fastas)

    output:
      path 'mash_output.phylip', emit: mash

    script:
    """
    mash triangle -p ${task.cpus} $fastas > mash_output.phylip
    """

}

process kraken {
    publishDir "${params.outdir}/kraken", mode: 'copy', overwrite: 'true'

     input:
        file(fasta)
    output:
        file("${fasta.baseName}_kraken.report")
    script:
         """
        kraken2 --db ${params.kraken_db}  --report ${fasta.baseName}_kraken.report ${fasta} 
         """
}