/*
module file for pangenome creation
Steps are:
    - Run prokka
    - create pangenome with Roary
*/

process prokka {
    publishDir "${params.outdir}/prokka", mode: 'copy'
    
    validExitStatus "-", 0, 5
    errorStrategy "retry"
    maxRetries 2
    
    input:
    file(fasta)
    
    output:
    path("${fasta.baseName}/${fasta.baseName}.gff", emit: gffs)
    path("${fasta.baseName}/${fasta.baseName}.faa", emit: annotated_faa)
    path("*", emit: prokka_all_out)
    

    script:
    """
    #!/bin/bash
    set -euxo pipefail
    prokka --outdir ${fasta.baseName} --prefix ${fasta.baseName} ${fasta} --centre C  \
    --kingdom $params.kingdom --force  --metagenome --locustag L --norrna --notrna --gram $params.gram
    """

}


process roary {
    label 'multithreaded'
    publishDir "${params.outdir}/pangenome", mode: 'copy'
    
    validExitStatus "-", 0

    input:
    file(gffs)

    output:
    path("*", emit: roary)
    path("pan_genome_reference.fa", emit: pan_genome)
    path("pan_genome_sequences/*", emit: alignment_files)
    path("gene_presence_absence.csv", emit: gene_presence)

    script:
    """
    roary -p ${task.cpus} -e -n -z $gffs
    """
}


process panaroo {
    label 'multithreaded'
    publishDir "${params.outdir}/pangenome", mode: 'copy'

    input:
    file(gffs)

    output:
    path("results/pan_genome_reference.fa", emit pan_genome)
    path("results/gene_presence_absence_roary.csv", emit: gene_presence)
    path("combined_DNA_CDS.fasta", emit: alignment_files)

    script:
    """
    panaroo -t $task.cpus -o results -i *.gff --clean-mode strict
    """
}

