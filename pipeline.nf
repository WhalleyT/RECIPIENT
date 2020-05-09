#!/usr/bin/env nextflow

/*
____________________________________________________________________________________________________
STEP 0
Setup
____________________________________________________________________________________________________
*/

def help_args() {
  log.info"""
  ==============================================
  Nextflow pipeline for bacterial vaccine design
  ==============================================

  Usage:

  Typically one would run the pipeline as follows:
  nextflow run pipeline.nf --fasta_folder

  Required arguments:
  --fasta_folder        Path to folder containing assembly files

  Optional arguments:
  --outdir              Output directory to which results will be saved
  --email               Email address to which summary will be sent to
  --name                Job identifier
  """.stripIndent()
}


if(params.help){
  help_args()
  exit 0
}


//grab args
fasta_path="${params.fasta_folder}/*.${params.fasta_extension}"
output_docs = file("$baseDir/docs/output.md")

kraken_db = params.kraken_db

//validate commandline inputs
fasta_dataset = Channel
  .fromPath(fasta_path)
  .ifEmpty{exit 1, "Fasta files not found not found: ${fasta_path}"}
  .map{file -> tuple(file.baseName, file)}

//assign custom run name if need be
custom_runname = params.name

if(!(workflow.runName ==~ /[a-z]+_[a-z]+/)){
  custom_runame = workflow.runName
}


//set up summary output
def summary = [:]
summary['Pipeline Name']  = 'bacterial-vaccine-design'
summary['Run Name']     = custom_runname ?: workflow.runName
summary['Fasta files']  = fasta_path
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile


fasta_dataset.into{
    fasta_dataset_prokka
    fasta_dataset_kraken
}


/*
____________________________________________________________________________________________________
STEP 1
Create annotations and pangenome
____________________________________________________________________________________________________
*/

process prokka {
    container 'file://singularity_recipes/images/prokka.simg'
    
    publishDir "${params.outdir}/prokka", mode: 'copy'
    
    validExitStatus 0, 5
    
    input:
    set fasta_prefix, file(fasta_file) from fasta_dataset_prokka
    
    output:
    file("${fasta_prefix}/${fasta_prefix}.gff") into gff
    file("${fasta_prefix}/${fasta_prefix}.faa") into faa_annots
    //set file("${fasta_prefix}/${fasta_prefix}.err"),
    //    file("${fasta_prefix}/${fasta_prefix}.ffn"),
    //    file("${fasta_prefix}/${fasta_prefix}.fsa"),
    //    file("${fasta_prefix}/${fasta_prefix}.log"),
    //    file("${fasta_prefix}/${fasta_prefix}.tsv"),
    //    file("${fasta_prefix}/${fasta_prefix}.fna"),
    //    file("${fasta_prefix}/${fasta_prefix}.tbl"),
    //    file("${fasta_prefix}/${fasta_prefix}.txt") into prokka
    

    script:
    """
    prokka --outdir ${fasta_prefix} --prefix ${fasta_prefix} ${fasta_file} --centre C  \
    --kingdom $params.kingdom --force  --metagenome --locustag L --norrna --notrna #--gram neg
    """

}


process roary {
    container 'file://singularity_recipes/images/roary.simg'
    label 'multithreaded'
    publishDir "${params.outdir}/roary", mode: 'copy'
    
    validExitStatus "-", 0

    input:
    file gff from gff.collect()

    output:
    file("*") into roary
    file("pan_genome_reference.fa") into pan_genome
    file("pan_genome_sequences/*") into alignment_files
    file("gene_presence_absence.Rtab") into gene_presence
    set file("*accessory*"), 
        file("*.Rtab"),
        file("_*"),
        file("*.txt"),
        file("*csv") into roary_out

    shell:
    '''
    roary -p $(nproc) -e -n -z !{gff}
    '''
}

/*
____________________________________________________________________________________________________
STEP 2
Create k-mer distance measures and QC with Kraken
____________________________________________________________________________________________________
*/

process mash {
    container 'file://singularity_recipes/images/mash.simg'    

    publishDir "${params.outdir}/mash", mode: 'copy'
    label 'multithreaded'
    input:
    file aa_fastas from faa_annots.collect()

    output:
    file("mash_output.phylip") into mash

    script:
    """
    mash triangle -p ${task.cpus}  -a  $aa_fastas > mash_output.phylip
    """
}

process kraken {
    container 'file://singularity_recipes/images/kraken2.simg'

    publishDir "${params.outdir}/kraken", mode: 'copy'

     input:
        set fasta_prefix, file(fasta_file) from fasta_dataset_kraken
    output:
        file("${fasta_prefix}_kraken.report")
    script:
         """
        kraken2 --db ${kraken_db}  --report ${fasta_prefix}_kraken.report ${fasta_file} 
         """
}

/*
____________________________________________________________________________________________________
STEP 3
Find core genes
____________________________________________________________________________________________________
*/

process core_genes {
    container 'file://singularity_recipes/images/python_R.simg'
    publishDir "${params.outdir}/core", mode: 'copy'

    input:
    file presence from gene_presence

    output:
    file("core_genes.txt") into core_gene_set

    script:
    """
    #!/usr/bin/env Rscript
    matrix <- read.table('$presence', header=TRUE, row.names = 1)

    presence_pcnt <- rowSums(matrix) / ncol(matrix) * 100

    core <- matrix[presence_pcnt >= 99,]
    core_genes <- rownames(core)
    write.table(core_genes, "core_genes.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
    """
}

/*
____________________________________________________________________________________________________
STEP 4
Calculate Tajima's D
____________________________________________________________________________________________________
*/



process tajima {
    container 'file://singularity_recipes/images/python_R.simg'
    publishDir "${params.outdir}/tajima", mode: 'copy'
    label 'multithreaded'


    input:
    file tajima_in from alignment_files.collect()

    output:
    file("tajima.txt") into tajima

    script:
    """
    python /opt/tajima.py -F $tajima_in
    """
}


process translate_reference {
    container 'file://singularity_recipes/images/python_R.simg'
    publishDir "${params.outdir}/protein_reference", mode: 'copy'
    
    input:
    file(reference_dna) from pan_genome
    
    output:
    file("protein_reference.faa") into reference_protein
    
    script:
    """
    python /opt/translate.py $reference_dna > protein_reference.faa
    """
}

process split_fasta {
    container 'file://singularity_recipes/images/python_R.simg'
    publishDir "${params.outdir}/loctree_in", mode: 'copy'

    input:
    file(ref) from reference_protein

    output:
    file("*fasta") into single_protein_files

    script:
    """
    python /opt/split_fasta.py $ref
    """
}

