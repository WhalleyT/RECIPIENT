nextflow.enable.dsl=2

include {qc} from "./workflows/quality_control.nf"
include {pangenome} from "./workflows/pangenome.nf"

workflow {

fasta_dataset = Channel
  .fromPath(params.input)
  .ifEmpty{exit 1, "Fasta files not found not found: ${params.input}"}

  

  main:
    qc(fasta_dataset)
    pangenome(fasta_dataset)

}