nextflow.enable.dsl=2

include {qc} from "./workflows/quality_control.nf"
include {pangenome} from "./workflows/pangenome.nf"



workflow {

fasta_dataset = Channel
  .fromPath(params.input)
  .ifEmpty{exit 1, "Fasta files not found not found: ${params.input}"}

def valid_params = [
    pangenome: ['roary', 'panaroo', 'pirate'],
    gram: ['pos', 'neg']

]

if(!valid_params["pangenome"].contains(params.pangenome)){
  exit 1, 'Invalid pangenome tool, must be roary, panaroo or pirate' 
}

if(!valid_params["gram"].contains(params.gram)){
  exit 1, 'Invalid Gram stain parameter for Prokka, must be pos or neg'
}

  main:
    qc(fasta_dataset)
    pangenome(fasta_dataset)

}