nextflow.enable.dsl = 2

include {kmer_distance; kraken} from "../modules/quality_control.nf" params(params)

workflow qc {
    take:
      fasta_files

    main:
      kmer_distance(fasta_files.collect())
      kraken(fasta_files)
}