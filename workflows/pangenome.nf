nextflow.enable.dsl = 2

include {prokka; roary; panaroo} from "../modules/pangenome.nf" params(params)

workflow pangenome {
    take:
      fasta_files
    
    main:
      prokka(fasta_files)
      
      if (params.pangenome == "roary") {
      roary(prokka.out.gffs.collect())
      pan_genome_reference = roary.out.pan_genome
      } else if (params.pangenome == "panaroo"){
        panaroo(prokka.out.gffs.collect())
        pan_genome_reference = panaroo.out.pan_genome
      }else if (params.pangenome == "pirate"){
        pirate(prokka.out.gffs.collect())
        pan_genome_reference = piirate.out.pan_genome
      }
    
    emit:
      prokka.out.prokka_all_out
      pan_genome_reference

}

