nextflow.enable.dsl = 2

include {prokka; roary; panaroo} from "../modules/pangenome.nf" params(params)

workflow pangenome {
    take:
      fasta_files
    
    main:
      prokka(fasta_files)
      
      if (params.pangenome == "roary") {
      roary(prokka.out.gffs.collect())
      } else if (params.pangenome == "panaroo"){
        panaroo(prokka.out.gffs.collect())
      }else if (params.pangenome == "pirate"){
        pirate(prokka.out.gffs.collect())
      }
    
    emit:
      prokka.out.prokka_all_out
}

