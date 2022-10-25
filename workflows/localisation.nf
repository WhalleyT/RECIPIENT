nextflow.enable.dsl = 2

include {psortb} from "../modules/localisation.nf" params(params)

workflow localisation {
    take:
      aa_reference

    main:
      psortb(aa_reference)
    
    emit:
    psortb.out.localisation
}