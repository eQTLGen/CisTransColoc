#!/bin/bash nextflow

process MakeLoci {
    //scratch true

    publishDir "${params.OutputDir}", mode: 'copy', overwrite: true, pattern: "cis_trans_info.txt"

    input:
        tuple path(sig_res), path(eqtls), path(ref), path(gtf), path(cis_filter), val(lead_variant_win), val(cis_win), val(trans_win), val(p_thresh)

    output:
       tuple path(sig_res), path(eqtls), path(ref), path(gtf), path(cis_filter), val(lead_variant_win), val(cis_win), val(trans_win), val(p_thresh), path("cis_trans_info.txt")

    script:
        """
        MakeLoci.R \
        --sig_res ${sig_res} \
        --eqtl_folder ${eqtls} \
        --reference ${ref} \
        --gtf ${gtf} \
        --lead_variant_win ${lead_variant_win} \
        --cis_win ${cis_win} \
        --trans_win ${trans_win} \
        --p_thresh ${p_thresh} \
        --cis_gene_filter ${cis_filter}
        """
}

process Coloc {
    input:
        tuple path(sig_res), path(eqtls), path(ref), path(gtf), path(cis_filter), val(lead_variant_win), val(cis_win), val(trans_win), val(p_thresh), path(loci), val(gene)

    output:
        path("*_coloc.txt")

    script:
        """
        Hyprcoloc.R \
        --gene_id ${gene} \
        --loci ${loci} \
        --eqtl_folder ${eqtls} \
        --output ${gene}_coloc.txt
        """
}

workflow MAKELOCI {
    take:
        data

    main:
        loci_ch = MakeLoci(data)
        
    emit:
        loci_ch

}

workflow COLOC {
    take:
        data

    main:
        coloc_output_ch = Coloc(data)

    emit:
        coloc_output_ch
}
