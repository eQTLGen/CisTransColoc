#!/bin/bash nextflow

process MakeLoci {
    //scratch true

    publishDir "${params.OutputDir}", mode: 'copy', overwrite: true, pattern: "cis_trans_info.txt"

    input:
        tuple path(sig_res), path(eqtls), path(ref), path(gtf), path(cis_filter), val(lead_variant_win), val(cis_win), val(trans_win), val(p_thresh), val(i2), val(maxN), val(minN), val(plot), val(WriteRegion)

    output:
       tuple path(sig_res), path(eqtls), path(ref), path(gtf), path(cis_filter), val(lead_variant_win), val(cis_win), val(trans_win), val(p_thresh), val(i2), val(maxN), val(minN), val(plot), val(WriteRegion), path("cis_trans_info.txt")

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
        --i2_thresh ${i2} \
        --maxN_thresh ${maxN} \
        --minN_thresh ${minN} \
        --cis_gene_filter ${cis_filter}
        """
}

process Coloc {
    tag "${gene}"
    input:
        tuple path(sig_res), path(eqtls), path(ref), path(gtf), path(cis_filter), val(lead_variant_win), val(cis_win), val(trans_win), val(p_thresh), val(i2), val(maxN), val(minN), val(plot), val(WriteRegion), path(loci), val(gene)

    output:
        path("ENSG*")

    if (params.plot) {
        publishDir "${params.OutputDir}", mode: 'copy', overwrite: true, pattern: "*.png"
    } else {}
    if (params.WriteRegionsOut) {
        publishDir "${params.OutputDir}", mode: 'copy', overwrite: true, pattern: "*region.Rds"
    } else {}

    script:
        """
        Hyprcoloc.R \
        --gene_id ${gene} \
        --i2_thresh ${i2} \
        --reference ${ref} \
        --maxN_thresh ${maxN} \
        --minN_thresh ${minN} \
        --loci ${loci} \
        --eqtl_folder ${eqtls} \
        --gtf ${gtf} \
        --output ${gene}_coloc.txt \
        --locusplot ${plot} \
        --WriteRegionOut ${WriteRegion}
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
