#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2


def helpmessage() {

log.info"""

CisTransColocalisation v${workflow.manifest.version}"
===========================================================
Pipeline for running HyprColoc colocalisation analyses (https://www.nature.com/articles/s41467-020-20885-8) between cis- and trans-eQTL loci in eQTLGen eQTL summary statistics.

Usage:

nextflow run main.nf 
--eqtl_files \
--allele_info \
--gtf \
--OutputDir

Mandatory arguments:
--eqtl_files                eQTLGen parquet dataset.
--sig_eqtls                 Text file with significant eQTLGen eQTLs.
--allele_info               Parquet file with alleles and SNP positions for eQTL dataset.
--gtf                       GTF file for gene annotation.

Optional arguments:
--OutputDir                 Output directory. Defaults to "results".
--leadvar_window            Window used in distance pruning to find locus lead variants. Defaults to 1000000.
--cis_window                cis-eQTL window. Defaults to 1000000.
--trans_window              trans-eQTL window. Defaults to 5000000.
--p_thresh                  P-value threshold for significant effects. Defaults to 5e-8.
--i2_thresh                 Heterogeneity threshold. Defaults to 40 (<40%).
--maxN_thresh               Per gene mimaximal sample size threshold. Defaults to 0.8 (SNPs with >=0.8*max(N))
--minN_thresh               Minimal sample size threshold. Defaults to 0 (no filtering)
--cis_gene_filter           Optional filter for the cis-eQTL genes. By default, all genes are queried.
""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

// Default parameters
params.OutputDir = 'results'
params.leadvar_window = 1000000
params.cis_window = 1000000
params.trans_window = 5000000
params.p_thresh = 5e-8
params.i2_thresh = 40
params.maxN_thresh = 0.8
params.minN_thresh = 0
params.cis_gene_filter = 'data/help_input.txt'


//Show parameter values
log.info """=======================================================
eQTLGen cis-trans colocalisation pipeline v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Version']                         = workflow.manifest.version
summary['Current user']                             = "$USER"
summary['Current home']                             = "$HOME"
summary['Current path']                             = "$PWD"
summary['Working dir']                              = workflow.workDir
summary['Script dir']                               = workflow.projectDir
summary['Config Profile']                           = workflow.profile
summary['Container Engine']                         = workflow.containerEngine
if(workflow.containerEngine) summary['Container']   = workflow.container
summary['Output directory']                         = params.OutputDir
summary['eQTL results']                             = params.eqtl_files
summary['Allele info file']                         = params.allele_info
summary['GTF file']                                 = params.gtf
summary['Pruning window']                           = params.leadvar_window
summary['cis-eQTL window']                          = params.cis_window
summary['trans-eQTL window']                        = params.trans_window
summary['P threshold']                              = params.p_thresh
summary['I2 threshold']                              = params.i2_thresh
summary['maxN threshold']                            = params.maxN_thresh
summary['minN threshold']                            = params.minN_thresh

// import modules
include { MAKELOCI; COLOC; MakeLoci; Coloc } from './modules/CisTransColocalization.nf'

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

// Define argument channels
// Get eQTL channel
empirical_ch = Channel.fromPath(params.eqtl_files, type: 'dir')
sig_ch = Channel.fromPath(params.sig_eqtls)
allele_ch = Channel.fromPath(params.allele_info)
gtf_ch = Channel.fromPath(params.gtf)
cis_filter_ch = Channel.fromPath(params.cis_gene_filter)

input_ch = sig_ch.combine(empirical_ch).combine(allele_ch).combine(gtf_ch).combine(cis_filter_ch)

// Define parameter channels
leadvar_window = Channel.value(params.leadvar_window)
cis_window = Channel.value(params.cis_window)
trans_window = Channel.value(params.trans_window)
p_thresh = Channel.value(params.p_thresh)
i2_thresh = Channel.value(params.i2_thresh)
maxN_thresh = Channel.value(params.maxN_thresh)
minN_thresh = Channel.value(params.minN_thresh)

input_ch = input_ch.combine(leadvar_window).combine(cis_window).combine(trans_window)
.combine(p_thresh).combine(i2_thresh).combine(maxN_thresh).combine(minN_thresh)


workflow {
        MAKELOCI(input_ch)
        
        genes = MAKELOCI.out.map { tuple -> tuple[-1] }
        
        genes = genes
           .splitCsv(header: true, sep: ',')
           .map { row ->
           def key = "${row.cis_gene}"
           return tuple(key)
           }

        coloc_input_ch = MAKELOCI.out.combine(genes)
        
        COLOC(coloc_input_ch)
        COLOC.out.flatten().collectFile(name: 'CisTransColocResults.txt', keepHeader: true, sort: true, storeDir: "${params.OutputDir}")
        }


workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
