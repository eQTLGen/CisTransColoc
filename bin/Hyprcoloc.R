#!/usr/bin/env Rscript

library(hyprcoloc)
library(argparse)
library(data.table)
library(arrow)
library(stringr)
library(IGUtilityPackage)
library(ggplot2)
library(tidyr)
library(rtracklayer)

parser <- ArgumentParser(description = 'Run HyprColoc for every cis-eQTL locus to detect colocalisation between cis- and trans-eQTL.')

parser$add_argument('--loci', metavar = 'file', type = 'character', 
help = 'R data file which contains cis genes, trans genes and variants.')
parser$add_argument('--gene_id', type = 'character', 
help = 'ENSEMBL Gene ID for which to run the analysis.')
parser$add_argument('--eqtl_folder', metavar = 'file', type = 'character', 
help = 'eQTLGen parquet folder format with per gene output files.')
parser$add_argument('--gtf', metavar = 'file', type = 'character',
                    help = "ENSEMBL .gtf file, needs to be hg38.")
parser$add_argument('--i2_thresh', type = 'numeric', default = 40,
                    help = 'Heterogeneity threshold. Defaults to <40%.')
parser$add_argument('--maxN_thresh', type = 'numeric', default = 0.8,
                    help = 'Per gene maximal sample size threshold. Defaults to 0.8 (SNPs with >=0.8*max(N))')
parser$add_argument('--minN_thresh', type = 'numeric', default = 0,
                    help = 'Minimal sample size threshold. Defaults to 0 (no filtering)')
parser$add_argument('--output', metavar = 'file', type = 'character', 
help = 'Output file.')
#parser$add_argument('--locusplot', type = 'logical', action = 'store_false')

args <- parser$parse_args()

# Functions
ParseInput <- function(inp_folder, loci, gene){
  
  message("Reading locus info...")
  
  loci <- fread(args$loci, sep2 = "|")
  setkey(loci, cis_gene)
  message("Reading locus info...done!")
  
  loci <- loci[cis_gene == gene]
  
  # Parse trans-eGenes
  message("Parsing trans-eQTL genes...")
  genes <- unlist(strsplit(loci$trans_genes, "|", fixed = TRUE))
  
  # Include only those trans-eQTL genes that are present in input folder
  available_genes <- str_replace(list.files(inp_folder), ".*phenotype=", "")
  genes <- genes[genes %in% available_genes]
  message("Parsing trans-eQTL genes...done!")

  # Parse SNP list
  message("Parsing SNP list...")
  snps <- unlist(strsplit(loci$SNPs, "|", fixed = TRUE))
  message("Parsing SNP list...done!")

  # Read in cis-eQTLs
  message("Reading parquet file...")
  #ds <- arrow::open_dataset(args$eqtl_folder, partitioning = "phenotype", hive_style = TRUE)
  eqtls <- arrow::open_dataset(list.files(paste0(args$eqtl_folder, "/phenotype=", args$gene_id), full.names = TRUE))

  eqtls <- eqtls %>% 
  filter(variant %in% snps) %>% 
  collect() %>%
  filter((i_squared < args$i2_thresh | is.na(i_squared)) & 
  sample_size >= args$maxN_thresh * max(sample_size) & 
  sample_size >= args$minN_thresh) %>% 
  as.data.table()

  message("cis-eQTLs read in!")
  message(paste(nrow(eqtls), "variants in the data!"))
  
  # Make beta and se matrix
  message("Make beta and se matrices...")
  eqtl_beta <- data.table(variant = eqtls$variant, eqtls$beta)
  eqtl_se <- data.table(variant = eqtls$variant, eqtls$standard_error)

  rm(eqtls)
  gc()

  if (nrow(eqtl_beta) >= 100){
  # Read iteratively in trans-eGenes
  message(paste(length(genes), "overlapping trans-eQTL genes"))
  genes_h <- genes
  for (i in 1:length(genes)){
  message(paste("Read", genes[i]))

  #eqtls2 <- read_parquet(list.files(paste0(args$eqtl_folder, "/phenotype=", genes[i]), full.names = TRUE))
  eqtls2 <- open_dataset(list.files(paste0(args$eqtl_folder, "/phenotype=", genes[i]), full.names = TRUE))

  eqtls2 <- eqtls2 %>% 
  filter(variant %in% snps) %>% 
  collect() %>% 
  filter((i_squared < args$i2_thresh | is.na(i_squared)) & 
  sample_size >= args$maxN_thresh * max(sample_size) & 
  sample_size >= args$minN_thresh) %>% 
  select(variant, beta, standard_error) %>% 
  collect() %>% 
  as.data.table()

  setkey(eqtls2, variant)
  
  if(nrow(eqtls2) > 100){ #So that there are at least some variants reasonable to run.

  colnames(eqtls2)[2:3] <- paste0(genes[i], "_",  colnames(eqtls2)[2:3])
  print(nrow(eqtls2))
  
  eqtl_beta <- merge(eqtl_beta, eqtls2[, -3, with = FALSE], by = "variant")
  eqtl_se <- merge(eqtl_se, eqtls2[, -2, with = FALSE], by = "variant")
  message(paste0(i, "/", length(genes)))
  } else {
    message("Does not make sense to include trans-gene: <100 variants after filters!")
    genes_h <- genes_h[genes_h != genes[i]]
    }
  gc()
  }
  colnames(eqtl_beta) <- c("SNP", gene, genes_h)
  colnames(eqtl_se) <- c("SNP", gene, genes_h)
  } else {
  message("Does not make sense to include cis-gene: <100 variants after filters!")
  colnames(eqtl_beta) <- c("SNP", gene)
  colnames(eqtl_se) <- c("SNP", gene)
  }
  snplist <- eqtl_beta$SNP
  
  eqtl_beta <- as.matrix(eqtl_beta[, -1, with = FALSE])
  rownames(eqtl_beta) <- snplist
  eqtl_se <- as.matrix(eqtl_se[, -1, with = FALSE])
  rownames(eqtl_se) <- snplist
  
  return(list(betas = eqtl_beta, standard_errors = eqtl_se))
  message("Make beta and se matrices...done!")
}

VisualiseLocus <- function(inp, reference, res){
  inputs <- inp
  ref <- arrow::open_dataset(reference) %>%
    select("ID", "CHR", "bp", "str_allele1", "str_allele2") %>% collect() %>% as.data.table()
  message("Reference loaded")
  setkey(ref, ID)
  
  P_table <- data.table(ZtoP(inputs$betas/inputs$standard_errors, largeZ = TRUE))
  P_table$ID <- rownames(inputs$betas)
  setkey(P_table, ID)
  
  ref2 <- ref[ID %in% P_table$ID, c(1, 3), with = FALSE]
  
  cis_gene <- colnames(P_table)[1]
  
  P_table <- merge(P_table, ref2, by = "ID")
  P_table <- P_table[order(P_table$ID)]
  P_table <- melt(P_table, id.vars = c("ID", "bp"))
  lead_var_pos <- P_table[variable == cis_gene, ]
  lead_var_pos <- lead_var_pos[value == max(value)]$bp
  
  ggplot(P_table, aes(x = bp, y = value)) + 
    geom_point(alpha = 0.4, shape = 21, colour = "black", fill = "grey") + theme_bw() +
    facet_wrap(~variable, ncol = 1, scales = "free") + ylab("-log10(P)") + 
    geom_vline(xintercept = lead_var_pos, colour = "red")
  
}

# Prepare inputs
inputs <- ParseInput(inp_folder = args$eqtl_folder, loci = args$loci, gene = args$gene)

if (ncol(inputs$betas < 2)) {

res <- data.table(cis_eQTL_gene = NA,
cis_eQTL_gene_name = NA,
type = NA,
iteration = NA,
trans_eQTL_gene = NA,
trans_eQTL_gene_name = NA,
posterior_prob = NA,
regional_prob = NA,
candidate_snp = NA,
posterior_explained_by_snp = NA,
dropped_trait = NA,
nr_snps_included = NA
)[-1]

message("Writing empty output file...")
fwrite(res, args$output, sep = "\t")
message("Writing empty output file...done!")
} else {

###########################
# Analyse cis-eQTL regions#
###########################

trait_names <- colnames(inputs$betas)
message("Running hyprcoloc analysis...")
res <- hyprcoloc(inputs$betas, 
                 inputs$standard_errors, 
                 trait.names = colnames(inputs$betas), 
                 snp.id = rownames(inputs$betas),
                 sample.overlap = matrix(1, ncol = ncol(inputs$betas), nrow = ncol(inputs$betas))
                 )  

print(res)

message("Running hyprcoloc analysis...done!")
res <- as.data.table(res$results)

res <- data.table(cis_eQTL_gene = trait_names[1], 
                       type = "trans-trans", 
                       as.data.table(res)
                  )

res <- res[res$traits != "None", ]

res[str_detect(traits, cis_eQTL_gene)]$type <- "cis-trans"
res[, traits := list(strsplit(traits, ",", fixed = TRUE))]
res$nr_snps_included <- nrow(inputs$betas)

message("Annotating results...")
gtf <- readGFF(args$gtf)
gtf <- as.data.table(unique(gtf[, c(9, 11)]))

res <- res %>% 
separate_rows(traits, sep = "\\| ") %>% 
mutate(traits = str_trim(traits))

res <- merge(res, gtf, by.x = "cis_eQTL_gene", by.y = "gene_id", all.x = TRUE)

res <- res[, c(1, 11, 2:10)]
colnames(res)[2] <- "cis_eQTL_gene_name"

res <- merge(res, gtf, by.x = "traits", by.y = "gene_id", all.x = TRUE)
res <- res[, c(2:5, 1, 12, 6:11)]
colnames(res)[c(5, 6)] <- c("trans_eQTL_gene", "trans_eQTL_gene_name")
res <- res[order(res$cis_eQTL_gene, res$type, res$iteration), ]
res <- res[res$cis_eQTL_gene != res$trans_eQTL_gene, ]

message("Annotating results...done!")
message("Writing output...")
fwrite(res, args$output, sep = "\t")
message("Writing output...done!")

#visualise
#if (isTRUE(args$locusplot)){
#VisualiseLocus(inp = inputs, reference = args$ref)
#ggsave(paste0(cis_eQTL_gene, ".pdf"), height = 15, width = 10, units = "in")
#}
}
