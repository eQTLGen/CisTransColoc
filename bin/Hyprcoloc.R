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
library(patchwork)

parser <- ArgumentParser(description = 'Run HyprColoc for every cis-eQTL locus to detect colocalisation between cis- and trans-eQTL.')

parser$add_argument('--loci', metavar = 'file', type = 'character', 
help = 'R data file which contains cis genes, trans genes and variants.')
parser$add_argument('--gene_id', type = 'character', 
help = 'ENSEMBL Gene ID for which to run the analysis.')
parser$add_argument('--eqtl_folder', metavar = 'file', type = 'character', 
help = 'eQTLGen parquet folder format with per gene output files.')
parser$add_argument('--gtf', metavar = 'file', type = 'character',
                    help = "ENSEMBL .gtf file, needs to be hg38.")
parser$add_argument('--reference', metavar = 'file', type = 'character',
                    help = 'eQTLGen SNP reference file in parquet format.')
parser$add_argument('--i2_thresh', type = 'numeric', default = 40,
                    help = 'Heterogeneity threshold. Defaults to <40%.')
parser$add_argument('--maxN_thresh', type = 'numeric', default = 0.8,
                    help = 'Per gene maximal sample size threshold. Defaults to 0.8 (SNPs with >=0.8*max(N))')
parser$add_argument('--minN_thresh', type = 'numeric', default = 0,
                    help = 'Minimal sample size threshold. Defaults to 0 (no filtering)')
parser$add_argument('--output', metavar = 'file', type = 'character', 
help = 'Output file.')
parser$add_argument('--locusplot', type = 'logical', help = "Whether to output locus plot.")
parser$add_argument('--WriteRegionOut', type = 'logical', help = "Whether to output region-specific sumstats.")

args <- parser$parse_args()

args$locusplot[args$locusplot == "true"] <- TRUE
args$locusplot[args$locusplot == "false"] <- FALSE

args$WriteRegionOut[args$WriteRegionOut == "true"] <- TRUE
args$WriteRegionOut[args$WriteRegionOut == "false"] <- FALSE

# Functions
ParseInput <- function(inp_folder, loci, gene){
  
  message("Reading locus info...")
  
  loci <- fread(loci, sep2 = "|")
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
  snps <- as.integer(unlist(strsplit(loci$SNPs, "|", fixed = TRUE)))
  message("Parsing SNP list...done!")

  # Read in cis-eQTLs
  message("Reading cis-eQTL parquet file...")
  eqtls <- arrow::open_dataset(list.files(paste0(inp_folder, "/phenotype=", gene), full.names = TRUE))

  eqtls <- eqtls %>% 
  filter(variant_index %in% snps) %>% 
  collect() %>%
  filter((i_squared <= args$i2_thresh | is.na(i_squared)) & 
  sample_size >= args$maxN_thresh * max(sample_size) & 
  sample_size >= args$minN_thresh) %>% 
  as.data.table()

  message("Reading cis-eQTL parquet file...done!")
  message(paste(nrow(eqtls), "variants in the data!"))
  
  # Make beta and se matrix
  message("Make beta and se matrices...")
  eqtl_beta <- data.table(variant = eqtls$variant_index, eqtls$beta)
  eqtl_se <- data.table(variant = eqtls$variant_index, eqtls$standard_error)

  rm(eqtls)
  gc()

  if (nrow(eqtl_beta) >= 100){
  # Read iteratively in trans-eGenes
  message(paste(length(genes), "overlapping trans-eQTL genes"))
  genes_h <- genes
  for (i in 1:length(genes)){
  message(paste("Read", genes[i]))

  eqtls2 <- open_dataset(list.files(paste0(args$eqtl_folder, "/phenotype=", genes[i]), full.names = TRUE))

  if(length(genes_h) < 5){

  eqtls2 <- eqtls2 %>% 
  filter(variant_index %in% snps) %>% 
  collect() %>% 
  filter((i_squared < args$i2_thresh | is.na(i_squared)) & 
  sample_size >= args$maxN_thresh * max(sample_size) & 
  sample_size >= args$minN_thresh) %>% 
  select(variant_index, beta, standard_error) %>% 
  collect() %>% 
  as.data.table()

  } else {

  eqtls2 <- eqtls2 %>% 
  filter(variant_index %in% snps) %>% 
  collect() %>% 
  select(variant_index, beta, standard_error) %>% 
  collect() %>% 
  as.data.table()

  }

  setkey(eqtls2, variant_index)
  
  if(nrow(eqtls2) > 100){ #So that there are at least some variants reasonable to run.

  colnames(eqtls2)[2:3] <- paste0(genes[i], "_",  colnames(eqtls2)[2:3])
  print(nrow(eqtls2))
  
  eqtl_beta <- merge(eqtl_beta, eqtls2[, -3, with = FALSE], by.x = "variant", by.y = "variant_index")
  eqtl_se <- merge(eqtl_se, eqtls2[, -2, with = FALSE], by.x = "variant", by.y = "variant_index")
  message(paste(nrow(eqtl_beta), "variants in the combined matrix"))
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


VisualiseLocus <- function(inputs, reference, gtf, ld = NULL){

  PlotRegionGenes <- function(chr, pos_start, pos_end, gtf, line, limit1, limit2){

    gtf <- gtf[gtf$type %in% c("gene"), c(3, 1, 4, 5, 7, 9, 11, 13, 22)]
    gtf <- as.data.table(gtf)

    gtf_f <- gtf[seqid %in% chr & start > pos_start & end < pos_end]
    gtf_f <- gtf_f[gtf_f$gene_biotype %in% c("protein_coding", "lncRNA", "miRNA"), ]
    gtf_f$line <- c(rep(seq(from = 1, to = 10), times = floor(nrow(gtf_f)/10)), 1:(nrow(gtf_f) - floor(nrow(gtf_f)/10) * 10))

    gtf_f$start2 <- gtf_f$start
    gtf_f$end2 <- gtf_f$end

    gtf_f[gtf_f$strand == "-", ]$start2 <- gtf_f[gtf_f$strand == "-", ]$end
    gtf_f[gtf_f$strand == "-", ]$end2 <- gtf_f[gtf_f$strand == "-", ]$start
    gtf_f[is.na(gtf_f$gene_name), ]$gene_name <- gtf_f[is.na(gtf_f$gene_name), ]$gene_id

    p <- ggplot(gtf_f, aes(x = start2, y = line, xend = end2, yend = line, colour = gene_biotype, label = gene_name)) +
      geom_segment(arrow = arrow(length = unit(0.20, "cm"), ends = "last", type = "closed")) +
      theme_bw() + guides(y = "none") + labs(y = NULL) + xlab(paste0("chr", unique(gtf_f$seqid), " position (bp)")) +
      geom_text(size = 1.5, nudge_y = 0.4) +
      scale_color_manual(values = c("protein_coding" = "salmon",
                                    "lncRNA" = "steelblue",
                                    "miRNA" = "darkgreen")) +
      scale_x_continuous(limits = c(limit1, limit2)) +
      guides(color = "none") +
      geom_vline(xintercept = line, colour = "red")

    return(p)

  }


  message("Loading variant reference...")

  snp_ind <- rownames(inputs$betas)

  ref <- arrow::open_dataset(reference) %>%
    filter(variant_index %in% snp_ind) %>%
    select("variant_index", "chromosome", "bp", "non_eff_allele", "eff_allele") %>%
    collect() %>%
    as.data.table()

  message("Loading variant reference...done!")

  message("Indexing variant reference...")
  setkey(ref, variant_index)
  message("Indexing variant reference...done!")

  message("Loading gene reference...")
  gtf <- readGFF(gtf)
  gtf2 <- as.data.table(unique(gtf[, c(9, 11)]))
  message("Loading gene reference...done!")

  P_table <- data.table(ZtoP(inputs$betas/inputs$standard_errors, largeZ = TRUE))
  P_table$variant_index <- rownames(inputs$betas)
  P_table$variant_index <- as.integer(P_table$variant_index)

  setkey(P_table, variant_index)

  ref2 <- ref[variant_index %in% P_table$variant_index, c(1, 2, 3), with = FALSE]
  ref2$variant_index <- as.integer(ref2$variant_index)

  cis_gene <- colnames(P_table)[1]

  P_table <- merge(P_table, ref2, by = "variant_index")
  P_table <- P_table[order(P_table$variant_index)]
  P_table <- melt(P_table, id.vars = c("variant_index", "bp", "chromosome"))
  lead_var_pos <- P_table[variable == cis_gene, ]
  lead_var_pos <- lead_var_pos[value == max(value)]$bp

  P_table <- merge(P_table, gtf2, by.x = "variable", by.y = "gene_id")
  rm(gtf2)
  P_table$type <- "trans"
  P_table[variable == cis_gene, ]$type <- "cis"
  P_table[, variable := paste0(gene_name, " ", variable)]

  if (!is.null(ld)){
    if (nrow(ld) < 1){message("Proxies file is empty!")}
    P_table <- merge(P_table, ld, by = "variant_index", all.x = TRUE)
    P_table[is.na(R2), ]$R2 <- 0

    p1 <- ggplot(P_table, aes(x = bp, y = value)) +
      geom_point(shape = 21, colour = "black", aes(fill = R2)) + theme_bw() +
      facet_wrap(~type * variable, ncol = 1, scales = "free") + ylab("-log10(P)") + xlab("") +
      geom_vline(xintercept = lead_var_pos, colour = "red") + guides(fill = "none") +
      scale_fill_gradient(low = "white", high = "red")

  } else {
  p1 <- ggplot(P_table, aes(x = bp, y = value)) +
    geom_point(alpha = 0.4, shape = 21, colour = "black", fill = "grey") + theme_bw() +
    facet_wrap(~type * variable, ncol = 1, scales = "free") + ylab("-log10(P)") + xlab("") +
    geom_vline(xintercept = lead_var_pos, colour = "red")
  }
  p2 <- PlotRegionGenes(chr = unique(P_table$chromosome), pos_start = min(P_table$bp), pos_end = max(P_table$bp),
                        gtf = gtf,
                        line = lead_var_pos,
                        limit1 = min(P_table$bp), limit2 = max(P_table$bp))

  p <- p1 / p2 + plot_layout(heights = c(6, 1))

  return(list(p, lead_var_pos))
}

# Prepare inputs
inputs <- ParseInput(inp_folder = args$eqtl_folder, loci = args$loci, gene = args$gene)

if (ncol(inputs$betas) < 2) {

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

# Remove rows where SE is 0
# TODO: check why such variants are in the results
any_row_contains_zero <- apply(inputs$standard_errors, 1, function(row) any(row != 0))


message("Running hyprcoloc analysis...")
res <- hyprcoloc(inputs$betas[any_row_contains_zero, ], 
                 inputs$standard_errors[any_row_contains_zero, ], 
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
gtf <- as.data.table(gtf)
gtf2 <- as.data.table(unique(gtf[, c(9, 11)]))

res <- res %>% 
separate_rows(traits, sep = "\\| ") %>% 
mutate(traits = str_trim(traits))

res <- merge(res, gtf2, by.x = "cis_eQTL_gene", by.y = "gene_id", all.x = TRUE)

res <- res[, c(1, 11, 2:10)]
colnames(res)[2] <- "cis_eQTL_gene_name"

res <- merge(res, gtf2, by.x = "traits", by.y = "gene_id", all.x = TRUE)
res <- res[, c(2:5, 1, 12, 6:11)]
colnames(res)[c(5, 6)] <- c("trans_eQTL_gene", "trans_eQTL_gene_name")
res <- res[order(res$cis_eQTL_gene, res$type, res$iteration), ]
res <- res[res$cis_eQTL_gene != res$trans_eQTL_gene, ]

message("Annotating results...done!")
message("Writing output...")
fwrite(res, args$output, sep = "\t")
message("Writing output...done!")

# Visualise
if (isTRUE(args$locusplot)){
  message("Writing out plot(s)...")
  if (ncol(inputs$betas) <= 10){
    p <- VisualiseLocus(inputs = inputs, reference = args$reference, gtf = gtf)
    ggsave(paste0(args$gene, ".png"), height = 15, width = 10, units = "in")
  } else {
      # Define the number of columns you want to select each time (excluding the first column)
      num_cols <- 9

      # Get the total number of columns in your dataframe
      total_cols <- ncol(inputs$betas)

      # Create a loop that selects 1 + num_cols at a time (first column + next num_cols columns)
      for (i in seq(2, total_cols, by = num_cols)) {
        # Calculate the range of columns for the current iteration
        col_range <- i:min(i + num_cols - 1, total_cols)  # Ensure we don't exceed the total number of columns
        inputs_temp <- list(
          betas = inputs$betas[, c(1, col_range)],
          standard_errors = inputs$standard_errors[, c(1, col_range)]
        )

      p <- VisualiseLocus(inputs = inputs_temp, reference = args$reference, gtf = gtf)
      ggsave(paste0(args$gene, "_", col_range, ".png"), height = 15, width = 10, units = "in")
      }

  }
  message("Writing out plot(s)...done!")
}
}

# Write out regional sumstats
if (isTRUE(args$WriteRegionOut)){
  message("Writing out region-specific-sumstats...")
  saveRDS(inputs, file = paste0(args$gene, "_region.Rds"))
  message("Writing out region-specific-sumstats...done!")
}
