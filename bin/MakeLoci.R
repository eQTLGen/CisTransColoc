#!/usr/bin/env Rscript

library(data.table)
library(stringr)
library(arrow)
library(dplyr)
library(rtracklayer)
library(IGUtilityPackage)
library(argparse)

setDTthreads(1)
parser <- ArgumentParser(description = 'Overlap cis and trans loci.')

parser$add_argument('--reference', metavar = 'file', type = 'character',
                    help = 'eQTLGen SNP reference file in parquet format.')
parser$add_argument('--sig_res', metavar = 'file', type = 'character',
                    help = 'Significant eQTL results in eQTLGen format. Can be gzipped.')
parser$add_argument('--eqtl_folder', metavar = 'file', type = 'character',
                    help = 'Parquet folder with all files.')
parser$add_argument('--gtf', metavar = 'file', type = 'character',
                    help = "ENSEMBL .gtf file, needs to be hg38.")
parser$add_argument('--lead_variant_win', type = 'numeric', default = 1000000,
                    help = 'Distance threshold for distance pruning.')
parser$add_argument('--cis_win', type = 'numeric', default = 1000000,
                    help = 'Distance threshold around cis eQTL lead variant. This window is used in coloc analysis.')
parser$add_argument('--trans_win', type = 'numeric', default = 5000000,
                    help = 'Distance threshold to declare variant trans.')
parser$add_argument('--p_thresh', type = 'numeric', default = 5e-8,
                    help = 'P-value threshold for significant effects.')
parser$add_argument('--i2_thresh', type = 'numeric', default = 40,
                    help = 'Heterogeneity threshold. Defaults to <40%.')
parser$add_argument('--maxN_thresh', type = 'numeric', default = 0.8,
                    help = 'Per gene maximal sample size threshold. Defaults to 0.8 (SNPs with >=0.8*max(N))')
parser$add_argument('--minN_thresh', type = 'numeric', default = 0,
                    help = 'Minimal sample size threshold. Defaults to 0 (no filtering)')
parser$add_argument('--cis_gene_filter', metavar = 'file', type = 'character',
                     help = 'File for filtering cis-eQTL genes included to analysis.')

args <- parser$parse_args()

# TEMP
args <- list(
    reference = "/Users/urmovosa/Documents/projects/2019/eQTLGenPhase2/projectfiles/data/derived_data/References/1000G-30x_index.parquet",
    sig_res = "/Users/urmovosa/Documents/projects/2019/eQTLGenPhase2_pipelines/eQTLGenCisTransColoc/tests/input/SigResults/subset_p5e8_hyprColocFormat_2024-09-05.csv.gz",
    eqtl_folder = "/Users/urmovosa/Documents/projects/2019/eQTLGenPhase2_pipelines/eQTLGenCisTransColoc/tests/input/sumstats/",
    gtf = "/Users/urmovosa/Documents/projects/2019/eQTLGenPhase2_pipelines/eQTLGenCisTransColoc/tests/Homo_sapiens.GRCh38.108.gtf.gz",
    lead_variant_win = 1000000,
    cis_win = 1000000,
    trans_win = 5000000,
    p_thresh = 5e-8,
    i2_thresh = 100,
    maxN_thresh = 0.8,
    minN_thresh = 0,
    cis_gene_filter = "/Users/urmovosa/Documents/projects/2019/eQTLGenPhase2_pipelines/eQTLGenCisTransColoc/data/help_input.txt"
)


message("Reading in sig. results...")
sig <- fread(args$sig_res, key = "variant_index")

message("Filter results to genes available in full files...")
eqtl_genes <- str_replace(list.files(args$eqtl), ".*phenotype=", "")
sig <- sig[phenotype %in% eqtl_genes]
message(paste(length(unique(sig$phenotype)), "genes in full results"))
if (length(unique(sig$phenotype)) < 2){stop("Less than two genes in the full results, terminating.")}
message("Filter results to genes available in full files...done!")

sig <- sig[P < args$p_thresh & (i_squared <= args$i2_thresh | is.na(i_squared))]

sig <- sig %>% 
    group_by(phenotype) %>% 
    filter(N >= args$maxN_thresh * max(N) & N >= args$minN_thresh) %>% 
    as.data.table()

message(paste(nrow(sig), "rows among significant results"))

message("Reading in sig. results...done!")

message("Reading in SNP list in the full files...")
snp_list  <- arrow::open_dataset(list.files(args$eqtl_folder, full.names = TRUE)[1])
snp_list <- snp_list %>% select(variant_index) %>% collect() %>% as.data.table()
message("Reading in SNP list in the full files...done!")

message("Reading in reference...")

ref <- arrow::open_dataset(args$reference)

ref <- ref %>% 
    filter(variant_index %in% !!snp_list$variant_index) %>% 
    collect()

ref <- data.table(ref, key = "variant_index")
ref <- ref[, c(6, 3, 2, 4, 5), with = FALSE]
message("Reading in reference...done!")

sig <- merge(sig, ref[, c(1:3), with = FALSE], by = "variant_index")
message(paste(nrow(sig), "rows among significant results, after merging with reference"))

message("Finding lead variants for each gene...")
LeadVariants <- sig %>% 
    group_by(phenotype) %>% 
    group_modify(~ IdentifyLeadSNPs(.x, 
    snp_id_col = "variant_index", 
    snp_chr_col = "chromosome", 
    snp_pos_col = "bp", 
    eff_all_col = "alt_all", 
    other_all_col = "ref_all",
    beta_col = "beta", 
    se_col = "se", 
    p_col = "P", 
    window = args$lead_variant_win))

message(paste(nrow(LeadVariants), "loci"))    
message("Finding lead variants for each gene...done!")

rm(sig)
gc()

message("Removing lead variants mapping to MHC region...")
message(paste(nrow(LeadVariants[LeadVariants$chr == 6 & LeadVariants$pos > 25000000 & LeadVariants$pos < 34000000, ]), "such variants removed."))
LeadVariants <- LeadVariants[!(LeadVariants$chr == 6 & LeadVariants$pos > 25000000 & LeadVariants$pos < 34000000), ]
message("Removing lead variants mapping to MHC region...done!")

# Annotate cis/trans
message("Annotating lead variants to cis/trans...")
ensg <- readGFF(args$gtf)
ensg <- as.data.table(ensg)
ensg <- ensg[type == "gene"]
ensg$tss <- ensg$start
ensg[strand == "-"]$tss <- ensg[strand == "-"]$end
ensg <- ensg[, c(9, 1, 27), with = FALSE]

Lead2 <-  merge(LeadVariants, ensg, by.x = "phenotype", by.y = "gene_id")
Lead2$type <- "interim"
Lead2[Lead2$chr == Lead2$seqid | abs(Lead2$pos - Lead2$tss) < args$cis_win, ]$type <- "cis"
Lead2[Lead2$chr != Lead2$seqid | abs(Lead2$pos - Lead2$tss) > args$trans_win, ]$type <- "trans"
Lead2 <- Lead2[Lead2$type == "cis" | (Lead2$type == "trans" & Lead2$P < args$p_thresh), ]

fwrite(Lead2, "eQtlLeadVariants.txt.gz", sep = "\t")

message("Annotating lead variants cis/trans...done!")

# For every cis locus find overlapping trans loci
message("Finding overlapping trans loci...")
cis_genes <- unique(Lead2[Lead2$type == "cis",]$phenotype)

cis_genes <- Lead2 %>% 
filter(type == "cis") %>%
group_by(phenotype) %>% 
filter(abs(Z) == max(abs(Z))) %>% 
filter(abs(pos - tss) == min(abs(pos - tss))) %>% 
mutate(region_start = tss - args$cis_win,
region_end = tss + args$cis_win) %>% as.data.table()

message("Filter input to predefined cis genes.")
cis_filter <- fread(args$cis_gene_filter, header = TRUE)

if(nrow(cis_filter) > 0) {

    cis_genes <- cis_genes[phenotype %in% cis_filter$V1]
    message(paste("Cis-gene filter active:", length(unique(cis_filter$V1)), "genes."))

}else{

    message("No cis-gene filtering done.")

}

Lead2 <- as.data.table(Lead2)
cis_genes <- as.data.table(cis_genes)

setkeyv(ref, c("chromosome", "bp"))

cis_genes <- data.table(
    cis_gene = cis_genes$phenotype, 
    cis_SNP = cis_genes$SNP, 
    chr = cis_genes$seqid, 
    start = cis_genes$region_start, 
    end = cis_genes$region_end)

cis_genes$chr <- as.factor(cis_genes$chr)
setkey(cis_genes, chr, start, end)

# overlap with trans genes
trans_genes <- Lead2 %>% 
filter(type == "trans") %>%
mutate(snp_start = pos,
snp_end = pos + 1) %>% as.data.table()

trans_genes <- data.table(
    trans_gene = trans_genes$phenotype, 
    trans_SNP = trans_genes$SNP, 
    chr = trans_genes$chr, 
    start = trans_genes$snp_start, 
    end = trans_genes$snp_end)

trans_genes$chr <- as.factor(trans_genes$chr)
setkey(trans_genes, chr, start, end)

print(head(cis_genes))
print(head(trans_genes))

cis_overlaps <- foverlaps(trans_genes, cis_genes, type = "within", nomatch = NULL)

cis_overlaps <- cis_overlaps[, .(
    cis_gene = unique(cis_gene), 
    cis_SNP = unique(cis_SNP), 
    trans_genes = list(unique(trans_gene)), 
    chr = unique(chr), 
    start = unique(start), 
    end = unique(end)), 
    by = cis_gene]

setkey(cis_overlaps, chr, start, end)
message("Finding overlapping trans loci...done!")

ref <- data.table(SNP = ref$variant_index, chr = ref$chromosome, start = ref$bp, 
end = ref$bp + 1)
ref$chr <- as.factor(ref$chr)
setkey(ref, chr, start, end)

if (length(cis_overlaps$cis_gene) <= 1000){
batches = 1
}else{
batches <- c(seq(from = 1, to = length(cis_overlaps$cis_gene), by = 1000), length(cis_overlaps$cis_gene))
}

res <- data.table(cis_gene = NA, cis_SNP = NA, SNPs = NA)[-1]

rm(LeadVariants)
gc()

message(paste0("Iteratively finding overlapping variants in ", length(batches), " batches..."))

for (i in 1:length(batches)){

message(paste0("Overlapping batch ", i, "..."))
if(i < length(batches) & length(batches) > 1 & i != 1){
        temp_cis_overlaps <- cis_overlaps[batches[i]:batches[i+1]]
        message("Multiple batches")
    } else if(length(batches) > 1 & i == 1){
        temp_cis_overlaps <- cis_overlaps[batches[i]]
        message("Multiple batches")
    } else if(length(batches) > 1 & i == length(batches)){
        temp_cis_overlaps <- cis_overlaps[batches[i]]
        message("Multiple batches")
    } else if(length(batches) == 1 & i == 1){
        temp_cis_overlaps <- cis_overlaps
        message("Single batch")
    } else {
        message("Debug!")
        }

ref_temp <- ref[chr %in% temp_cis_overlaps$chr]

head(temp_cis_overlaps[, -c(1, 4), with = FALSE])

overlap_temp <- foverlaps(ref_temp, temp_cis_overlaps[, -c(1, 4), with = FALSE], type = "within", nomatch = NULL)

overlap_temp <- overlap_temp[, .(cis_SNP = unique(cis_SNP), SNPs = list(unique(SNP))), by = cis_gene]
res <- rbind(res, overlap_temp)

message(paste0("Done"))

}

res <- merge(cis_overlaps[, c(3, 2, 4), with = FALSE], res, by = c("cis_SNP", "cis_gene"))
#TODO: debug why few duplicated rows were introduced
res <- unique(as.data.frame(res))

if(nrow(res[duplicated(res$cis_gene), ]) == 0){

message("Saving results...")
fwrite(res, file = "cis_trans_info.txt")
message("Saving results...done!")

} else {message("Debug!")}
