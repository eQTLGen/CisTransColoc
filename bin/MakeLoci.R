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
 parser$add_argument('--cis_gene_filter', metavar = 'file', type = 'character',
                     help = 'File for filtering cis-eQTL genes included to analysis.')
# parser$add_argument('--max_lead_distance', type = 'numeric', default = 250000,
#                     help = 'Maximum distance between primary cis and trans lead variants.')

args <- parser$parse_args()


message("Reading in reference...")
ref <- arrow::open_dataset(args$reference) %>% 
select("ID", "CHR", "bp", "str_allele1", "str_allele2") %>% collect() %>% as.data.table()

message(paste("Reference read"))

message("Reading in sig. results...")
sig <- fread(args$sig_res, key = "SNP")

message("Filter input to predefined cis genes.")
cis_filter <- fread(args$cis_gene_filter, header = FALSE)

if(nrow(cis_filter) > 0) {

    sig <- sig[phenotype %in% cis_filter$V1]
    message(paste("Cis-gene filter active:", length(unique(cis_filter$V1)), "genes."))

}else{

    message("No cis-gene filtering done.")

}

sig <- merge(sig, ref[, c(1:3), with = FALSE], by.x = "SNP", by.y = "ID")
message(paste("Sig. results read"))


message("Done")
message("Filter results to genes available in full files...")
eqtl_genes <- str_replace(list.files(args$eqtl), ".*phenotype=", "")
sig <- sig[phenotype %in% eqtl_genes]
message(length(unique(sig$phenotype)))
message("Done!")

message("Finding lead variants for each gene...")
LeadVariants <- sig %>% 
    group_by(phenotype) %>% 
    group_modify(~ IdentifyLeadSNPs(.x, 
    snp_id_col = "SNP", 
    snp_chr_col = "CHR", 
    snp_pos_col = "bp", 
    eff_all_col = "alt_all", 
    other_all_col = "ref_all",
    beta_col = "beta", 
    se_col = "se", 
    p_col = "P", 
    window = args$lead_variant_win))
message(paste("Lead variants found"))

rm(sig)
gc()

# Annotate cis/trans
message("Annotating lead variants cis/trans...")
ensg <- readGFF(args$gtf)
ensg <- as.data.table(ensg)
ensg <- ensg[type == "gene"]
ensg$tss <- ensg$start
ensg[strand == "-"]$tss <- ensg[strand == "-"]$end
ensg <- ensg[, c(9, 1, 27), with = FALSE]

Lead2 <-  merge(LeadVariants, ensg, by.x = "phenotype", by.y = "gene_id")
Lead2$type <- "cis"
Lead2[Lead2$chr != Lead2$seqid | abs(Lead2$pos - Lead2$tss) > args$trans_win, ]$type <- "trans"
Lead2 <- Lead2[Lead2$type == "cis" | (Lead2$type == "trans" & Lead2$P < 1e-12), ]

message(paste("Lead variants annotated"))

# For every cis locus find overlapping trans loci
message("Finding overlapping trans loci...")
cis_genes <- unique(Lead2[Lead2$type == "cis",]$phenotype)

cis_genes <- Lead2 %>% 
filter(type == "cis") %>%
group_by(phenotype) %>% 
filter(abs(Z) == max(abs(Z))) %>% 
mutate(region_start = tss - args$cis_win,
region_end = tss + args$cis_win) %>% as.data.table()

Lead2 <- as.data.table(Lead2)
cis_genes <- as.data.table(cis_genes)

setkeyv(ref, c("CHR", "bp"))

cis_genes <- data.table(cis_gene = cis_genes$phenotype, 
cis_SNP = cis_genes$SNP, chr = cis_genes$seqid, start = cis_genes$region_start, 
end = cis_genes$region_end)

cis_genes$chr <- as.factor(cis_genes$chr)
setkey(cis_genes, chr, start, end)

# overlap with trans genes
trans_genes <- Lead2 %>% 
filter(type == "trans") %>%
mutate(snp_start = pos,
snp_end = pos + 1) %>% as.data.table()

trans_genes <- data.table(trans_gene = trans_genes$phenotype, 
trans_SNP = trans_genes$SNP, chr = trans_genes$chr, start = trans_genes$snp_start, 
end = trans_genes$snp_end)

trans_genes$chr <- as.factor(trans_genes$chr)
setkey(trans_genes, chr, start, end)

cis_overlaps <- foverlaps(trans_genes, cis_genes, type = "within", nomatch = NULL)

cis_overlaps <- cis_overlaps[, .(cis_gene = unique(cis_gene), cis_SNP = unique(cis_SNP), 
                                 trans_genes = list(unique(trans_gene)), chr = unique(chr), 
                                 start = unique(start), end = unique(end)), by = cis_gene]

setkey(cis_overlaps, chr, start, end)
message("Overlapping trans loci found")

ref <- data.table(SNP = ref$ID, chr = ref$CHR, start = ref$bp, 
end = ref$bp + 1)
ref$chr <- as.factor(ref$chr)
setkey(ref, chr, start, end)


if (length(cis_overlaps$cis_gene) <= 1000){batches = 1}else{
batches <- c(seq(from = 1, to = length(cis_overlaps$cis_gene), by = 1000), length(cis_overlaps$cis_gene))
}

res <- data.table(cis_SNP = NA, SNPs = NA)[-1]

rm(LeadVariants)
gc()

message(paste0("Iteratively finding overlapping variants in ", length(batches), " batches..."))

for (i in 1:length(batches)){

message(paste0("Overlapping batch ", i, "..."))
if(i < length(batches) & i > 1){
    temp_cis_overlaps <- cis_overlaps[batches[i]:batches[i+1]]
    } else if(length(batches) == 1 & i == 1){
        temp_cis_overlaps <- cis_overlaps
    } else {
        message("Debug!")
        }

print(head(temp_cis_overlaps))

ref_temp <- ref[chr %in% temp_cis_overlaps$chr]

overlap_temp <- foverlaps(ref_temp, temp_cis_overlaps[, -c(1, 4), with = FALSE], type = "within", nomatch = NULL)

overlap_temp <- overlap_temp[, .(cis_SNP = unique(cis_SNP), SNPs = list(unique(SNP))), by = cis_gene]
res <- rbind(res, overlap_temp[, -1, with = FALSE])

message(paste0("Done"))

}

res <- merge(cis_overlaps[, c(3, 2, 4), with = FALSE], res, by = "cis_SNP")
message("Saving results...")
fwrite(res, file = "cis_trans_info.txt")
message("Done")
