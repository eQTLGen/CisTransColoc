library(data.table)
library(rtracklayer)
library(tidyr)

res <- fread("/Users/urmovosa/Documents/projects/2019/eQTLGenPhase2_pipelines/eQTLGenCisTransColoc/tests/output/CisTransColocResults.txt")
gtf <- readGFF("/Users/urmovosa/Documents/projects/2019/eQTLGenPhase2_pipelines/eQTLGenCisTransColoc/tests/Homo_sapiens.GRCh38.108.gtf.gz")
gtf <- as.data.table(unique(gtf[, c(9, 11)]))

res <- res %>% 
separate_rows(traits, sep = "\\| ")

res <- merge(res, gtf, by.x = "cis_eQTL_gene", by.y = "gene_id", all.x = TRUE)

res <- res[, c(1, 11, 2:10)]
colnames(res)[2] <- "cis_eQTL_gene_name"

res <- merge(res, gtf, by.x = "traits", by.y = "gene_id", all.x = TRUE)
res <- res[, c(2:5, 1, 12, 6:11)]
colnames(res)[c(5, 6)] <- c("trans_eQTL_gene", "trans_eQTL_gene_name")
res <- res[order(res$cis_eQTL_gene, res$type, res$iteration), ]
res <- res[res$cis_eQTL_gene != res$trans_eQTL_gene, ]