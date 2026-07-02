library(DESeq2)
library(dplyr)
library(tximport)

################### Preble
txi <- readRDS("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/PREB_txi_subset.rds")
samp <- as.data.frame(read.csv("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/preb_metadata.csv", head = T, row.names=1))

# design is Temp because if using Population design, we have to set a reference population
ddsTxi_temp <- DESeqDataSetFromTximport(txi,
                                   colData = samp,
                                   design = ~ Temperature)

ddsTxi_temp$Temperature <- relevel(ddsTxi_temp$Temperature, ref = "10C")

dds <- DESeq(ddsTxi_temp)

results <- results(dds)

write.csv(results,"/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/PREB_20vs10_deseq_results.csv", quote=FALSE)

save(dds, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/PREB_20vs10_deseq_object.RData")

upreg<-row.names(results[!is.na(results$padj) & results$log2FoldChange >= 2 & results$padj<0.01,])

downreg<-row.names(results[!is.na(results$padj) & results$log2FoldChange <= -2 & results$padj<0.01,])

write.table(upreg, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/degs/PREB_20vs10_up_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = F, row.names = F )

write.table(downreg, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/degs/PREB_20vs10_down_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = F, row.names = F )

###################### Boone
txi <- readRDS("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/BOON_txi_subset.rds")
samp <- as.data.frame(read.csv("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/boon_metadata.csv", head = T, row.names=1))

# design is Temp because if using Population design, we have to set a reference population
ddsTxi_temp <- DESeqDataSetFromTximport(txi,
                                   colData = samp,
                                   design = ~ Temperature)

ddsTxi_temp$Temperature <- relevel(ddsTxi_temp$Temperature, ref = "10C")

dds <- DESeq(ddsTxi_temp)

results <- results(dds)

write.csv(results,"/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/BOON_20vs10_deseq_results.csv", quote=FALSE)

save(dds, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/BOON_20vs10_deseq_object.RData")

upreg<-row.names(results[!is.na(results$padj) & results$log2FoldChange >= 2 & results$padj<0.01,])

downreg<-row.names(results[!is.na(results$padj) & results$log2FoldChange <= -2 & results$padj<0.01,])

write.table(upreg, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/degs/BOON_20vs10_up_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = F, row.names = F )

write.table(downreg, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/degs/BOON_20vs10_down_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = F, row.names = F )

###################### Scott
txi <- readRDS("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/SCOT_txi_subset.rds")
samp <- as.data.frame(read.csv("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/scot_metadata.csv", head = T, row.names=1))

# design is Temp because if using Population design, we have to set a reference population
ddsTxi_temp <- DESeqDataSetFromTximport(txi,
                                   colData = samp,
                                   design = ~ Temperature)

ddsTxi_temp$Temperature <- relevel(ddsTxi_temp$Temperature, ref = "10C")

dds <- DESeq(ddsTxi_temp)

results <- results(dds)

write.csv(results,"/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/SCOT_20vs10_deseq_results.csv", quote=FALSE)

save(dds, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/SCOT_20vs10_deseq_object.RData")

upreg<-row.names(results[!is.na(results$padj) & results$log2FoldChange >= 2 & results$padj<0.01,])

downreg<-row.names(results[!is.na(results$padj) & results$log2FoldChange <= -2 & results$padj<0.01,])

write.table(upreg, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/degs/SCOT_20vs10_up_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = F, row.names = F )
write.table(downreg, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/degs/SCOT_20vs10_down_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = F, row.names = F )

##################### WILS FG
txi <- readRDS("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/WILS_FG_txi_subset.rds")
samp <- as.data.frame(read.csv("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wils_fg_metadata.csv", head = T, row.names=1))

# design is Temp because if using Population design, we have to set a reference population
ddsTxi_temp <- DESeqDataSetFromTximport(txi,
                                   colData = samp,
                                   design = ~ Temperature)

ddsTxi_temp$Temperature <- relevel(ddsTxi_temp$Temperature, ref = "10C")

dds <- DESeq(ddsTxi_temp)

results <- results(dds)

write.csv(results,"/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/WILS_FG_20vs10_deseq_results.csv", quote=FALSE)

save(dds, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/WILS_FG_20vs10_deseq_object.RData")

upreg<-row.names(results[!is.na(results$padj) & results$log2FoldChange >= 2 & results$padj<0.01,])

downreg<-row.names(results[!is.na(results$padj) & results$log2FoldChange <= -2 & results$padj<0.01,])

write.table(upreg, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/degs/WILS_FG_20vs10_up_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = F, row.names = F )

write.table(downreg, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/degs/WILS_FG_20vs10_down_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = F, row.names = F )

##################### WILS RW
txi <- readRDS("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/WILS_RW_txi_subset.rds")
samp <- as.data.frame(read.csv("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wils_rw_metadata.csv", head = T, row.names=1))

# design is Temp because if using Population design, we have to set a reference population
ddsTxi_temp <- DESeqDataSetFromTximport(txi,
                                   colData = samp,
                                   design = ~ Temperature)

ddsTxi_temp$Temperature <- relevel(ddsTxi_temp$Temperature, ref = "10C")

dds <- DESeq(ddsTxi_temp)

results <- results(dds)

write.csv(results,"/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/WILS_RW_20vs10_deseq_results.csv", quote=FALSE)

save(dds, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/WILS_RW_20vs10_deseq_object.RData")

upreg<-row.names(results[!is.na(results$padj) & results$log2FoldChange >= 2 & results$padj<0.01,])

downreg<-row.names(results[!is.na(results$padj) & results$log2FoldChange <= -2 & results$padj<0.01,])

write.table(upreg, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/degs/WILS_RW_20vs10_up_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = F, row.names = F )

write.table(downreg, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/degs/WILS_RW_20vs10_down_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = F, row.names = F )

##################### RUTH
txi <- readRDS("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/RUTH_txi_subset.rds")
samp <- as.data.frame(read.csv("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/ruth_metadata.csv", head = T, row.names=1))

# design is Temp because if using Population design, we have to set a reference population
ddsTxi_temp <- DESeqDataSetFromTximport(txi,
                                   colData = samp,
                                   design = ~ Temperature)

ddsTxi_temp$Temperature <- relevel(ddsTxi_temp$Temperature, ref = "10C")

dds <- DESeq(ddsTxi_temp)

results <- results(dds)

write.csv(results,"/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/RUTH_20vs10_deseq_results.csv", quote=FALSE)

save(dds, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/RUTH_20vs10_deseq_object.RData")

upreg<-row.names(results[!is.na(results$padj) & results$log2FoldChange >= 2 & results$padj<0.01,])

downreg<-row.names(results[!is.na(results$padj) & results$log2FoldChange <= -2 & results$padj<0.01,])

write.table(upreg, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/degs/RUTH_20vs10_up_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = F, row.names = F )

write.table(downreg, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/degs/RUTH_20vs10_down_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = F, row.names = F )
