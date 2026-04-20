library(DESeq2)
library(dplyr)

genes <- as.data.frame(read.csv("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/temp/input_subsets/RUTH_counts_matrix.csv", head = T, row.names=1))

samp <- as.data.frame(read.csv("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/temp/input_subsets/ruth_metadata.csv", head = T, row.names=1))

## Define your model 
pre_dds <- DESeqDataSetFromMatrix(countData = genes,
                              colData = samp,
                              design= ~ Temperature)
                              
## Relevel the treatment factor to ensure the correct reference level
pre_dds$Temperature <- relevel(pre_dds$Temperature, ref = "10C")
                              
## Run the thing
dds <- DESeq(pre_dds)

#names <- print(resultsNames(dds))
#write.csv(names, file="/work/projects/5-39268/mgaupp42/deseq_post_evigene/gene_level/result_names1.csv")

#### get the results ####
results <- results(dds)

write.csv(results,"/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/temp/results_deseq/RUTH_20vs10_deseq_results.csv", quote=FALSE)

save(dds, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/temp/results_deseq/RUTH_20vs10_deseq_object.RData")

#library(DESeq2)
library(ggplot2)

## Load the DESeq2 object
#load("/work/projects/5-39268/mgaupp42/deseq_post_evigene/gene_levelpop_deseq_object.RData")

## Perform PCA
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Clutch", "Temperature"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

## Ensure Temperature is treated as a factor
pcaData$Temperature <- as.factor(pcaData$Temperature)

## Plot PCA
p <- ggplot(pcaData, aes(PC1, PC2, color=Clutch, shape=Temperature)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  coord_fixed(ratio = 2) 

## Save the plot
#ggsave("/work/projects/5-39268/mgaupp42/deseq_post_evigene/population/pca_plot1.png", plot=p)
ggsave("/work/projects/5-39268/mgaupp42//diff_express/deseq/gene_level/temp/deseq_pcas/ruth_pca_plot.png", plot=p, width=5, height=5)

####################### Export the up and down regulated genes ################
#library(DESeq2)
library(ggplot2)
#library(dplyr)

# Load the results
results <- read.csv("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/temp/results_deseq/RUTH_20vs10_deseq_results.csv", row.names=1, head=T)

upreg<-row.names(results[!is.na(results$padj) & results$log2FoldChange >= 2 & results$padj<0.01,])

downreg<-row.names(results[!is.na(results$padj) & results$log2FoldChange <= -2 & results$padj<0.01,])

write.table(upreg, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/temp/DEGs_deseq/RUTH_20vs10_up_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = F, row.names = F )

write.table(downreg, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/temp/DEGs_deseq/RUTH_20vs10_down_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = F, row.names = F )
