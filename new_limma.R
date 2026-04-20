# Load necessary libraries
library(limma)
library(edgeR)
library(ggplot2)
library(dplyr)

# Read expression matrix and metadata
genes <- as.data.frame(read.csv("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/temp/input_subsets/RUTH_counts_matrix.csv", head = TRUE, row.names = 1))
samp_df <- as.data.frame(read.csv("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/temp/input_limma/ruth_metadata.csv", head = TRUE, row.names = 1))

# Convert to factors and relevel
samp_df$Temperature <- factor(samp_df$Temperature)
samp_df$Temperature <- relevel(samp_df$Temperature, ref = "Temperature10C")  # Make 10C the reference
samp_df$Clutch <- factor(samp_df$Clutch)

# Create design matrix (includes intercept)
design <- model.matrix(~ Temperature, data = samp_df)

# Create DGE object and normalize
dge <- DGEList(counts = genes)
dge <- calcNormFactors(dge)

# Initial voom transformation
v <- voom(dge, design)

# Estimate within-clutch correlation
corfit <- duplicateCorrelation(v, design, block = samp_df$Clutch)

# Recalculate voom with correlation
v <- voom(dge, design, block = samp_df$Clutch)

# Fit model using clutch as a random effect
fit <- lmFit(v, design, block = samp_df$Clutch, correlation = corfit$consensus)
fit <- eBayes(fit)

print(colnames(fit$coefficients))

# Extract results for Temperature20C vs 10C (10C is intercept)
all_results <- topTable(fit, coef = "TemperatureTemperature20C", adjust = "fdr", number = Inf)
write.csv(all_results, file = "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/temp/results_limma/NEW_RUTH_20vs10_limma_results.csv", quote = FALSE)

# Filter DE genes (logFC >= 2 or <= -2 and adj.P.Val < 0.01)
de_genes <- all_results[
  (all_results$logFC >= 2 & all_results$adj.P.Val < 0.01) |
  (all_results$logFC <= -2 & all_results$adj.P.Val < 0.01), 
]
write.csv(de_genes, file = "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/temp/DEGs_limma/RUTH_20vs10_DEGs.csv", quote = FALSE)

# Export upregulated and downregulated genes
upreg <- row.names(de_genes[de_genes$logFC >= 2, ])
downreg <- row.names(de_genes[de_genes$logFC <= -2, ])

write.table(upreg, file = "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/temp/DEGs_limma/NEW_RUTH_20vs10_up_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downreg, file = "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/temp/DEGs_limma/NEW_RUTH_20vs10_down_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
