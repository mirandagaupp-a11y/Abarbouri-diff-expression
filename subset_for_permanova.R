# Load necessary libraries
library(tximport)

# Read the tximport object
txi <- readRDS(
  "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/txi_genelevel.rds"
)

# Read the list of genes to include
genes_to_include <- readLines(
  "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/degs/all_degs.txt"
)

# Find genes present in txi
genes_keep <- intersect(rownames(txi$counts), genes_to_include)

# Subset all tximport matrices consistently
txi_subset <- txi

txi_subset$counts <- txi$counts[genes_keep, , drop = FALSE]
txi_subset$abundance <- txi$abundance[genes_keep, , drop = FALSE]
txi_subset$length <- txi$length[genes_keep, , drop = FALSE]

# Save subsetted tximport object
saveRDS(
  txi_subset,
  "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/permanova/genelevel_permanova_subset.rds"
)