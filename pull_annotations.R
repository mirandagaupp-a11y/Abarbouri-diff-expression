library(readr)
library(dplyr)

# Input files
deg_file <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/RUTH_top20_genes.txt"
trinotate_file <- "/work/projects/5-39268/mgaupp42/transcriptome_assembly/trinotate_post_ncbi/one_trans_per_gene_annotation_filtered.tsv"

# Read DEG IDs (one per line)
deg_ids <- read_lines(deg_file)

# Read Trinotate annotations
trinotate <- read_tsv(
  trinotate_file,
  show_col_types = FALSE
)

# Filter
deg_annotations <- trinotate %>%
  filter(gene_id %in% deg_ids)

# Write output
write_tsv(deg_annotations, "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/annotations/rutherford_top20_genes_annotations.tsv")

# Print summary
cat("DEG IDs:", length(deg_ids), "\n")
cat("Annotations found:", nrow(deg_annotations), "\n")