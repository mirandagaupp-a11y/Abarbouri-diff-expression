# Output directory
out_dir <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results"

# Input files and population names
files <- c(
  RUTH = file.path(out_dir, "RUTH_20vs10_deseq_results.csv"),
  RW   = file.path(out_dir, "WILS_RW_20vs10_deseq_results.csv"),
  FG   = file.path(out_dir, "WILS_FG_20vs10_deseq_results.csv"),
  SCOT = file.path(out_dir, "SCOT_20vs10_deseq_results.csv"),
  BOON = file.path(out_dir, "BOON_20vs10_deseq_results.csv"),
  PREB = file.path(out_dir, "PREB_20vs10_deseq_results.csv")
)

for (pop in names(files)) {
  
  # Load DESeq2 results
  data <- read.csv(files[pop], row.names = 1)
  
  # Add gene IDs as a column
  data$gene_id <- rownames(data)
  
  # Keep only significant genes and sort by absolute log2FoldChange
  sig_large_logFC <- data[!is.na(data$padj) & data$padj <= 0.01, ]
  sig_large_logFC <- sig_large_logFC[order(-abs(sig_large_logFC$log2FoldChange)), ]
  
  # Save all significant genes sorted by absolute log2FoldChange
  write.csv(
    sig_large_logFC,
    file.path(out_dir, paste0(pop, "_padj0.01_sorted_abs_log2FoldChange.csv")),
    row.names = FALSE
  )
  
  # Save top 20 significant genes by absolute log2FoldChange
  top100 <- head(sig_large_logFC, 20)
  
  write.csv(
    top100,
    file.path(out_dir, paste0(pop, "_top20_padj0.01_abs_log2FoldChange.csv")),
    row.names = FALSE
  )
  
  # Print summary
  cat(pop, "\n")
  cat("Significant genes padj <= 0.01:", nrow(sig_large_logFC), "\n")
  cat("Top abs(log2FoldChange):", max(abs(sig_large_logFC$log2FoldChange), na.rm = TRUE), "\n\n")
}