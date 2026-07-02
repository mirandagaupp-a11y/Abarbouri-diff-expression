#!/usr/bin/env Rscript

# ============================
# User-defined file paths
# ============================
genes_txt   <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/wollenberg_genes.txt"         # text file with one gene ID per line
input_csv   <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/PREB_20vs10_deseq_results.csv" # CSV file with gene IDs in the first column
output_tsv  <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/PREB_wb_matched_genes.tsv" # output file
missing_out <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/PREB_genes_MISSING.txt" # missing IDs file

# ============================
# Libraries
# ============================
suppressWarnings(suppressMessages({
  library(data.table)
}))

# --- 1) Read gene list ---
genes <- unique(trimws(readLines(genes_txt, warn = FALSE)))
genes <- genes[nzchar(genes)]
if (length(genes) == 0) stop("No gene IDs found in the txt file.")

# --- 2) Read CSV (fast, memory-friendly; handles .gz too) ---
dt <- fread(input_csv)
if (is.na(names(dt)[1]) || names(dt)[1] == "" || names(dt)[1] == "V1") {
  setnames(dt, 1, "gene_id")
} else {
  setnames(dt, 1, "gene_id")
}

# --- 3) Filter to requested genes ---
dt_sub <- dt[gene_id %in% genes]

# --- 4) Reorder rows to match the order of the input gene list ---
dt_sub[, gene_order := match(gene_id, genes)]
setorder(dt_sub, gene_order)
dt_sub[, gene_order := NULL]

# --- 5) Write outputs ---
fwrite(dt_sub, output_tsv, sep = "\t", quote = FALSE, na = "NA")

not_found <- setdiff(genes, dt_sub$gene_id)
if (length(not_found) > 0) {
  writeLines(not_found, con = missing_out)
  message(sprintf("Wrote %d matched rows to %s", nrow(dt_sub), output_tsv))
  message(sprintf("Wrote %d missing gene IDs to %s", length(not_found), missing_out))
} else {
  if (file.exists(missing_out)) unlink(missing_out)
  message(sprintf("Wrote %d matched rows to %s; all gene IDs were found.", nrow(dt_sub), output_tsv))
}

################## Boone
genes_txt   <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/wollenberg_genes.txt"         # text file with one gene ID per line
input_csv   <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/BOON_20vs10_deseq_results.csv" # CSV file with gene IDs in the first column
output_tsv  <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/BOON_wb_matched_genes.tsv" # output file
missing_out <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/BOON_genes_MISSING.txt" # missing IDs file

# ============================
# Libraries
# ============================
suppressWarnings(suppressMessages({
  library(data.table)
}))

# --- 1) Read gene list ---
genes <- unique(trimws(readLines(genes_txt, warn = FALSE)))
genes <- genes[nzchar(genes)]
if (length(genes) == 0) stop("No gene IDs found in the txt file.")

# --- 2) Read CSV (fast, memory-friendly; handles .gz too) ---
dt <- fread(input_csv)
if (is.na(names(dt)[1]) || names(dt)[1] == "" || names(dt)[1] == "V1") {
  setnames(dt, 1, "gene_id")
} else {
  setnames(dt, 1, "gene_id")
}

# --- 3) Filter to requested genes ---
dt_sub <- dt[gene_id %in% genes]

# --- 4) Reorder rows to match the order of the input gene list ---
dt_sub[, gene_order := match(gene_id, genes)]
setorder(dt_sub, gene_order)
dt_sub[, gene_order := NULL]

# --- 5) Write outputs ---
fwrite(dt_sub, output_tsv, sep = "\t", quote = FALSE, na = "NA")

not_found <- setdiff(genes, dt_sub$gene_id)
if (length(not_found) > 0) {
  writeLines(not_found, con = missing_out)
  message(sprintf("Wrote %d matched rows to %s", nrow(dt_sub), output_tsv))
  message(sprintf("Wrote %d missing gene IDs to %s", length(not_found), missing_out))
} else {
  if (file.exists(missing_out)) unlink(missing_out)
  message(sprintf("Wrote %d matched rows to %s; all gene IDs were found.", nrow(dt_sub), output_tsv))
}

############## Scott
genes_txt   <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/wollenberg_genes.txt"         # text file with one gene ID per line
input_csv   <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/SCOT_20vs10_deseq_results.csv" # CSV file with gene IDs in the first column
output_tsv  <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/SCOT_wb_matched_genes.tsv" # output file
missing_out <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/SCOT_genes_MISSING.txt" # missing IDs file

# ============================
# Libraries
# ============================
suppressWarnings(suppressMessages({
  library(data.table)
}))

# --- 1) Read gene list ---
genes <- unique(trimws(readLines(genes_txt, warn = FALSE)))
genes <- genes[nzchar(genes)]
if (length(genes) == 0) stop("No gene IDs found in the txt file.")

# --- 2) Read CSV (fast, memory-friendly; handles .gz too) ---
dt <- fread(input_csv)
if (is.na(names(dt)[1]) || names(dt)[1] == "" || names(dt)[1] == "V1") {
  setnames(dt, 1, "gene_id")
} else {
  setnames(dt, 1, "gene_id")
}

# --- 3) Filter to requested genes ---
dt_sub <- dt[gene_id %in% genes]

# --- 4) Reorder rows to match the order of the input gene list ---
dt_sub[, gene_order := match(gene_id, genes)]
setorder(dt_sub, gene_order)
dt_sub[, gene_order := NULL]

# --- 5) Write outputs ---
fwrite(dt_sub, output_tsv, sep = "\t", quote = FALSE, na = "NA")

not_found <- setdiff(genes, dt_sub$gene_id)
if (length(not_found) > 0) {
  writeLines(not_found, con = missing_out)
  message(sprintf("Wrote %d matched rows to %s", nrow(dt_sub), output_tsv))
  message(sprintf("Wrote %d missing gene IDs to %s", length(not_found), missing_out))
} else {
  if (file.exists(missing_out)) unlink(missing_out)
  message(sprintf("Wrote %d matched rows to %s; all gene IDs were found.", nrow(dt_sub), output_tsv))
}


######################## WILS FG
genes_txt   <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/wollenberg_genes.txt"         # text file with one gene ID per line
input_csv   <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/WILS_FG_20vs10_deseq_results.csv" # CSV file with gene IDs in the first column
output_tsv  <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/WILS_FG_wb_matched_genes.tsv" # output file
missing_out <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/WILS_FG_genes_MISSING.txt" # missing IDs file

# ============================
# Libraries
# ============================
suppressWarnings(suppressMessages({
  library(data.table)
}))

# --- 1) Read gene list ---
genes <- unique(trimws(readLines(genes_txt, warn = FALSE)))
genes <- genes[nzchar(genes)]
if (length(genes) == 0) stop("No gene IDs found in the txt file.")

# --- 2) Read CSV (fast, memory-friendly; handles .gz too) ---
dt <- fread(input_csv)
if (is.na(names(dt)[1]) || names(dt)[1] == "" || names(dt)[1] == "V1") {
  setnames(dt, 1, "gene_id")
} else {
  setnames(dt, 1, "gene_id")
}

# --- 3) Filter to requested genes ---
dt_sub <- dt[gene_id %in% genes]

# --- 4) Reorder rows to match the order of the input gene list ---
dt_sub[, gene_order := match(gene_id, genes)]
setorder(dt_sub, gene_order)
dt_sub[, gene_order := NULL]

# --- 5) Write outputs ---
fwrite(dt_sub, output_tsv, sep = "\t", quote = FALSE, na = "NA")

not_found <- setdiff(genes, dt_sub$gene_id)
if (length(not_found) > 0) {
  writeLines(not_found, con = missing_out)
  message(sprintf("Wrote %d matched rows to %s", nrow(dt_sub), output_tsv))
  message(sprintf("Wrote %d missing gene IDs to %s", length(not_found), missing_out))
} else {
  if (file.exists(missing_out)) unlink(missing_out)
  message(sprintf("Wrote %d matched rows to %s; all gene IDs were found.", nrow(dt_sub), output_tsv))
}


##################### WILS_RW
genes_txt   <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/wollenberg_genes.txt"         # text file with one gene ID per line
input_csv   <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/WILS_RW_20vs10_deseq_results.csv" # CSV file with gene IDs in the first column
output_tsv  <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/WILS_RW_wb_matched_genes.tsv" # output file
missing_out <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/WILS_RW_genes_MISSING.txt" # missing IDs file

# ============================
# Libraries
# ============================
suppressWarnings(suppressMessages({
  library(data.table)
}))

# --- 1) Read gene list ---
genes <- unique(trimws(readLines(genes_txt, warn = FALSE)))
genes <- genes[nzchar(genes)]
if (length(genes) == 0) stop("No gene IDs found in the txt file.")

# --- 2) Read CSV (fast, memory-friendly; handles .gz too) ---
dt <- fread(input_csv)
if (is.na(names(dt)[1]) || names(dt)[1] == "" || names(dt)[1] == "V1") {
  setnames(dt, 1, "gene_id")
} else {
  setnames(dt, 1, "gene_id")
}

# --- 3) Filter to requested genes ---
dt_sub <- dt[gene_id %in% genes]

# --- 4) Reorder rows to match the order of the input gene list ---
dt_sub[, gene_order := match(gene_id, genes)]
setorder(dt_sub, gene_order)
dt_sub[, gene_order := NULL]

# --- 5) Write outputs ---
fwrite(dt_sub, output_tsv, sep = "\t", quote = FALSE, na = "NA")

not_found <- setdiff(genes, dt_sub$gene_id)
if (length(not_found) > 0) {
  writeLines(not_found, con = missing_out)
  message(sprintf("Wrote %d matched rows to %s", nrow(dt_sub), output_tsv))
  message(sprintf("Wrote %d missing gene IDs to %s", length(not_found), missing_out))
} else {
  if (file.exists(missing_out)) unlink(missing_out)
  message(sprintf("Wrote %d matched rows to %s; all gene IDs were found.", nrow(dt_sub), output_tsv))
}


#################### RUTH
genes_txt   <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/wollenberg_genes.txt"         # text file with one gene ID per line
input_csv   <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/results/RUTH_20vs10_deseq_results.csv" # CSV file with gene IDs in the first column
output_tsv  <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/RUTH_wb_matched_genes.tsv" # output file
missing_out <- "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/wollenberg/RUTH_genes_MISSING.txt" # missing IDs file

# ============================
# Libraries
# ============================
suppressWarnings(suppressMessages({
  library(data.table)
}))

# --- 1) Read gene list ---
genes <- unique(trimws(readLines(genes_txt, warn = FALSE)))
genes <- genes[nzchar(genes)]
if (length(genes) == 0) stop("No gene IDs found in the txt file.")

# --- 2) Read CSV (fast, memory-friendly; handles .gz too) ---
dt <- fread(input_csv)
if (is.na(names(dt)[1]) || names(dt)[1] == "" || names(dt)[1] == "V1") {
  setnames(dt, 1, "gene_id")
} else {
  setnames(dt, 1, "gene_id")
}

# --- 3) Filter to requested genes ---
dt_sub <- dt[gene_id %in% genes]

# --- 4) Reorder rows to match the order of the input gene list ---
dt_sub[, gene_order := match(gene_id, genes)]
setorder(dt_sub, gene_order)
dt_sub[, gene_order := NULL]

# --- 5) Write outputs ---
fwrite(dt_sub, output_tsv, sep = "\t", quote = FALSE, na = "NA")

not_found <- setdiff(genes, dt_sub$gene_id)
if (length(not_found) > 0) {
  writeLines(not_found, con = missing_out)
  message(sprintf("Wrote %d matched rows to %s", nrow(dt_sub), output_tsv))
  message(sprintf("Wrote %d missing gene IDs to %s", length(not_found), missing_out))
} else {
  if (file.exists(missing_out)) unlink(missing_out)
  message(sprintf("Wrote %d matched rows to %s; all gene IDs were found.", nrow(dt_sub), output_tsv))
}

