# Load necessary libraries
library(tximport)

# Load tximport object
txi <- readRDS(
  "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/txi_genelevel.rds"
)

# Define the individuals you want to include
   
individuals_to_include <- c("WILS_RW_1_3", "WILS_RW_1_4", "WILS_RW_1_6", "WILS_RW_2_3", "WILS_RW_2_4", "WILS_RW_2_5", "WILS_RW_1_7", "WILS_RW_1_8", "WILS_RW_1_9",
    "WILS_RW_2_12", "WILS_RW_2_13", "WILS_RW_2_14")

# Subset tximport object
txi_subset <- txi

txi_subset$counts <-
  txi$counts[, individuals_to_include, drop = FALSE]

txi_subset$abundance <-
  txi$abundance[, individuals_to_include, drop = FALSE]

txi_subset$length <-
  txi$length[, individuals_to_include, drop = FALSE]

# Save subsetted tximport object
saveRDS(
  txi_subset,
  "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/WILS_RW_txi_subset.rds"
)

### PREB
individuals_to_include <- c("PREB_HW_1_3", "PREB_HW_1_5", "PREB_HW_1_7", "PREB_HW_1_8", "PREB_HW_1_9", "PREB_HW_2_3", "PREB_HW_2_5",
    "PREB_HW_2_2", "PREB_HW_2_8", "PREB_HW_2_9", "PREB_HW_2_10")

# Subset tximport object
txi_subset <- txi

txi_subset$counts <-
  txi$counts[, individuals_to_include, drop = FALSE]

txi_subset$abundance <-
  txi$abundance[, individuals_to_include, drop = FALSE]

txi_subset$length <-
  txi$length[, individuals_to_include, drop = FALSE]

# Save subsetted tximport object
saveRDS(
  txi_subset,
  "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/PREB_txi_subset.rds"
)

individuals_to_include <- c("BOON_WALT_2_2", "BOON_WALT_2_4", "BOON_WALT_2_5", "BOON_WALT_2_8", "BOON_WALT_2_10", "BOON_WALT_2_11", "BOON_WALT_3_3", "BOON_WALT_3_5",
    "BOON_WALT_3_6", "BOON_WALT_3_7", "BOON_WALT_3_8", "BOON_WALT_3_10")
    
# Subset tximport object
txi_subset <- txi

txi_subset$counts <-
  txi$counts[, individuals_to_include, drop = FALSE]

txi_subset$abundance <-
  txi$abundance[, individuals_to_include, drop = FALSE]

txi_subset$length <-
  txi$length[, individuals_to_include, drop = FALSE]

# Save subsetted tximport object
saveRDS(
  txi_subset,
  "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/BOON_txi_subset.rds"
)

individuals_to_include <- c("SCOT_JAM_1_1", "SCOT_JAM_1_2", "SCOT_JAM_1_5", "SCOT_JAM_1_7", "SCOT_JAM_1_8", "SCOT_JAM_1_11", "SCOT_JAM_2_2", "SCOT_JAM_2_5", "SCOT_JAM_2_6",
    "SCOT_JAM_2_8", "SCOT_JAM_2_10", "SCOT_JAM_2_11")

# Subset tximport object
txi_subset <- txi

txi_subset$counts <-
  txi$counts[, individuals_to_include, drop = FALSE]

txi_subset$abundance <-
  txi$abundance[, individuals_to_include, drop = FALSE]

txi_subset$length <-
  txi$length[, individuals_to_include, drop = FALSE]

# Save subsetted tximport object
saveRDS(
  txi_subset,
  "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/SCOT_txi_subset.rds"
)
   
individuals_to_include <- c("WILS_FG_1_3", "WILS_FG_1_8", "WILS_FG_1_9", "WILS_FG_2_9", "WILS_FG_2_10", "WILS_FG_2_15", "WILS_FG_3_10", "WILS_FG_3_11", "WILS_FG_3_14",
    "WILS_FG_1_13", "WILS_FG_1_15", "WILS_FG_1_17", "WILS_FG_2_1", "WILS_FG_2_2", "WILS_FG_2_3", "WILS_FG_3_2", "WILS_FG_3_3", "WILS_FG_3_4")

# Subset tximport object
txi_subset <- txi

txi_subset$counts <-
  txi$counts[, individuals_to_include, drop = FALSE]

txi_subset$abundance <-
  txi$abundance[, individuals_to_include, drop = FALSE]

txi_subset$length <-
  txi$length[, individuals_to_include, drop = FALSE]

# Save subsetted tximport object
saveRDS(
  txi_subset,
  "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/WILS_FG_txi_subset.rds"
)

individuals_to_include <- c("RUTH_LCC_2_2", "RUTH_LCC_2_3", "RUTH_LCC_2_4", "RUTH_LCC_2_8", "RUTH_LCC_1_3", "RUTH_LCC_1_5", "RUTH_LCC_1_8", "RUTH_LCC_1_10", "RUTH_LCC_1_12",
    "RUTH_LCC_1_14", "RUTH_LCC_2_10", "RUTH_LCC_2_11", "RUTH_LCC_2_12", "RUTH_LCC_2_13")

# Subset tximport object
txi_subset <- txi

txi_subset$counts <-
  txi$counts[, individuals_to_include, drop = FALSE]

txi_subset$abundance <-
  txi$abundance[, individuals_to_include, drop = FALSE]

txi_subset$length <-
  txi$length[, individuals_to_include, drop = FALSE]

# Save subsetted tximport object
saveRDS(
  txi_subset,
  "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/RUTH_txi_subset.rds"
)
