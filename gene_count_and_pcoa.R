#create a gene counts matrix istead of a transcript counts matrix and run deseq

# Load necessary libraries
library(tximport)
library(readr)
library(DESeq2)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(vegan)


# Define the directory containing the Salmon quant files
quant_dir <- "/work/projects/5-39268/mgaupp42/diff_express/salmon/"

# List all .sf files recursively
files <- list.files(quant_dir, pattern = "\\.sf$", full.names = TRUE, recursive = TRUE)

# Create a sample table
# Assuming you have a metadata file with sample names and conditions
metadata <- read.csv("/work/projects/5-39268/mgaupp42/diff_express/salmon/metadata.csv") # Adjust the path to your metadata file

# Ensure the metadata file has columns: sampleName, condition
# sampleName should match the names of the .sf files (without the .sf extension)

# Create a named vector of files
names(files) <- gsub("\\.sf$", "", basename(files))

# Read the transcript-to-gene mapping file
# Adjust the path and delimiter as needed
tx2gene <- read.table("/work/projects/5-39268/mgaupp42/diff_express/salmon/gene_to_trans_map.txt", header = FALSE, sep = "\t") # Use sep = "\t" for tab-delimited files

# Swap the columns to ensure the first column is transcripts and the second is genes
tx2gene <- tx2gene[, c(2, 1)]

# Assign column names manually
colnames(tx2gene) <- c("transcript", "gene")

# Import Salmon quant files using tximport with tx2gene
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# Save txi for reuse
saveRDS(txi, file = "/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/txi_genelevel.rds")

#### Deseq for PCoA ####
txi <- readRDS("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/txi_genelevel.rds")
samp <- as.data.frame(read.csv("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/metadata.csv", head = T, row.names=1))

# design is Temp because if using Population design, we have to set a reference population
ddsTxi_temp <- DESeqDataSetFromTximport(txi,
                                   colData = samp,
                                   design = ~ Temperature)

ddsTxi_temp$Temperature <- relevel(ddsTxi_temp$Temperature, ref = "10C")

dds_txi <- DESeq(ddsTxi_temp)
save(dds_txi, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/population/temp_deseq_object.RData")

# Variance stabilizing transformation
vsd <- vst(dds_txi, blind = TRUE)
vst_mat <- assay(vsd)

samp$Sample <- rownames(samp)

# Euclidean distance among samples
sample_dist <- dist(t(vst_mat), method = "euclidean")

# PCoA
pcoa_sample <- capscale(sample_dist ~ 1)
pcoa_scores <- scores(pcoa_sample)$sites

eigvals <- eigenvals(pcoa_sample)
var_explained <- eigvals / sum(eigvals) * 100

# Prepare plotting dataframe
pcoa_df <- data.frame(
  Sample = rownames(pcoa_scores),
  PCoA1 = pcoa_scores[,1],
  PCoA2 = pcoa_scores[,2]
) %>%
  left_join(samp, by = "Sample")

# Rename populations for plotting
pcoa_df <- pcoa_df %>% mutate(Population = recode(Population, "BOON" = "Boone", "PREB" = "Preble", "RUTH" = "Rutherford", "SCOT" = "Scott", "WILS_FG" = "Wilson1", "WILS_RW" = "Wilson2"))

# Set population order in legend
pcoa_df$Population <- factor(pcoa_df$Population, levels = c("Preble", "Boone", "Scott", "Wilson1", "Wilson2", "Rutherford"))

# Color matching 
population_colors <- c("Preble" = "#E41A1C","Boone" = "#377EB8", "Scott" = "#4DAF4A", "Wilson1" = "#984EA3", "Wilson2" = "#FF7F00", "Rutherford" = "#999999")

# Export as JPEG
jpeg("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/population/sample_pcoa_pop_TEST.jpeg", width = 2400, height = 2000, res = 300)

ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Population, shape = Temperature)) +
  geom_point(size = 4) + scale_color_manual(values = population_colors) +
 scale_shape_manual(values = c("10C" = 16, "20C" = 17)) +
  labs(x = paste0("PCoA1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_bw(base_size = 12)

dev.off()
