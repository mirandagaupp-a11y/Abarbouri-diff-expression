#create a gene counts matrix istead of a transcript counts matrix and run deseq

# Load necessary libraries
library(tximport)
library(readr)
library(DESeq2)

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

# Round the NumReads to the nearest integer
rounded_counts <- round(txi$counts)

# Ensure the columns of the counts matrix match the sample names in the metadata
rounded_counts <- rounded_counts[, metadata$sampleName]

# Save the matrix of gene counts to a file
write.csv(rounded_counts, file = "/work/projects/5-39268/mgaupp42/diff_express/salmon/gene_matrix1.csv", quote=FALSE)

samp <- as.data.frame(read.csv("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/population/metadata.csv", head = T, row.names=1))

# Construct a DESeqDataSet from the txi object and sample information in metadata
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samp,
                                   design = ~ Population)
                                                                      
## Run the thing
dds <- DESeq(ddsTxi)

names <- print(resultsNames(dds))
write.csv(names, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/population/result_names1.csv")

#### get the results ####
results <- results(dds)

write.csv(results,"/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/population/deseq_results.csv", quote=FALSE)

save(dds, file="/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/population/deseq_object.RData")

library(ggplot2)

## Perform PCA
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Population", "Temperature"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

## Ensure Temperature is treated as a factor
pcaData$Temperature <- as.factor(pcaData$Temperature)

## Plot PCA
p <- ggplot(pcaData, aes(PC1, PC2, color=Population, shape=Temperature)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  coord_fixed(ratio = 2) 

## Save the plot
#ggsave("/work/projects/5-39268/mgaupp42/deseq_post_evigene/population/pca_plot1.png", plot=p)
ggsave("/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/population/pca_plot1.png", plot=p, width=5, height=5)








