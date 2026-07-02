library(DESeq2)
library(tximport)
library(vegan)

# Load subsetted tximport object
txi <- readRDS(
"/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/permanova/genelevel_permanova_subset.rds"
)

# Metadata
samp <- read.csv(
"/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/metadata.csv",
header = TRUE,
row.names = 1
)

# Match order
txi$counts <- txi$counts[, rownames(samp)]

# Construct DESeq object
dds <- DESeqDataSetFromTximport(
txi,
colData = samp,
design = ~ 1
)

# Estimate size factors only
dds <- estimateSizeFactors(dds)

# Variance-stabilized expression
vsd <- vst(dds, blind = TRUE)

expr <- assay(vsd)

# Factors
samp$Population <- factor(samp$Population)
samp$Temperature <- factor(samp$Temperature)
samp$Clutch <- factor(samp$Clutch)
samp$RangeUnit <- factor(samp$RangeUnit)

# Distance
sample_dist <- dist(t(expr), method = "euclidean")

# PERMANOVA
permanova_res <- adonis2(
sample_dist ~ RangeUnit + Population + Temperature +
Population:Temperature + Clutch,
data = samp,
by = "terms",
permutations = 999
)

write.csv(
as.data.frame(permanova_res),
"/work/projects/5-39268/mgaupp42/diff_express/deseq/gene_level/length_correction/temp/permanova/permanova_results.csv",
quote = FALSE
)