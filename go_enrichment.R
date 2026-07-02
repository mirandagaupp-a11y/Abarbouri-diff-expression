install.packages("tidyverse")
library(tidyverse)
install.packages("readxl")
library(readxl)
library(xlsx)
library(dplyr)
install.packages("xlsx")

# Read in the file
data <- read_tsv("gene_gos2.tsv")

# Extract the relevant columns
gene_id <- data$gene_id
#transcript_id <- data$transcript_id
gene_ontology_blast <- data$gene_ontology_BLASTX
gene_ontology_egg <- data$EggNM.GOs

# Define a function to extract GO terms
extract_GO_terms <- function(text) {
  # Extract all instances of GO terms followed by exactly 7 digits
  go_terms <- regmatches(text, gregexpr("GO:\\d{7}", text))[[1]]
  # Concatenate the extracted GO terms with "|"
  cleaned <- paste(go_terms, collapse = "|")
  return(cleaned)
}

# Apply the function to the gene_ontology columns
cleaned_gene_ontology <- sapply(gene_ontology_blast, extract_GO_terms)
cleaned_egg <- sapply(gene_ontology_egg, extract_GO_terms)

# Create a new data frame with the cleaned columns
cleaned_data <- data.frame(gene_id = gene_id, gene_ontology_BLASTX = cleaned_gene_ontology, EggNM.GOs = cleaned_egg)

# Write the updated data to a new Excel file
write.csv(cleaned_data, "cleaned_gos.csv", row.names = FALSE)

# Merge the two columns into one
merged <- cleaned_data %>%
  mutate(merged_go_ids = paste(gene_ontology_BLASTX, EggNM.GOs, sep = "|")) %>%
  select(-gene_ontology_BLASTX, -EggNM.GOs)


#remove rows that don't have GO terms 
# Remove rows with "." in any column
cleaned_godata <- merged %>%
  filter(!apply(., 1, function(row) any(row == "|")))

write.csv(cleaned_godata, "merged_cleaned_gos.csv", row.names = FALSE)
#### edit the data to make it 1 gene 1 go term per line
# Read the CSV file
data <- read.csv("deduped_merged_gos.csv", header = TRUE)

# Function to split GO terms and create a new data frame
split_go_terms <- function(gene_id, go_terms) {
  go_list <- unlist(strsplit(go_terms, "\\|"))
  data.frame(gene_id = gene_id, go_id = go_list, stringsAsFactors = FALSE)
}

# Apply the function to each row and combine the results
expanded_data <- data %>%
  rowwise() %>%
  do(split_go_terms(.$gene_id, .$merged_go_ids)) %>%
  ungroup()

# Write the expanded data to a new CSV file
write.csv(expanded_data, "exp_go_terms.csv", row.names = FALSE, quote = FALSE)



# Load necessary libraries
library(clusterProfiler)
library(DOSE)
library(GO.db)
library(ggplot2)


# Get all GO terms
go_terms <- keys(GO.db, keytype = "GOID")

# Get GO term names and ontology
go_info <- select(GO.db, keys = go_terms, columns = c("TERM", "ONTOLOGY"), keytype = "GOID")

# Create a data frame with GO terms, their names, and ontology
go_pathway_names <- data.frame(GOterm = go_info$GOID, PathwayName = go_info$TERM, Ontology = go_info$ONTOLOGY, stringsAsFactors = FALSE)

# Write the data frame to a CSV file
write.csv(go_pathway_names, "go_pathway_names.csv", row.names = FALSE)

#read in term2gene file
TERM2GENE <-read.csv("exp_go_terms.csv", header = TRUE)

go_pathway_names <-read.csv("go_pathway_names.csv", header = TRUE)


##### Comparing across up/down and across pops
kegg_degs <-read.csv("kegg_degs.csv", header = TRUE)
#I don't think this works bc of the ont part
all_go <- compareCluster(gene_id~direction+pop,
                         data = kegg_degs, fun = enricher,
                         #pvalueCutoff = 1, 
                         #qvalueCutoff = 1,
                         TERM2GENE = TERM2GENE, 
                         TERM2NAME = go_pathway_names,
                         minGSSize = 1,
                         maxGSSize=10000) 

all_go_df <- as.data.frame(all_go)
write.xlsx(all_go_df, file = 'go_enrichment_comparison.xlsx')

dotplot(ekdegs)

dotplot(
  all_go,
  x = "pop",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "",
  by = "geneRatio",
  size = "geneRatio",
  #includeAll = TRUE,
  label_format = 45,
  facet = "direction"
)

term2gene_bp <- TERM2GENE[TERM2GENE$Pathway %in% go_pathway_names$GOterm[go_pathway_names$Ontology == "BP"], ]
term2gene_mf <- TERM2GENE[TERM2GENE$Pathway %in% go_pathway_names$GOterm[go_pathway_names$Ontology == "MF"], ]
term2gene_cc <- TERM2GENE[TERM2GENE$Pathway %in% go_pathway_names$GOterm[go_pathway_names$Ontology == "CC"], ]

term2name_bp <- go_pathway_names[go_pathway_names$Ontology == "BP", c("GOterm", "PathwayName")]
term2name_mf <- go_pathway_names[go_pathway_names$Ontology == "MF", c("GOterm", "PathwayName")]
term2name_cc <- go_pathway_names[go_pathway_names$Ontology == "CC", c("GOterm", "PathwayName")]

write.table(term2gene_bp, file = "term2gene_bp", sep = "\t", quote = F, row.names = F)
write.table(term2gene_mf, file = "term2gene_mf", sep = "\t", quote = F, row.names = F)
write.table(term2gene_cc, file = "term2gene_cc", sep = "\t", quote = F, row.names = F)

write.table(term2name_bp, file = "term2name_bp", sep = "\t", quote = F, row.names = F)
write.table(term2name_mf, file = "term2name_mf", sep = "\t", quote = F, row.names = F)
write.table(term2name_cc, file = "term2name_cc", sep = "\t", quote = F, row.names = F)

#### Read in term2gene and term2name files
term2gene_bp <- read.table("term2gene_bp", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
term2gene_mf <- read.table("term2gene_mf", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
term2gene_cc <- read.table("term2gene_cc", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)

term2name_bp <- read.table("term2name_bp", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
term2name_mf <- read.table("term2name_mf", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
term2name_cc <- read.table("term2name_cc", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)



# Run compareCluster separately for each
go_bp <- compareCluster(
  gene_id ~ direction + pop,
  data = kegg_degs,
  fun = enricher,
  TERM2GENE = term2gene_bp,
  TERM2NAME = term2name_bp,
  minGSSize = 1,
  maxGSSize = 10000)

go_mf <- compareCluster(
  gene_id ~ direction + pop,
  data = kegg_degs,
  fun = enricher,
  TERM2GENE = term2gene_mf,
  TERM2NAME = term2name_mf,
  minGSSize = 1,
  maxGSSize = 10000)

go_cc <- compareCluster(
  gene_id ~ direction + pop,
  data = kegg_degs,
  fun = enricher,
  TERM2GENE = term2gene_cc,
  TERM2NAME = term2name_cc,
  minGSSize = 1,
  maxGSSize = 10000)

#reorder the clusters so they aren't in alphabetical order when you plot them
# Define your desired order

desired_order <- c("Preble", "Boone", "Scott", "Wilson1", 
                   "Wilson2", "Rutherford")

# Set the factor levels of Cluster to your desired order
go_mf@compareClusterResult$pop <- factor(go_mf@compareClusterResult$pop, levels = desired_order)
go_bp@compareClusterResult$pop <- factor(go_bp@compareClusterResult$pop, levels = desired_order)
go_cc@compareClusterResult$pop <- factor(go_cc@compareClusterResult$pop, levels = desired_order)

desired_order2 <- c("up", "down")

go_mf@compareClusterResult$direction <- factor(go_mf@compareClusterResult$direction, levels = desired_order2)
go_bp@compareClusterResult$direction <- factor(go_bp@compareClusterResult$direction, levels = desired_order2)
go_cc@compareClusterResult$direction <- factor(go_cc@compareClusterResult$direction, levels = desired_order2)


dotplot(go_bp, x = "pop", color = "p.adjust", showCategory = 10,
 font.size = 12, title = "Go Enrichment Biological Processes",
  by = "geneRatio", size = "count", label_format = 60,
  facet = "direction") + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0))

dotplot(
  go_mf,
  x = "pop",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "GO Enrchment Molecular Functions",
  by = "geneRatio",
  size = "count",
  #includeAll = TRUE,
  label_format = 60,
  facet = "direction") + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0))

dotplot(
  go_cc,
  x = "pop",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "Go Enrichment Cellular Components",
  by = "geneRatio",
  size = "count",
  #includeAll = TRUE,
  label_format = 60,
  facet = "direction")  + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0))

#convert to dataframes and export 
bp_df <- as.data.frame(go_bp)
write.csv(bp_df, file = 'go_enrichment_comparison_bp_final.csv')

mf_df <- as.data.frame(go_mf)
write.csv(mf_df, file = 'go_enrichment_comparison_mf_final.csv')

cc_df <- as.data.frame(go_cc)
write.csv(cc_df, file = 'go_enrichment_comparison_cc_final.csv')

#export dotplots as jpegs
p <- dotplot(go_bp, x = "pop", color = "p.adjust", showCategory = 10,
             font.size = 12, title = "Go Enrichment Biological Processes",
             by = "geneRatio", size = "count", label_format = 60,
             facet = "direction") + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        panel.grid = element_blank(),
        axis.title.x = element_blank())

ggsave("GO_BP.jpeg", plot = p,
       device = "jpeg", dpi = 600,
       width = 8, height =8, units = "in")

#export dotplots as jpegs
q <- dotplot(go_mf, x = "pop", color = "p.adjust", showCategory = 10,
             font.size = 12, title = "Go Enrichment Molecular Function",
             by = "geneRatio", size = "count", label_format = 60,
             facet = "direction") + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        panel.grid = element_blank(),
        axis.title.x = element_blank())

ggsave("GO_MF.jpeg", plot = q,
       device = "jpeg", dpi = 600,
       width = 8, height =10, units = "in")

r <- dotplot(go_cc, x = "pop", color = "p.adjust", showCategory = 10,
             font.size = 12, title = "Go Enrichment Cellular Component",
             by = "geneRatio", size = "count", label_format = 60,
             facet = "direction") + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        panel.grid = element_blank(),
        axis.title.x = element_blank())

ggsave("GO_CC.jpeg", plot = r,
       device = "jpeg", dpi = 600,
       width = 8, height =8, units = "in")


### see how many genes have go BP assignments

# all unique input DEG genes
input_genes <- unique(kegg_degs$gene_id)

# genes with GO Biological Process assignments
go_bp_genes <- unique(term2gene_bp$GID)

# genes in your input list that have BP assignments
assigned_bp_genes <- intersect(input_genes, go_bp_genes)

cat("Input genes:", length(input_genes), "\n")
cat("Genes with GO BP assignments:", length(assigned_bp_genes), "\n")
cat("Percent with GO BP assignments:",
    round(100 * length(assigned_bp_genes) / length(input_genes), 1), "%\n")

# this is for the comparison lists for each pop
# I wanted to see if running the enrichment turned up different results running
#it this way and the way above. it shouldn't but we shall see

#read in up/down comparison for each pop
preb_list <-read.csv("kegg_degs_preb.csv", header = TRUE)
boon_list <-read.csv("kegg_degs_boon.csv", header = TRUE)
scot_list <-read.csv("kegg_degs_scot.csv", header = TRUE)
fg_list <-read.csv("kegg_degs_wils_fg.csv", header = TRUE)
rw_list <-read.csv("kegg_degs_wils_rw.csv", header = TRUE)
ruth_list <-read.csv("kegg_degs_ruth.csv", header = TRUE)


ruth_bp <- compareCluster(gene_id~direction+pop,
                       data = ruth_list, fun = enricher,
                       #pvalueCutoff = 1, 
                       #qvalueCutoff = 1,
                       TERM2GENE = term2gene_bp,
                       TERM2NAME = term2name_bp,
                       minGSSize = 1,
                       maxGSSize=1000)

ruth_mf <- compareCluster(gene_id~direction+pop,
                       data = ruth_list, fun = enricher,
                       #pvalueCutoff = 1, 
                       #qvalueCutoff = 1,
                       TERM2GENE = term2gene_mf,
                       TERM2NAME = term2name_mf,
                       minGSSize = 1,
                       maxGSSize=1000)

ruth_cc <- compareCluster(gene_id~direction+pop,
                       data = ruth_list, fun = enricher,
                       #pvalueCutoff = 1, 
                       #qvalueCutoff = 1,
                       TERM2GENE = term2gene_cc,
                       TERM2NAME = term2name_cc,
                       minGSSize = 1,
                       maxGSSize=1000)

rw_bp <- compareCluster(gene_id~direction+pop,
                          data = rw_list, fun = enricher,
                          #pvalueCutoff = 1, 
                          #qvalueCutoff = 1,
                          TERM2GENE = term2gene_bp,
                          TERM2NAME = term2name_bp,
                          minGSSize = 1,
                          maxGSSize=1000)

rw_mf <- compareCluster(gene_id~direction+pop,
                          data = rw_list, fun = enricher,
                          #pvalueCutoff = 1, 
                          #qvalueCutoff = 1,
                          TERM2GENE = term2gene_mf,
                          TERM2NAME = term2name_mf,
                          minGSSize = 1,
                          maxGSSize=1000)

rw_cc <- compareCluster(gene_id~direction+pop,
                          data = rw_list, fun = enricher,
                          #pvalueCutoff = 1, 
                          #qvalueCutoff = 1,
                          TERM2GENE = term2gene_cc,
                          TERM2NAME = term2name_cc,
                          minGSSize = 1,
                          maxGSSize=1000)

fg_bp <- compareCluster(gene_id~direction+pop,
                          data = fg_list, fun = enricher,
                          #pvalueCutoff = 1, 
                          #qvalueCutoff = 1,
                          TERM2GENE = term2gene_bp,
                          TERM2NAME = term2name_bp,
                          minGSSize = 1,
                          maxGSSize=1000)

fg_mf <- compareCluster(gene_id~direction+pop,
                          data = fg_list, fun = enricher,
                          #pvalueCutoff = 1, 
                          #qvalueCutoff = 1,
                          TERM2GENE = term2gene_mf,
                          TERM2NAME = term2name_mf,
                          minGSSize = 1,
                          maxGSSize=1000)

fg_cc <- compareCluster(gene_id~direction+pop,
                          data = fg_list, fun = enricher,
                          #pvalueCutoff = 1, 
                          #qvalueCutoff = 1,
                          TERM2GENE = term2gene_cc,
                          TERM2NAME = term2name_cc,
                          minGSSize = 1,
                          maxGSSize=1000)

scot_bp <- compareCluster(gene_id~direction+pop,
                          data = scot_list, fun = enricher,
                          #pvalueCutoff = 1, 
                          #qvalueCutoff = 1,
                          TERM2GENE = term2gene_bp,
                          TERM2NAME = term2name_bp,
                          minGSSize = 1,
                          maxGSSize=1000)

scot_mf <- compareCluster(gene_id~direction+pop,
                          data = scot_list, fun = enricher,
                          #pvalueCutoff = 1, 
                          #qvalueCutoff = 1,
                          TERM2GENE = term2gene_mf,
                          TERM2NAME = term2name_mf,
                          minGSSize = 1,
                          maxGSSize=1000)

scot_cc <- compareCluster(gene_id~direction+pop,
                          data = scot_list, fun = enricher,
                          #pvalueCutoff = 1, 
                          #qvalueCutoff = 1,
                          TERM2GENE = term2gene_cc,
                          TERM2NAME = term2name_cc,
                          minGSSize = 1,
                          maxGSSize=1000)

boon_bp <- compareCluster(gene_id~direction+pop,
                          data = boon_list, fun = enricher,
                          #pvalueCutoff = 1, 
                          #qvalueCutoff = 1,
                          TERM2GENE = term2gene_bp,
                          TERM2NAME = term2name_bp,
                          minGSSize = 1,
                          maxGSSize=1000)

boon_mf <- compareCluster(gene_id~direction+pop,
                          data = boon_list, fun = enricher,
                          #pvalueCutoff = 1, 
                          #qvalueCutoff = 1,
                          TERM2GENE = term2gene_mf,
                          TERM2NAME = term2name_mf,
                          minGSSize = 1,
                          maxGSSize=1000)

boon_cc <- compareCluster(gene_id~direction+pop,
                          data = boon_list, fun = enricher,
                          #pvalueCutoff = 1, 
                          #qvalueCutoff = 1,
                          TERM2GENE = term2gene_cc,
                          TERM2NAME = term2name_cc,
                          minGSSize = 1,
                          maxGSSize=1000)

preb_bp <- compareCluster(gene_id~direction+pop,
                          data = preb_list, fun = enricher,
                          #pvalueCutoff = 1, 
                          #qvalueCutoff = 1,
                          TERM2GENE = term2gene_bp,
                          TERM2NAME = term2name_bp,
                          minGSSize = 1,
                          maxGSSize=1000)

pren_mf <- compareCluster(gene_id~direction+pop,
                          data = preb_list, fun = enricher,
                          #pvalueCutoff = 1, 
                          #qvalueCutoff = 1,
                          TERM2GENE = term2gene_mf,
                          TERM2NAME = term2name_mf,
                          minGSSize = 1,
                          maxGSSize=1000)

preb_cc <- compareCluster(gene_id~direction+pop,
                          data = preb_list, fun = enricher,
                          #pvalueCutoff = 1, 
                          #qvalueCutoff = 1,
                          TERM2GENE = term2gene_cc,
                          TERM2NAME = term2name_cc,
                          minGSSize = 1,
                          maxGSSize=1000)




dotplot(
  ruth_bp,
  x = "direction",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "Ruth BP",
  by = "geneRatio",
  size = "RichFactor",
  #includeAll = TRUE,
  label_format = 45,
  #facet = "pop"
)

dotplot(
  ruth_mf,
  x = "direction",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "Ruth MF",
  by = "geneRatio",
  size = "RichFactor",
  #includeAll = TRUE,
  label_format = 45,
  #facet = "pop"
)

dotplot(
  ruth_cc,
  x = "direction",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "Ruth CC",
  by = "geneRatio",
  size = "RichFactor",
  #includeAll = TRUE,
  label_format = 45,
  #facet = "pop"
)


dotplot(
  rw_bp,
  x = "direction",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "Rw BP",
  by = "geneRatio",
  size = "RichFactor",
  #includeAll = TRUE,
  label_format = 45,
  #facet = "pop"
)

dotplot(
  rw_mf,
  x = "direction",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "Rw MF",
  by = "geneRatio",
  size = "RichFactor",
  #includeAll = TRUE,
  label_format = 45,
  #facet = "pop"
)

dotplot(
  rw_cc,
  x = "direction",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "Rw CC",
  by = "geneRatio",
  size = "RichFactor",
  #includeAll = TRUE,
  label_format = 45,
  #facet = "pop"
)


dotplot(
  fg_bp,
  x = "direction",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "fg BP",
  by = "geneRatio",
  size = "RichFactor",
  #includeAll = TRUE,
  label_format = 45,
  #facet = "pop"
)

dotplot(
  fg_mf,
  x = "direction",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "fg MF",
  by = "geneRatio",
  size = "RichFactor",
  #includeAll = TRUE,
  label_format = 45,
  #facet = "pop"
)

dotplot(
  fg_cc,
  x = "direction",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "fg CC",
  by = "geneRatio",
  size = "RichFactor",
  #includeAll = TRUE,
  label_format = 45,
  #facet = "pop"
)


dotplot(
  scot_bp,
  x = "direction",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "scot BP",
  by = "geneRatio",
  size = "RichFactor",
  #includeAll = TRUE,
  label_format = 45,
  #facet = "pop"
)

dotplot(
  scot_mf,
  x = "direction",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "scot MF",
  by = "geneRatio",
  size = "RichFactor",
  #includeAll = TRUE,
  label_format = 45,
  #facet = "pop"
)

dotplot(
  scot_cc,
  x = "direction",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "scot CC",
  by = "geneRatio",
  size = "RichFactor",
  #includeAll = TRUE,
  label_format = 45,
  #facet = "pop"
)


dotplot(
  boon_mf,
  x = "direction",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "boon MF",
  by = "geneRatio",
  size = "RichFactor",
  #includeAll = TRUE,
  label_format = 45,
  #facet = "pop"
)

dotplot(
  boon_cc,
  x = "direction",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "boon CC",
  by = "geneRatio",
  size = "RichFactor",
  #includeAll = TRUE,
  label_format = 45,
  #facet = "pop"
)



dotplot(
  pren_mf,
  x = "direction",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "preb MF",
  by = "geneRatio",
  size = "RichFactor",
  #includeAll = TRUE,
  label_format = 45,
  #facet = "pop"
)

dotplot(
  preb_cc,
  x = "direction",
  color = "p.adjust",
  showCategory = 10,
  #split = "pop",
  font.size = 12,
  title = "preb CC",
  by = "geneRatio",
  size = "RichFactor",
  #includeAll = TRUE,
  label_format = 45,
  #facet = "pop"
)

ruth_bp_df <- as.data.frame(ruth_bp)
write.csv(ruth_bp_df, file = 'ruth_bp_comparison_enrichment.csv')

ruth_mf_df <- as.data.frame(ruth_mf)
write.csv(ruth_mf_df, file = 'ruth_mf_comparison_enrichment.csv')

ruth_cc_df <- as.data.frame(ruth_cc)
write.csv(ruth_cc_df, file = 'ruth_cc_comparison_enrichment.csv')

rw_bp_df <- as.data.frame(rw_bp)
write.csv(rw_bp_df, file = 'rw_bp_comparison_enrichment.csv')

rw_mf_df <- as.data.frame(rw_mf)
write.csv(rw_mf_df, file = 'rw_mf_comparison_enrichment.csv')

rw_cc_df <- as.data.frame(rw_cc)
write.csv(rw_cc_df, file = 'rw_cc_comparison_enrichment.csv')

fg_bp_df <- as.data.frame(fg_bp)
write.csv(fg_bp_df, file = 'fg_bp_comparison_enrichment.csv')

fg_mf_df <- as.data.frame(fg_mf)
write.csv(fg_mf_df, file = 'fg_mf_comparison_enrichment.csv')

fg_cc_df <- as.data.frame(fg_cc)
write.csv(fg_cc_df, file = 'fg_cc_comparison_enrichment.csv')

scot_bp_df <- as.data.frame(scot_bp)
write.csv(scot_bp_df, file = 'scot_bp_comparison_enrichment.csv')

scot_mf_df <- as.data.frame(scot_mf)
write.csv(scot_mf_df, file = 'scot_mf_comparison_enrichment.csv')

scot_cc_df <- as.data.frame(scot_cc)
write.csv(scot_cc_df, file = 'scot_cc_comparison_enrichment.csv')

boon_mf_df <- as.data.frame(boon_mf)
write.csv(boon_mf_df, file = 'boon_mf_comparison_enrichment.csv')

boon_cc_df <- as.data.frame(boon_cc)
write.csv(boon_cc_df, file = 'boon_cc_comparison_enrichment.csv')

preb_mf_df <- as.data.frame(pren_mf)
write.csv(preb_mf_df, file = 'preb_mf_comparison_enrichment.csv')

preb_cc_df <- as.data.frame(preb_cc)
write.csv(preb_cc_df, file = 'preb_cc_comparison_enrichment.csv')
