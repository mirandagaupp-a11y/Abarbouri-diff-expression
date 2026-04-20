
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(openxlsx)

# Load KEGG annotation tables
pathway2gene <- read.csv("kegg_pathway2gene.csv", header = TRUE)
pathway2name <- read.table("kegg_pathway2name_gene", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Human disease pathway codes — add 'ko' prefix
human_disease_pathways <- c(
  "05200", "05202", "05206", "05205", "05204", "05207", "05208", "05203", "05230", "05231", "05235",
  "05210", "05212", "05225", "05226", "05214", "05216", "05221", "05220", "05217", "05218", "05211",
  "05219", "05215", "05213", "05224", "05222", "05223", "05166", "05170", "05161", "05160", "05171",
  "05164", "05162", "05168", "05163", "05167", "05169", "05165", "05110", "05120", "05130", "05132",
  "05131", "05135", "05133", "05134", "05150", "05152", "05100", "05146", "05144", "05145", "05140",
  "05142", "05143", "05310", "05322", "05323", "05320", "05321", "05330", "05332", "05340", "05010",
  "05012", "05014", "05016", "05017", "05020", "05022", "05030", "05031", "05032", "05033", "05034",
  "05417", "05418", "05410", "05412", "05414", "05415", "05416", "04930", "04940", "04950", "04936",
  "04932", "04931", "04933", "04934", "01501", "01502", "01503", "01521", "01524", "01523", "01522"
)

human_disease_pathways <- paste0("ko", human_disease_pathways)

# Filter out human disease pathways from pathway2name
pathway2name_filtered2 <- pathway2name_filtered %>%
  filter(!Pathway %in% human_disease_pathways)

# Filter pathway2gene to match
pathway2gene_filtered2 <- pathway2gene_filtered %>%
  filter(Pathway %in% pathway2name_filtered2$Pathway)

# Save filtered versions if needed
write.table(pathway2name_filtered2, file = "kegg_pathway2name_wo_hd.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pathway2gene_filtered2, file = "kegg_pathway2gene_wo_hd.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Read in filtered files if starting from here 
pathway2gene_filtered <- read.table("kegg_pathway2gene_filtered.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pathway2name_filtered <- read.table("kegg_pathway2name_filtered.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

deg_combined <- read.csv("./deg_lists/kegg_degs.csv")

ekdegs <- compareCluster(gene_id ~ direction + pop,
                         data = deg_combined,
                         fun = enricher,
                         TERM2GENE = pathway2gene_filtered2,
                         TERM2NAME = pathway2name_filtered2,
                         minGSSize = 1,
                         maxGSSize = 1000)

desired_order <- c("Preble", "Boone", "Scott", "Wils_FG", 
                   "Wils_RW", "Rutherford")

ekdegs@compareClusterResult$pop <- factor(ekdegs@compareClusterResult$pop, levels = desired_order)

desired_order2 <- c("up", "down")

ekdegs@compareClusterResult$direction <- factor(ekdegs@compareClusterResult$direction, levels = desired_order2)


write.xlsx(as.data.frame(ekdegs), file = "kegg_comparison_enrichment_wo_hd.xlsx")


dotplot(ekdegs, x = "pop", color = "p.adjust", showCategory = 8,
        title = "KEGG Comparison", by = "geneRatio", size = "Count",
        facet = "direction", label_format = 60, font.size = 12)

dotplot(ekdegs, x = "pop", color = "p.adjust", showCategory = 3,
        title = "KEGG Comparison", by = "geneRatio", size = "Count",
        facet = "direction", label_format = 60, font.size = 12) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))

# If not already a dataframe:
kegg_df <- as.data.frame(ekdegs)
spliceosome_rows <- kegg_df %>%
  filter(grepl("spliceosome", Description, ignore.case = TRUE))

print(spliceosome_rows)

# Remove any spliceosome-related entries EXCEPT for ko03041
kegg_df_filtered <- kegg_df %>%
  filter(!(grepl("spliceosome", Description, ignore.case = TRUE) & ID != "ko03041"))

ekdegs@compareClusterResult <- kegg_df_filtered

write.xlsx(kegg_df_filtered, "kegg_enrichment_filtered.xlsx")


#### Running each pop separately ####

deg_ruth <- read.csv("kegg_degs_ruth.csv")
deg_rw <- read.csv("kegg_degs_wils_rw.csv")
deg_fg <- read.csv("kegg_degs_wils_fg.csv")
deg_scot <- read.csv("kegg_degs_scot.csv")
deg_boon <- read.csv("kegg_degs_boon.csv")
deg_preb <- read.csv("kegg_degs_preb.csv")

ekdegs <- compareCluster(gene_id ~ direction + pop,
                         data = deg_preb,
                         fun = enricher,
                         TERM2GENE = pathway2gene_filtered2,
                         TERM2NAME = pathway2name_filtered2,
                         minGSSize = 1,
                         maxGSSize = 1000)

desired_order2 <- c("up", "down")

ekdegs@compareClusterResult$direction <- factor(ekdegs@compareClusterResult$direction, levels = desired_order2)


write.xlsx(as.data.frame(ekdegs), file = "kegg_enrichment_woHD_preb.xlsx")


dotplot(ekdegs, x = "pop", color = "p.adjust", showCategory = 20,
        title = "KEGG Comparison", by = "geneRatio", size = "Count",
        facet = "direction", label_format = 60, font.size = 12)

dotplot(ekdegs, x = "pop", color = "p.adjust", showCategory = 8,
        title = "KEGG Comparison", by = "geneRatio", size = "Count",
        facet = "direction", label_format = 60, font.size = 12) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))

# Remove any spliceosome-related entries EXCEPT for ko03041
# If not already a dataframe:
kegg_df1 <- as.data.frame(ekdegs)
spliceosome_rows <- kegg_df1 %>%
  filter(grepl("spliceosome", Description, ignore.case = TRUE))


kegg_df_filtered1 <- kegg_df1 %>%
  filter(!(grepl("spliceosome", Description, ignore.case = TRUE) & ID != "ko03041"))

ekdegs@compareClusterResult <- kegg_df_filtered1
