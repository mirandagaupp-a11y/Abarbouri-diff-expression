##### CLEANED FULL KEGG ENRICHMENT SCRIPT #####

# Load required libraries
library(jsonlite)
library(tibble)
library(stringr)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(openxlsx)
library(ggplot2)

# ----------------------------
# SECTION 1 — Build pathway2gene and pathway2name from KEGG JSON
# ----------------------------

update_kegg <- function(json = "ko00001.json") {
  pathway2name <- tibble(Pathway = character(), Name = character())
  ko2pathway <- tibble(Ko = character(), Pathway = character())
  kegg <- fromJSON(json)
  
  for (a in seq_along(kegg[["children"]][["children"]])) {
    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
        pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
        pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
        pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "")
        pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
        
        kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
        kos <- str_match(kos_info, "K[0-9]+")[,1]
        ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
      }
    }
  }
  
  save(pathway2name, ko2pathway, file = "kegg_info.RData")
  write.table(pathway2name, file = "kegg_pathway2name_gene.tsv", sep = "\t", quote = F, row.names = F)
  write.table(ko2pathway, file = "kegg_pathway2gene.tsv", sep = "\t", quote = F, row.names = F)
}

# Run once to build pathway files
update_kegg("ko00001.json")

# ----------------------------
# SECTION 2 — Build pathway2gene from gene2ko
# ----------------------------

# Load gene2ko
gene2ko <- read.csv("eggnog_kegg_one_per_line.csv") %>%
  dplyr::rename(GID = gene_id, Ko = KEGG_ko) %>%
  dplyr::filter(Ko != "-" & !is.na(Ko))

# Load pathway mappings
load("kegg_info.RData")  # loads ko2pathway, pathway2name

# Build pathway2gene
pathway2gene <- gene2ko %>%
  left_join(ko2pathway, by = "Ko") %>%
  dplyr::filter(!is.na(Pathway)) %>%
  dplyr::select(Pathway, GID)

# Save pathway2gene
write.table(pathway2gene, file = "kegg_pathway2gene.tsv", sep = "\t", quote = F, row.names = F)

# ----------------------------
# SECTION 3 — Enrichment function
# ----------------------------

run_kegg_enrichment <- function(gene_list, name) {
  enrich <- enricher(
    gene_list,
    TERM2GENE = pathway2gene,
    TERM2NAME = pathway2name,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 1
  )
  
  write.xlsx(as.data.frame(enrich), file = paste0(name, "_kegg_enrichment.xlsx"))
  
  p <- dotplot(enrich, x = "RichFactor", color = "p.adjust", showCategory = 20,
               title = paste(name, "KEGG enrichment"), size = "GeneRatio") +
    theme(axis.text.y = element_text(size = 10))
  
  print(p)
}

# ----------------------------
# SECTION 4 — Global DEGs (up/down)
# ----------------------------

# Load pathway2gene + pathway2name
pathway2gene <- read.csv("kegg_pathway2gene.csv", header = TRUE)
pathway2name <- read.table("kegg_pathway2name_gene", header = TRUE, sep = "\t")

# Load DEG lists
up_list <- read.csv("deg_list_up.csv")
down_list <- read.csv("deg_list_down.csv")

# Run enrichment for up/down DEGs
run_kegg_enrichment(up_list$gene_id, "DEG_up")
run_kegg_enrichment(down_list$gene_id, "DEG_down")

# ----------------------------
# SECTION 5 — Loop over population DEG lists
# ----------------------------

pop_lists <- list(
  preb = read.csv("kegg_degs_preb.csv"),
  boon = read.csv("kegg_degs_boon.csv"),
  scot = read.csv("kegg_degs_scot.csv"),
  fg   = read.csv("kegg_degs_wils_fg.csv"),
  rw   = read.csv("kegg_degs_wils_rw.csv"),
  ruth = read.csv("kegg_degs_ruth.csv")
)

for (pop in names(pop_lists)) {
  pop_list <- pop_lists[[pop]]
  run_kegg_enrichment(pop_list$gene_id, paste0("DEG_", pop))
}

# ----------------------------
# SECTION 6 — CompareCluster across pops
# ----------------------------

deg_combined <- read.csv("./deg_lists/kegg_degs.csv")

ekdegs <- compareCluster(gene_id ~ direction + pop,
                         data = deg_combined,
                         fun = enricher,
                         TERM2GENE = pathway2gene,
                         TERM2NAME = pathway2name,
                         minGSSize = 1,
                         maxGSSize = 1000)

desired_order <- c("Preble", "Boone", "Scott", "Wils_FG", 
                   "Wils_RW", "Rutherford")

ekdegs@compareClusterResult$pop <- factor(ekdegs@compareClusterResult$pop, levels = desired_order)

desired_order2 <- c("up", "down")

ekdegs@compareClusterResult$direction <- factor(ekdegs@compareClusterResult$direction, levels = desired_order2)


write.xlsx(as.data.frame(ekdegs), file = "kegg_comparison_enrichment.xlsx")


dotplot(ekdegs, x = "pop", color = "p.adjust", showCategory = 8,
        title = "KEGG Comparison", by = "geneRatio", size = "Count",
        facet = "direction", label_format = 60, font.size = 12)

dotplot(ekdegs, x = "pop", color = "p.adjust", showCategory = 8,
        title = "KEGG Comparison", by = "geneRatio", size = "Count",
        facet = "direction", label_format = 60, font.size = 12) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))
