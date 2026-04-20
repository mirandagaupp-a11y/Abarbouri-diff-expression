# Load dplyr & tidyr
library(dplyr)
library(tidyr)

# Read your full annotation file
annotations <- read.delim("/work/projects/5-39268/mgaupp42/transcriptome_assembly/trinotate/final_trinotate.out",
                          header = TRUE, sep = "\t")

# Clean and un-nest the KEGG pathways
kegg_data_clean <- annotations %>%
  select(gene_id, EggNM.KEGG_Pathway) %>%
  filter(EggNM.KEGG_Pathway != "." & EggNM.KEGG_Pathway != "") %>%
  separate_rows(EggNM.KEGG_Pathway, sep = ",") %>%
  mutate(EggNM.KEGG_Pathway = trimws(EggNM.KEGG_Pathway))

# Write the cleaned version to file (this is what you copy to desktop)
write.table(kegg_data_clean, "/work/projects/5-39268/mgaupp42/transcriptome_assembly/trinotate/kegg_annotation_clean.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
