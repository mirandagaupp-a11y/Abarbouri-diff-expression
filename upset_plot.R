install.packages("UpSetR")
library(UpSetR)

degs <- read.csv("degs_by_pop.csv", stringsAsFactors = FALSE)

gene_sets <- list(
  Rutherford = na.omit(degs$Rutherford),
  Wilson2 = na.omit(degs$Wilson2),
  Wilson1 = na.omit(degs$Wilson1),
  Boone = na.omit(degs$Boone),
  Scott = na.omit(degs$Scott),
  Preble = na.omit(degs$Preble)
)

shared_all_6 <- Reduce(intersect, gene_sets)
shared_all_6_df <- data.frame(Gene = shared_all_6)

write.csv(shared_all_6_df,
          "genes_shared_all_6_populations.csv",
          row.names = FALSE)

jpeg("upset_plot.jpg",
     width = 13, height = 7,
     units = "in", res = 600)

upset(
  fromList(gene_sets),
  sets = c("Rutherford","Wilson2", "Wilson1", "Boone", "Scott", "Preble"),
  nsets = 6,
  nintersects = 35,
  keep.order = TRUE,
  order.by = "freq",
  decreasing = TRUE,
  main.bar.color = "dodgerblue",
  sets.bar.color = "orange",
  matrix.color = "gray20",
  point.size = 2.5,
  line.size = 0.7,
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size",
  text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5)
)

dev.off()

upset(
  fromList(gene_sets),
  sets = c("Rutherford","Wilson2", "Wilson1", "Boone", "Scott", "Preble"),
  nsets = 6,
  nintersects = 70,
  keep.order = TRUE,
  order.by = "freq",
  decreasing = TRUE,
  main.bar.color = "dodgerblue",
  sets.bar.color = "orange",
  matrix.color = "gray20",
  point.size = 2.5,
  line.size = 0.7,
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size",
  text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5)
)
