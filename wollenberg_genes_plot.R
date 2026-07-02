library(tidyverse)
library(scales)
library(ggplot2)

dat <- read_csv("wollenberg_genes_final.csv", show_col_types = FALSE)

pop_order <- c("Preble","Boone","Scott","Wilson1","Wilson2","Rutherford")

prep_plot_df <- function(raw, direction = c("up","down")) {
  
  direction <- match.arg(direction)
  
  base <- raw %>%
    mutate(
      pop = factor(pop, levels = pop_order),
      size_val = abs(log2FoldChange),
      gene_label = gene_info,
      sig_up = !is.na(padj) & padj <= 0.05 & log2FoldChange > 0,
      sig_down = !is.na(padj) & padj <= 0.05 & log2FoldChange < 0
    )
  
  sig_col <- if (direction == "up") "sig_up" else "sig_down"
  
  genes_any <- base %>%
    group_by(gene_id) %>%
    summarise(any_sig = any(.data[[sig_col]], na.rm = TRUE),
              .groups = "drop") %>%
    filter(any_sig)
  
  base %>%
    inner_join(genes_any, by = "gene_id") %>%
    mutate(
      sig = .data[[sig_col]],
      padj_plot = if_else(sig, padj, NA_real_)
    ) %>%
    filter(!is.na(size_val))
}

df_up <- prep_plot_df(dat, "up")
df_down <- prep_plot_df(dat, "down")

size_limits <- range(c(df_up$size_val, df_down$size_val), na.rm = TRUE)

make_plot <- function(d, title_txt) {
  
  ggplot(d, aes(pop, gene_label, size = size_val)) +
    
    geom_point(
      data = filter(d, sig),
      aes(color = pmin(padj_plot, 0.05)),
      alpha = 0.9
    ) +
    
    geom_point(
      data = filter(d, !sig),
      color = "grey75",
      alpha = 0.4
    ) +
    
  scale_size_continuous(name = "|log2FC|", range = c(2,10), breaks = c(0.5,1,2), 
                        limits = size_limits) +
    
  scale_color_gradientn(
  name = "padj (≤ 0.05)",
  colours = c("#B2182B", "#EF8A62", "#FDAE61", "#FFD000"),
  limits = c(0, 0.05),
  oob = squish,
  labels = label_number(accuracy = 0.001)
) +
    
    labs(
      title = title_txt,
      x = "pop",
      y = NULL
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = -45, hjust = 0, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA)
    )
}

p_up <- make_plot(df_up, "Candidate Stress Genes: Up")
p_down <- make_plot(df_down, "Candidate Stress Genes: Down")

p_up
p_down

ggsave("wb_DE_gene_dotplot_UP.jpg", p_up,
       width = 6, height = 8, units = "in",
       dpi = 600, device = "jpeg", bg = "white")

ggsave("wb_DE_gene_dotplot_DOWN.jpg", p_down,
       width = 6, height = 8, units = "in",
       dpi = 600, device = "jpeg", bg = "white")

install.packages("patchwork")
library(patchwork)


combined_plot <- (p_up + p_down) + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right", legend.box = "vertical", plot.margin = margin(5, 5, 5, 5))

ggsave("wb_DE_gene_dotplots_combined.jpg", combined_plot,
       width = 11, height = 8, units = "in", dpi = 600, device = "jpeg", bg = "white")
