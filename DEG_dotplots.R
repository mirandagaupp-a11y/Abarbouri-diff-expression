# ================================
# Dotplot of DE genes by population
# size = |log2FC| ; color = padj ; padj <= 0.05 only
# ================================


library(dplyr)
library(readr)
library(ggplot2)
library(forcats)
library(scales)

# -------- Inputs --------
in_csv  <- "wollenberg_genes_final.csv"   # change if needed
#in_csv  <- "hsp_DEs_for_dotplot.csv" 
out_png <- "wb_DE_gene_dotplot_by_pop.png"
pop_order <- c("Preble","Boone","Scott","Wilson1","Wilson2","Rutherford")

# -------- Read + normalize expected columns --------
# Expected columns (case-insensitive): gene_id, gene_info, pop, log2FoldChange, padj
df_raw <- read_csv(in_csv, show_col_types = FALSE)

# Try to map common variants to standard names
std_map <- c(
  "gene_id"        = "gene_id|geneid|gene",
  "gene_info"      = "gene_info|gene.info|description|desc|annotation",
  "pop"            = "pop|population",
  "log2FoldChange" = "log2foldchange|log2fc|log2_change|log2_fc",
  "padj"           = "padj|p_adj|p.adjust|p.adjusted|fdr|qvalue|q_value"
)

for (std in names(std_map)) {
  if (!std %in% names(df_raw)) {
    alt <- grep(std_map[[std]], names(df_raw), ignore.case = TRUE, value = TRUE)
    if (length(alt) == 1) names(df_raw)[match(alt, names(df_raw))] <- std
  }
}

needed <- c("gene_id","gene_info","pop","log2FoldChange","padj")
missing <- setdiff(needed, names(df_raw))
if (length(missing) > 0) {
  stop(paste0("Missing required column(s): ", paste(missing, collapse = ", "),
              ". Make sure your CSV has these fields (case-insensitive)."))
}

# -------- Filter + derive plotting variables --------
df <- df_raw %>%
  mutate(
    pop  = factor(pop, levels = pop_order),
    direction = ifelse(log2FoldChange >= 0, "up", "down"),
    direction = factor(direction, levels = c("up","down")),
    size_val  = abs(log2FoldChange),
    gene_label = ifelse(is.na(gene_info) | gene_info == "",
                        gene_id, gene_info)
  ) %>%
  filter(!is.na(padj), padj <= 0.05, !is.na(size_val), !is.na(pop))

# Optional: cap extreme sizes at 99th percentile so one big point doesn't dwarf the rest
cap <- quantile(df$size_val, 0.99, na.rm = TRUE)
df <- df %>% mutate(size_val = pmin(size_val, cap))

# If you only want the top N per facet, uncomment the next lines and set N:
# N <- 30
# df <- df %>%
#   group_by(direction) %>%
#   slice_max(order_by = size_val, n = N, with_ties = FALSE) %>%
#   ungroup()



padj_limits <- c(0, 0.05)
padj_cols   <- c("#D73027", "#FC8D59", "#91BFDB", "#4575B4")  # red -> blue, 0.05 darkest

p <- ggplot(df, aes(x = pop, y = gene_label, size = size_val, color = padj)) +
  geom_point(alpha = 0.9) +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y",
             switch = "y",                                   # put strips on the left
             labeller = labeller(direction = c(up = "UP", down = "DOWN"))) +
  scale_size_continuous(name = "|log2FC|",
                        range = c(2, 10),
                        breaks = c(0.5, 1, 2, 3)) +           # continuous sizing; labeled ticks
  scale_color_gradientn(name = "padj",
                        colours = padj_cols,
                        limits  = padj_limits,
                        oob     = squish,
                        labels  = label_number(accuracy = 0.001)) +
  labs(x = "pop", y = NULL, title = "WB DE Genes by Population (padj \u2264 0.05)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x        = element_text(angle = -45, hjust = 0),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA),   # box around each panel
    strip.placement    = "inside",                                   # <— place strips after labels
    strip.background.y = element_rect(fill = "grey85", color = "grey40"),
    strip.text.y.left  = element_text(angle = 90, face = "bold", vjust = 0.5)
  ) +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2))

p



library(dplyr)
library(readr)
library(ggplot2)
library(scales)

# -------- Inputs --------
in_csv   <- "wollenberg_DEs_for_dotplot2.csv"
out_up   <- "wb_DE_gene_dotplot_by_pop_UP.png"
out_down <- "wb_DE_gene_dotplot_by_pop_DOWN.png"
pop_order <- c("Preble","Boone","Scott","Wilson1","Wilson2","Rutherford")

# -------- Read + normalize expected columns --------
df_raw <- read_csv(in_csv, show_col_types = FALSE)

std_map <- c(
  "gene_id"        = "gene_id|geneid|gene",
  "gene_info"      = "gene_info|gene.info|description|desc|annotation",
  "pop"            = "pop|population",
  "log2FoldChange" = "log2foldchange|log2fc|log2_change|log2_fc",
  "padj"           = "padj|p_adj|p.adjust|p.adjusted|fdr|qvalue|q_value"
)
for (std in names(std_map)) {
  if (!std %in% names(df_raw)) {
    alt <- grep(std_map[[std]], names(df_raw), ignore.case = TRUE, value = TRUE)
    if (length(alt) == 1) names(df_raw)[match(alt, names(df_raw))] <- std
  }
}
needed <- c("gene_id","gene_info","pop","log2FoldChange","padj")
missing <- setdiff(needed, names(df_raw))
if (length(missing) > 0) {
  stop(paste0("Missing required column(s): ", paste(missing, collapse = ", "),
              ". Make sure your CSV has these fields (case-insensitive)."))
}

# -------- Helpers: make plot data per direction --------
prep_plot_df <- function(raw, direction = c("up","down")) {
  direction <- match.arg(direction)
  
  # Base table with flags for each direction
  base <- raw %>%
    dplyr::mutate(
      pop        = factor(pop, levels = pop_order),
      size_val   = abs(log2FoldChange),
      gene_label = ifelse(is.na(gene_info) | gene_info == "", gene_id, gene_info),
      sig_up     = !is.na(padj) & padj <= 0.05 & log2FoldChange >  0,
      sig_down   = !is.na(padj) & padj <= 0.05 & log2FoldChange <  0
    )
  
  if (direction == "up") {
    genes_any <- base %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(any_sig = any(sig_up, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(any_sig)
    
    out <- base %>%
      dplyr::inner_join(genes_any, by = "gene_id") %>%
      dplyr::mutate(
        sig       = sig_up,
        padj_plot = ifelse(sig, padj, NA_real_),
        dir_plot  = "UP"
      )
  } else {
    genes_any <- base %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(any_sig = any(sig_down, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(any_sig)
    
    out <- base %>%
      dplyr::inner_join(genes_any, by = "gene_id") %>%
      dplyr::mutate(
        sig       = sig_down,
        padj_plot = ifelse(sig, padj, NA_real_),
        dir_plot  = "DOWN"
      )
  }
  
  out %>%
    dplyr::filter(!is.na(size_val), !is.na(pop))
}


df_up   <- prep_plot_df(df_raw, "up")
df_down <- prep_plot_df(df_raw, "down")

# Cap extreme sizes (apply same cap to both so scales are comparable)
cap <- quantile(c(df_up$size_val, df_down$size_val), 0.99, na.rm = TRUE)
df_up   <- df_up   %>% mutate(size_val = pmin(size_val, cap))
df_down <- df_down %>% mutate(size_val = pmin(size_val, cap))

# -------- Plot factory (with darker grey for non-sig) --------
make_plot <- function(d, title_txt) {
  grey_col <- "#666666"
  padj_limits <- c(0, 0.05)
  padj_cols   <- c("#D73027", "#FC8D59", "#91BFDB", "#4575B4")
  
  ggplot(d, aes(x = pop, y = gene_label, size = size_val)) +
    # significant in this direction (colored)
    geom_point(data = d %>% filter(sig),
               aes(color = pmin(padj_plot, 0.05)),
               alpha = 0.9) +
    # not significant in this direction (grey)
    geom_point(data = d %>% filter(!sig),
               color = grey_col, alpha = 0.35) +
    scale_size_continuous(name = "|log2FC|", range = c(2, 10), breaks = c(0.5, 1, 2, 3)) +
    scale_color_gradientn(
      name   = "padj (≤ 0.05)",
      colours = padj_cols,
      limits  = padj_limits,
      oob     = squish,
      labels  = label_number(accuracy = 0.001)
    ) +
    labs(x = "pop", y = NULL, title = title_txt) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x        = element_text(angle = -45, hjust = 0),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "black", fill = NA)
    ) +
    guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2))
}

# -------- Make & save --------
p_up   <- make_plot(df_up,   "WB DE Genes by Population — UP (colored where padj ≤ 0.05 and LFC > 0)")
p_down <- make_plot(df_down, "WB DE Genes by Population — DOWN (colored where padj ≤ 0.05 and LFC < 0)")

p_up
p_down

