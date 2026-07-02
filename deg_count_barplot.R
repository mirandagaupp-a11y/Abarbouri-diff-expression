library(ggplot2)

# Create data frame
counts <- data.frame(
  County = c("Preble", "Boone", "Scott", "Wilson1", "Wilson2", "Rutherford"),
  Up = c(404, 222, 204, 690, 703, 682),
  Down = c(193, 684, 519, 583, 1564, 1383)
)

# Make downregulated values negative
counts$County <- factor(counts$County, levels = c("Preble", "Boone", "Scott", "Wilson1", "Wilson2", "Rutherford"))
counts$Down_plot <- -counts$Down

q <- ggplot(counts) +
  geom_col(aes(County, Up, fill = "Up"), color = "black", linewidth = 0.2, width = 0.9) +
  geom_col(aes(County, Down_plot, fill = "Down"), color = "black", linewidth = 0.2, width = 0.9) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "black") +
  scale_fill_manual(values = c("Down" = "blue", "Up" = "grey75"),
                    name = "Regulation") +
  scale_y_continuous(labels = function(x) abs(x),
                     limits = c(-1700, 1000),
                     breaks = seq(-1500, 1000, 500)) +
  labs(x = NULL, y = "Number of Differentially Expressed Genes") +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "right",
    axis.text = element_text(color = "black")
  )

q 

# Save the figure
ggsave("degs_count.jpg", q,
       width = 9, height = 6, units = "in", dpi = 600, device = "jpeg", bg = "white")
