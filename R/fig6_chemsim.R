# ============================================================================
# fig6_chemsim.R - Figure 6: Chemical Similarity Control
# Defensive figure: Close the "maybe it's just chemistry" escape hatch
# Simple threshold plot, intentionally boring
# ============================================================================

source("R/00_setup_pub.R")

# --- Load Data ---
chemsim <- read_csv(here("results", "tables", "chemical_similarity_summary.csv"))

# --- Prepare Plot Data ---
plot_data <- chemsim %>%
  select(compound, max_sim_DILI_positive, max_sim_DILI_negative) %>%
  pivot_longer(
    cols = c(max_sim_DILI_positive, max_sim_DILI_negative),
    names_to = "Reference",
    values_to = "Max_Tanimoto"
  ) %>%
  mutate(
    Reference = case_when(
      Reference == "max_sim_DILI_positive" ~ "DILI-positive drugs",
      Reference == "max_sim_DILI_negative" ~ "DILI-negative drugs"
    ),
    Reference = factor(Reference, levels = c("DILI-positive drugs", "DILI-negative drugs")),
    compound = factor(compound, levels = c("Hyperforin", "Quercetin"))
  )

# --- Threshold ---
similarity_threshold <- 0.4

# --- Plot ---
p <- ggplot(plot_data, aes(x = compound, y = Max_Tanimoto, fill = Reference)) +
  
  # Threshold line
  geom_hline(
    yintercept = similarity_threshold,
    linetype = "dashed", color = "#888888", linewidth = 0.7
  ) +
  
  # Threshold label (using geom_label for discrete x-axis)
  geom_label(
    data = data.frame(compound = "Hyperforin", y = similarity_threshold),
    aes(x = compound, y = y),
    label = "Threshold (0.4)",
    inherit.aes = FALSE,
    hjust = 0, vjust = -0.3, nudge_x = -0.3,
    size = 2.8, fontface = "italic", family = "Arial", color = "#666666",
    fill = "white", label.size = 0
  ) +
  
  # Bars
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.6, color = "#2E3440", linewidth = 0.4
  ) +
  
  # Value labels
  geom_text(
    aes(label = sprintf("%.2f", Max_Tanimoto)),
    position = position_dodge(width = 0.7),
    vjust = -0.5, size = 3.2, fontface = "bold", family = "Arial", color = "#2E3440"
  ) +
  
  # Scales
  scale_y_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, 0.1),
    expand = expansion(mult = c(0, 0.1))
  ) +
  
  scale_fill_manual(
    values = c("DILI-positive drugs" = "#666666", "DILI-negative drugs" = "#D8D8D8"),
    name = "Reference set"
  ) +
  
  labs(
    title = "Chemical similarity control",
    subtitle = "Neither compound resembles known hepatotoxins",
    x = NULL,
    y = "Maximum Tanimoto similarity",
    caption = str_wrap(
      "[ORTHOGONAL EXCLUSION] Maximum Tanimoto similarity to DILIrank drugs. Threshold 0.4 is conservative for structural analog detection. Both compounds fall well below. Reference: DILIrank (542 DILI+, 365 DILI−).",
      width = 85
    )
  ) +
  
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    # Axes
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.line.x = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_blank(),
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_text(face = "bold", size = 11),
    axis.title.y = element_text(face = "bold", size = 11, margin = margin(r = 8)),
    
    # Titles
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#4A4A4A", margin = margin(b = 15)),
    plot.caption = element_text(size = 8.5, hjust = 0, color = "#5A5A5A", 
                                lineheight = 1.3, margin = margin(t = 12)),
    
    # Legend
    legend.position = "top",
    legend.justification = "center",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.5, "cm"),
    legend.margin = margin(b = 10),
    
    # Grid
    panel.grid.major.y = element_line(color = "#F0F0F0", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    
    # Margins
    plot.margin = margin(15, 20, 15, 25)
  ) +
  
  coord_cartesian(clip = "off")

# --- Display & Save ---
print(p)

ggsave(
  filename = here("figures", "main", "fig6_chemsim.pdf"),
  plot = p,
  width = 160, height = 130, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = here("figures", "main", "fig6_chemsim.tiff"),
  plot = p,
  width = 160, height = 130, units = "mm",
  dpi = 300, compression = "lzw"
)

message("✓ Figure 6 (Chemical Similarity Control) saved - Defensive closure")

