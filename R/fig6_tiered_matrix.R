# ============================================================================
# fig6_tiered_matrix.R - Figure 6: Tiered Alignment Table
# Synthesis and convergence check with visual grouping
# Row-level labels, tier-block shading, minimal glyphs
# ============================================================================

source("R/00_setup_pub.R")

# --- Matrix Data (Grouped by Epistemic Role) ---
matrix_data <- tibble(
  Block = c(
    "Descriptive context",
    "Descriptive context",
    "Core inference",
    "Core inference",
    "Core inference",
    "Validation & robustness",
    "Validation & robustness"
  ),
  Row_Label = c(
    "Target count",
    "Proximity to DILI genes",
    "Network influence (RWI)",
    "Expression-weighted (EWI)",
    "Per-target efficiency (PTNI)",
    "Bootstrap sensitivity",
    "Chemical similarity"
  ),
  Hyperforin = c("Fewer (10)", "Farther", "Higher", "Higher", "Higher", "Robust", "None"),
  Quercetin = c("More (62)", "Closer", "Lower", "Lower", "Lower", "—", "None"),
  Hyp_Glyph = c("↓", "↓", "↑", "↑", "↑", "✔", "✔"),
  Que_Glyph = c("↑", "↑", "↓", "↓", "↓", "—", "✔")
)

# Assign row positions (bottom to top: 1-7)
matrix_data <- matrix_data %>%
  mutate(Row_Pos = 8 - row_number())

# --- Block shading coordinates ---
block_shades <- tibble(
  Block = c("Descriptive context", "Core inference", "Validation & robustness"),
  ymin = c(5.5, 2.5, 0.5),
  ymax = c(7.5, 5.5, 2.5),
  fill = c("#F0F0F0", "#FFFFFF", "#F8F8F8")
)

# --- Plot ---
p <- ggplot() +
  
  # Block shading
  geom_rect(
    data = block_shades,
    aes(xmin = 0.5, xmax = 3.5, ymin = ymin, ymax = ymax, fill = fill),
    color = NA
  ) +
  scale_fill_identity() +
  
  # Horizontal block separators
  geom_hline(yintercept = c(2.5, 5.5), color = "#888888", linewidth = 0.7) +
  
  # Cell borders
  geom_tile(
    data = matrix_data %>% 
      pivot_longer(c(Hyperforin, Quercetin), names_to = "Compound", values_to = "Value") %>%
      mutate(x = ifelse(Compound == "Hyperforin", 2, 3)),
    aes(x = x, y = Row_Pos),
    width = 1, height = 1,
    fill = NA, color = "#CCCCCC", linewidth = 0.3
  ) +
  
  # Row labels (column 1)
  geom_text(
    data = matrix_data,
    aes(x = 1, y = Row_Pos, label = Row_Label),
    hjust = 1, size = 3.2, family = "Arial", color = "#2E3440"
  ) +
  
  # Hyperforin values (column 2)
  geom_text(
    data = matrix_data,
    aes(x = 2, y = Row_Pos, label = Hyperforin),
    hjust = 0.5, size = 3.2, family = "Arial", color = "#2E3440"
  ) +
  
  # Hyperforin glyphs (top-right of cell)
  geom_text(
    data = matrix_data,
    aes(x = 2.42, y = Row_Pos + 0.38, label = Hyp_Glyph),
    hjust = 1, vjust = 1, size = 2.8, family = "Arial", fontface = "bold", color = "#555555"
  ) +
  
  # Quercetin values (column 3)
  geom_text(
    data = matrix_data,
    aes(x = 3, y = Row_Pos, label = Quercetin),
    hjust = 0.5, size = 3.2, family = "Arial", color = "#2E3440"
  ) +
  
  # Quercetin glyphs (top-right of cell)
  geom_text(
    data = matrix_data,
    aes(x = 3.42, y = Row_Pos + 0.38, label = Que_Glyph),
    hjust = 1, vjust = 1, size = 2.8, family = "Arial", fontface = "bold", color = "#555555"
  ) +
  
  # Column headers
  annotate("text", x = 2, y = 7.9, label = "Hyperforin", 
           size = 4, fontface = "bold", family = "Arial", color = "#2E3440") +
  annotate("text", x = 3, y = 7.9, label = "Quercetin", 
           size = 4, fontface = "bold", family = "Arial", color = "#2E3440") +
  
  # Block labels (left margin)
  annotate("text", x = -0.1, y = 6.5, label = "Descriptive\ncontext",
           hjust = 1, size = 2.6, fontface = "italic", family = "Arial", 
           color = "#666666", lineheight = 0.85) +
  annotate("text", x = -0.1, y = 4, label = "Core\ninference",
           hjust = 1, size = 2.6, fontface = "italic", family = "Arial", 
           color = "#666666", lineheight = 0.85) +
  annotate("text", x = -0.1, y = 1.5, label = "Validation &\nrobustness",
           hjust = 1, size = 2.6, fontface = "italic", family = "Arial", 
           color = "#666666", lineheight = 0.85) +
  
  # Scales
  scale_x_continuous(limits = c(-0.5, 3.6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, 8.3), expand = c(0, 0)) +
  
  labs(
    title = "Tiered evidence alignment",
    subtitle = "Directional agreement across independent analytical layers",
    caption = str_wrap(
      "Comparative synthesis across analytical tiers. Tiers represent distinct epistemic roles: descriptive context establishes baseline differences; core inference provides causal evidence; validation confirms robustness. Matrix is not used for inference or aggregation. Glyphs: ↑ = Higher/Stronger, ↓ = Lower/Weaker, ✔ = Pass.",
      width = 95
    )
  ) +
  
  theme_void(base_size = 11, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#4A4A4A", margin = margin(b = 12)),
    plot.caption = element_text(size = 8.5, hjust = 0, color = "#5A5A5A", 
                                lineheight = 1.3, margin = margin(t = 15)),
    plot.margin = margin(15, 20, 15, 60)
  ) +
  
  coord_fixed(ratio = 0.7, clip = "off")

# --- Display & Save ---
print(p)

ggsave(
  filename = here("figures", "main", "fig6_tiered_matrix.pdf"),
  plot = p,
  width = 180, height = 150, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = here("figures", "main", "fig6_tiered_matrix.tiff"),
  plot = p,
  width = 180, height = 150, units = "mm",
  dpi = 300, compression = "lzw"
)

message("✓ Figure 6 (Tiered Alignment Table) saved")
