# ============================================================================
# fig4_ptni_phase.R - Figure 4: Efficiency-Coverage Phase Plot (Faceted)
# PTNI as geometric slope, not just a ratio
# x = Target count, y = Total influence, iso-PTNI diagonal lines
# Faceted: RWI | EWI with identical axes
# ============================================================================

source("R/00_setup_pub.R")
source("R/01_load_data.R")

# --- Data Preparation (Both RWI and EWI) ---
phase_data <- bind_rows(
  tibble(
    Compound = c("Hyperforin", "Quercetin"),
    Targets = c(10, 62),
    Influence = c(
      rwr_900 %>% filter(compound == "Hyperforin") %>% pull(observed_influence),
      rwr_900 %>% filter(compound == "Quercetin") %>% pull(observed_influence)
    ),
    Method = "RWI"
  ),
  tibble(
    Compound = c("Hyperforin", "Quercetin"),
    Targets = c(10, 62),
    Influence = c(
      ewr_900 %>% filter(compound == "Hyperforin") %>% pull(observed_influence),
      ewr_900 %>% filter(compound == "Quercetin") %>% pull(observed_influence)
    ),
    Method = "EWI"
  )
) %>%
  mutate(
    PTNI = Influence / Targets,
    Method = factor(Method, levels = c("RWI", "EWI"))
  )

# --- Calculate Iso-PTNI Lines ---
# PTNI = Influence / Targets
# Therefore: Influence = PTNI × Targets
# For diagonal reference lines

ptni_levels <- c(0.001, 0.005, 0.01, 0.02, 0.05)  # Efficiency contours
x_range <- seq(0, 70, length.out = 100)

# Duplicate iso-lines for both facets
iso_lines <- expand_grid(
  Method = factor(c("RWI", "EWI"), levels = c("RWI", "EWI")),
  PTNI = ptni_levels
) %>%
  rowwise() %>%
  mutate(data = list(tibble(x = x_range, y = PTNI * x_range))) %>%
  unnest(data)

# Iso-line labels (right edge only)
iso_labels <- iso_lines %>%
  filter(x == max(x)) %>%
  distinct(Method, PTNI, x, y)

# --- Faceted Phase Plot ---
p <- ggplot() +
  # Iso-PTNI reference lines (lightly drawn)
  geom_line(
    data = iso_lines,
    aes(x = x, y = y, group = factor(PTNI)),
    color = "#D0D0D0", linewidth = 0.4, linetype = "dashed"
  ) +
  
  # Iso-PTNI labels (at right edge)
  geom_text(
    data = iso_labels,
    aes(x = x, y = y, label = sprintf("%.3f", PTNI)),
    hjust = -0.05, size = 2.5, color = "#808080", family = "Arial", fontface = "italic"
  ) +
  
  # Data points (compounds)
  geom_point(
    data = phase_data,
    aes(x = Targets, y = Influence),
    size = 5, color = "#2E3440", alpha = 0.95
  ) +
  
  # Compound labels
  geom_text_repel(
    data = phase_data,
    aes(x = Targets, y = Influence, 
        label = paste0(Compound, "\n(PTNI = ", sprintf("%.4f", PTNI), ")")),
    size = 3.5, fontface = "bold", family = "Arial",
    box.padding = 1.2, point.padding = 0.5,
    force = 2, max.overlaps = Inf,
    segment.color = "#666666", segment.size = 0.4,
    lineheight = 0.9
  ) +
  
  # Facet by method
  facet_wrap(~ Method, ncol = 2) +
  
  # Scales (identical for both panels)
  scale_x_continuous(
    limits = c(0, 70),
    breaks = seq(0, 70, 20),
    expand = expansion(mult = c(0.02, 0.2)) # More space for right labels
  ) +
  
  scale_y_continuous(
    limits = c(0, 0.12),
    breaks = seq(0, 0.12, 0.04),
    expand = expansion(mult = c(0.02, 0.1)) # More space for top labels
  ) +
  
  labs(
    title = "Per-target efficiency: coverage versus influence",
    subtitle = "Hyperforin occupies a higher efficiency contour despite fewer targets",
    x = "Target count",
    y = "Total influence",
    caption = str_wrap(
      "[EFFICIENCY NORMALIZATION] PTNI = influence / target count (derived metric, not formally tested). This reframes polypharmacology as an efficiency problem: Hyperforin achieves greater influence per target. PTNI ratio differences are consistent with the Influence Z-score hierarchy (Figs 1–3). Diagonal lines show iso-efficiency contours. Data: STRING v12.0 (≥900), n = 1,000 permutations.",
      width = 120
    )
  ) +
  
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    # Axes
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 11),
    
    # Titles
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#4A4A4A", margin = margin(b = 15)),
    plot.caption = element_text(size = 8.5, hjust = 0, color = "#5A5A5A", 
                                lineheight = 1.3, margin = margin(t = 12)),
    
    # Facet strips
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold", color = "#2E3440"),
    
    # Grid (minimal)
    panel.grid.major = element_line(color = "#F5F5F5", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.5, "cm"),
    
    # Margins
    plot.margin = margin(15, 20, 15, 15)
  )

# --- Display & Save ---
print(p)

ggsave(
  filename = here("figures", "main", "fig4_ptni_phase.pdf"),
  plot = p,
  width = 280, height = 140, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = here("figures", "main", "fig4_ptni_phase.tiff"),
  plot = p,
  width = 280, height = 140, units = "mm",
  dpi = 300, compression = "lzw"
)

message("✓ Figure 4 (PTNI Phase Plot - Faceted RWI|EWI) saved")

