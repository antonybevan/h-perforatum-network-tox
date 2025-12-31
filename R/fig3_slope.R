# ============================================================================
# fig3_slope.R - Figure 3: Cleveland-style Slopegraph
# RWI → EWI: Shows ranking stability under biological constraint
# ============================================================================

source("R/00_setup_pub.R")
source("R/01_load_data.R")

# --- Data Preparation ---
slope_data <- tibble(
  Compound = c("Hyperforin", "Quercetin"),
  RWI = c(
    rwr_900 %>% filter(compound == "Hyperforin") %>% pull(z_score),
    rwr_900 %>% filter(compound == "Quercetin") %>% pull(z_score)
  ),
  EWI = c(
    ewr_900 %>% filter(compound == "Hyperforin") %>% pull(z_score),
    ewr_900 %>% filter(compound == "Quercetin") %>% pull(z_score)
  )
) %>%
  mutate(
    Compound = factor(Compound, levels = c("Hyperforin", "Quercetin"))
  )

# --- Cleveland Slopegraph ---
p <- ggplot(slope_data) +
  # Connecting lines (slopes)
  geom_segment(
    aes(x = 1, xend = 2, y = RWI, yend = EWI),
    color = "#3E3E3E", linewidth = 1.2
  ) +
  
  # Left points (RWI)
  geom_point(aes(x = 1, y = RWI), size = 4.5, color = "#2E3440") +
  
  # Right points (EWI)
  geom_point(aes(x = 2, y = EWI), size = 4.5, color = "#2E3440") +
  
  # Left labels (compound names + values)
  geom_text(
    aes(x = 1, y = RWI, label = paste0(Compound, "\n(Z = +", round(RWI, 1), ")")),
    hjust = 1.15, size = 3.5, fontface = "bold", family = "Arial", color = "#2E3440",
    lineheight = 0.9
  ) +
  
  # Right labels (values only)
  geom_text(
    aes(x = 2, y = EWI, label = paste0("Z = +", round(EWI, 1))),
    hjust = -0.15, size = 3.5, fontface = "bold", family = "Arial", color = "#2E3440"
  ) +
  
  # Column headers (full form with abbreviations, small and aligned)
  annotate("text", x = 1, y = 13.5, label = "Random walk\nwith restart (RWR)", 
           hjust = 0.5, size = 3.5, fontface = "bold", 
           family = "Arial", color = "#2E3440", lineheight = 0.9) +
  annotate("text", x = 2, y = 13.5, label = "Expression-weighted\ninfluence (EWI)", 
           hjust = 0.5, size = 3.5, fontface = "bold", 
           family = "Arial", color = "#2E3440", lineheight = 0.9) +
  
  # Scales
  scale_x_continuous(
    limits = c(0.2, 2.8),  # More horizontal space
    breaks = NULL
  ) +
  scale_y_continuous(
    limits = c(3, 15),  # Extended upper limit for headers
    breaks = seq(4, 14, 2),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  labs(
    title = "Ranking stability under biological constraint",
    subtitle = "Expression weighting preserves the influence hierarchy",
    x = NULL,
    y = "Influence Z-score",
    caption = str_wrap(
      "[BIOLOGICAL VALIDATION] Expression weighting tests whether influence rankings are robust to tissue-specific constraints. Parallel slopes indicate ranking stability. Both methods use independent permutation null models (n = 1,000 each). GTEx liver expression (TPM ≥1) weights edge transitions in EWI. Data: STRING v12.0 (≥900).",
      width = 100
    )
  ) +
  
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    # Axes
    axis.line.x = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.y = element_text(face = "bold", size = 11),
    
    # Titles
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, 
                              margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#4A4A4A", 
                                  margin = margin(b = 20)),
    plot.caption = element_text(size = 8.5, hjust = 0, color = "#5A5A5A", 
                                lineheight = 1.3, margin = margin(t = 15)),
    
    # Grid
    panel.grid.major.y = element_line(color = "#F0F0F0", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    
    # Margins (extra room for labels)
    plot.margin = margin(15, 20, 15, 20)
  )

# --- Display & Save ---
print(p)

ggsave(
  filename = here("figures", "main", "fig3_slope.pdf"),
  plot = p,
  width = 160, height = 140, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = here("figures", "main", "fig3_slope.tiff"),
  plot = p,
  width = 160, height = 140, units = "mm",
  dpi = 300, compression = "lzw"
)

message("✓ Figure 3 (Cleveland Slopegraph) saved")
