# ============================================================================
# fig2_dumbbell.R - Figure 2: Proximity ≠ Influence (Dumbbell Chart)
# Shows rank reversal: Quercetin wins proximity, loses influence
# ============================================================================

source("R/00_setup_pub.R")
source("R/01_load_data.R")

# --- Data Preparation ---
# Need to reshape data for dumbbell plot
data_long <- tibble(
  Compound = rep(c("Hyperforin", "Quercetin"), each = 2),
  Metric = rep(c("Proximity", "Influence"), 2),
  Z_score = c(
    # Hyperforin
    sp_900 %>% filter(compound == "Hyperforin") %>% pull(z_score),
    rwr_900 %>% filter(compound == "Hyperforin") %>% pull(z_score),
    # Quercetin
    sp_900 %>% filter(compound == "Quercetin") %>% pull(z_score),
    rwr_900 %>% filter(compound == "Quercetin") %>% pull(z_score)
  )
) %>%
  mutate(
    Compound = factor(Compound, levels = c("Hyperforin", "Quercetin")), # Consistent order: Hyperforin first
    Metric = factor(Metric, levels = c("Proximity", "Influence"))
  )

# Wide format for connecting lines
data_wide <- data_long %>%
  pivot_wider(names_from = Metric, values_from = Z_score)

# --- Dumbbell Plot ---
p <- ggplot(data_wide, aes(y = Compound)) +
  # Connecting lines (the "dumbbells")
  geom_segment(
    aes(x = Proximity, xend = Influence, y = Compound, yend = Compound),
    color = "#3E3E3E", linewidth = 1.5, lineend = "round"
  ) +
  
  # Points for Proximity (left)
  geom_point(
    aes(x = Proximity), 
    size = 5, color = "#2E3440", alpha = 0.95
  ) +
  
  # Points for Influence (right)
  geom_point(
    aes(x = Influence), 
    size = 5, color = "#2E3440", alpha = 0.95
  ) +
  
  # Value labels for Proximity
  geom_text(
    aes(x = Proximity, label = sprintf("Z = %.1f", Proximity)),
    hjust = 1.3, size = 3.3, fontface = "bold", family = "Arial", color = "#2E3440"
  ) +
  
  # Value labels for Influence
  geom_text(
    aes(x = Influence, label = sprintf("Z = +%.1f", Influence)),
    hjust = -0.3, size = 3.3, fontface = "bold", family = "Arial", color = "#2E3440"
  ) +
  
  # Reference line at zero
  geom_vline(xintercept = 0, linetype = "dashed", color = "#AAAAAA", linewidth = 0.5) +
  
  # Column headers (with explicit directional cues)
  annotate("text", x = -5.5, y = 2.6, label = "Proximity\n← more negative = closer", 
           hjust = 0.5, vjust = 0, size = 3.5, fontface = "bold", 
           family = "Arial", color = "#4A4A4A", lineheight = 0.9) +
  annotate("text", x = 9, y = 2.6, label = "Influence\nmore positive = stronger →", 
           hjust = 0.5, vjust = 0, size = 3.5, fontface = "bold", 
           family = "Arial", color = "#4A4A4A", lineheight = 0.9) +
  
  # Scales
  scale_x_continuous(
    limits = c(-7, 12),
    breaks = seq(-6, 12, 3),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  labs(
    title = "Proximity does not predict influence",
    subtitle = "Proximity ranking is threshold-dependent; influence ranking is stable",
    x = "Z-score",
    y = NULL,
    caption = str_wrap(
      "[CORE INFERENCE] The rank reversal demonstrates that shortest-path proximity does not predict functional influence. Lines connect each compound's proximity Z-score with its influence Z-score (random walk with restart, RWI). Both metrics derived from degree-matched permutation null models (n = 1,000). Data: STRING v12.0 (≥900).",
      width = 100
    )
  ) +
  
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    # Axes
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 10)),
    
    # Titles
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, 
                              margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#4A4A4A", 
                                 margin = margin(b = 15)),
    plot.caption = element_text(size = 8.5, hjust = 0, color = "#5A5A5A", 
                                lineheight = 1.3, margin = margin(t = 12)),
    
    # Grid
    panel.grid.major.x = element_line(color = "#F0F0F0", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    
    # Margins
    plot.margin = margin(15, 15, 15, 15)
  ) +
  
  # Expand y-axis to fit headers
  coord_cartesian(ylim = c(0.5, 2.7), clip = "off")

# --- Display & Save ---
print(p)

ggsave(
  filename = here("figures", "main", "fig2_dumbbell.pdf"),
  plot = p,
  width = 180, height = 110, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = here("figures", "main", "fig2_dumbbell.tiff"),
  plot = p,
  width = 180, height = 110, units = "mm",
  dpi = 300, compression = "lzw"
)

message("✓ Figure 2 (Dumbbell) saved - Rank reversal clearly visualized")
