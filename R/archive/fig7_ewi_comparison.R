# ============================================================================
# fig7_ewi_comparison.R - Figure 7: RWR vs EWI Comparison (Paired Bar Chart)
# Shows Hyperforin's persistent advantage despite expression-weighted narrowing
# ============================================================================

source("R/00_setup_pub.R")
source("R/01_load_data.R")

# --- Data Preparation ---
comparison_data <- tibble(
  Compound = rep(c("Hyperforin", "Quercetin"), 2),
  Method = rep(c("Random walk with restart\n(RWR)", "Expression-weighted influence\n(EWI)"), each = 2),
  Z_score = c(
    # RWR
    rwr_900 %>% filter(compound == "Hyperforin") %>% pull(z_score),
    rwr_900 %>% filter(compound == "Quercetin") %>% pull(z_score),
    # EWI
    ewr_900 %>% filter(compound == "Hyperforin") %>% pull(z_score),
    ewr_900 %>% filter(compound == "Quercetin") %>% pull(z_score)
  )
) %>%
  mutate(
    Compound = factor(Compound, levels = c("Hyperforin", "Quercetin")),
    Method = factor(Method, levels = c("Random walk with restart\n(RWR)", "Expression-weighted influence\n(EWI)"))
  )

# --- Calculate Gaps ---
rwr_gap <- comparison_data %>%
  filter(Method == "Random walk with restart\n(RWR)") %>%
  summarise(gap = diff(rev(Z_score))) %>%
  pull(gap)

ewi_gap <- comparison_data %>%
  filter(Method == "Expression-weighted influence\n(EWI)") %>%
  summarise(gap = diff(rev(Z_score))) %>%
  pull(gap)

# --- Paired Bar Chart ---
p <- ggplot(comparison_data, aes(x = Method, y = Z_score, fill = Compound)) +
  
  # Bars
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.6, color = "#2E3440", linewidth = 0.4
  ) +
  
  # Value labels on bars
  geom_text(
    aes(label = sprintf("Z = +%.1f", Z_score)),
    position = position_dodge(width = 0.7),
    vjust = -0.5, size = 3.5, fontface = "bold", family = "Arial", color = "#2E3440"
  ) +
  
  # Gap annotations
  annotate(
    "segment",
    x = 0.65, xend = 1.35,
    y = max(comparison_data$Z_score[comparison_data$Method == "Random walk with restart\n(RWR)"]) + 1.5,
    yend = max(comparison_data$Z_score[comparison_data$Method == "Random walk with restart\n(RWR)"]) + 1.5,
    color = "#666666", linewidth = 0.5
  ) +
  annotate(
    "text",
    x = 1, y = max(comparison_data$Z_score[comparison_data$Method == "Random walk with restart\n(RWR)"]) + 2.2,
    label = sprintf("Gap: +%.1f", rwr_gap),
    size = 3.5, fontface = "bold", family = "Arial", color = "#2E3440"
  ) +
  
  annotate(
    "segment",
    x = 1.65, xend = 2.35,
    y = max(comparison_data$Z_score[comparison_data$Method == "Expression-weighted influence\n(EWI)"]) + 1.5,
    yend = max(comparison_data$Z_score[comparison_data$Method == "Expression-weighted influence\n(EWI)"]) + 1.5,
    color = "#666666", linewidth = 0.5
  ) +
  annotate(
    "text",
    x = 2, y = max(comparison_data$Z_score[comparison_data$Method == "Expression-weighted influence\n(EWI)"]) + 2.2,
    label = sprintf("Gap: +%.1f", ewi_gap),
    size = 3.5, fontface = "bold", family = "Arial", color = "#2E3440"
  ) +
  
  # Scales
  scale_y_continuous(
    limits = c(0, 14),
    breaks = seq(0, 12, 2),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  scale_fill_manual(
    values = c("Hyperforin" = "#2E3440", "Quercetin" = "#B0B0B0"),
    name = NULL
  ) +
  
  labs(
    title = "Hyperforin's advantage persists under expression weighting",
    subtitle = "The gap narrows but Hyperforin still leads in both methods",
    x = NULL,
    y = "Influence Z-score",
    caption = str_wrap(
      "[BIOLOGICAL VALIDATION] Expression-weighted influence (EWI) tests whether Hyperforin's advantage reflects liver-specific biology or network topology alone. GTEx v8 liver expression (TPM ≥1) weights edge transitions, channeling signal preferentially through liver-expressed neighbors. The gap narrows from +5.8 (RWR) to +3.5 (EWI) as Quercetin's high-expression targets (e.g., CFB at 1115 TPM) improve its signal propagation—yet Hyperforin's lead persists. Data: STRING v12.0 (≥900), n = 1,000 degree-matched permutations per method.",
      width = 105
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
    axis.text.x = element_text(face = "bold", size = 10, lineheight = 0.9),
    axis.title.y = element_text(face = "bold", size = 11, margin = margin(r = 8)),
    
    # Titles
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#4A4A4A", margin = margin(b = 15)),
    plot.caption = element_text(size = 8.5, hjust = 0, color = "#5A5A5A", 
                                lineheight = 1.3, margin = margin(t = 12)),
    
    # Legend
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 10, face = "bold"),
    legend.key.size = unit(0.5, "cm"),
    legend.margin = margin(b = 10),
    
    # Grid
    panel.grid.major.y = element_line(color = "#F0F0F0", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    
    # Margins
    plot.margin = margin(15, 20, 15, 15)
  )

# --- Display & Save ---
print(p)

ggsave(
  filename = here("figures", "main", "fig7_ewi_comparison.pdf"),
  plot = p,
  width = 180, height = 130, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = here("figures", "main", "fig7_ewi_comparison.tiff"),
  plot = p,
  width = 180, height = 130, units = "mm",
  dpi = 300, compression = "lzw"
)

message("✓ Figure 7 (EWI Comparison - Paired Bar Chart) saved")
