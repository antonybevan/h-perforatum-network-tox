# ============================================================================
# fig8_waterfall.R - Figure 8: Waterfall Chart - Gap Decomposition
# Shows WHY the gap narrows: expression effect on each compound
# ============================================================================

source("R/00_setup_pub.R")
source("R/01_load_data.R")

# --- Extract Z-scores ---
hyp_rwr <- rwr_900 %>% filter(compound == "Hyperforin") %>% pull(z_score)
quer_rwr <- rwr_900 %>% filter(compound == "Quercetin") %>% pull(z_score)
hyp_ewi <- ewr_900 %>% filter(compound == "Hyperforin") %>% pull(z_score)
quer_ewi <- ewr_900 %>% filter(compound == "Quercetin") %>% pull(z_score)

# --- Calculate Components ---
rwr_gap <- hyp_rwr - quer_rwr  # Starting gap
hyp_change <- hyp_ewi - hyp_rwr  # Hyperforin's change (negative = lost)
quer_change <- quer_ewi - quer_rwr  # Quercetin's change (positive = gained)
ewi_gap <- hyp_ewi - quer_ewi  # Final gap

# --- Waterfall Data ---
waterfall_data <- tibble(
  Step = c(
    "RWR advantage",
    "Hyperforin\nexpression effect",
    "Quercetin\nexpression effect",
    "EWI advantage"
  ),
  Value = c(rwr_gap, hyp_change, -quer_change, ewi_gap),
  Type = c("start", "change", "change", "end"),
  # Calculate cumulative positions for waterfall
  End = cumsum(c(rwr_gap, hyp_change, -quer_change, 0)),
  Start = c(0, rwr_gap, rwr_gap + hyp_change, 0)
)

# Adjust the last bar to show final value from 0
waterfall_data$Start[4] <- 0
waterfall_data$End[4] <- ewi_gap

# Set factor order
waterfall_data <- waterfall_data %>%
  mutate(
    Step = factor(Step, levels = Step),
    Fill = case_when(
      Type == "start" ~ "Initial",
      Type == "end" ~ "Final",
      Value < 0 ~ "Decrease",
      TRUE ~ "Increase"
    )
  )

# --- Waterfall Chart ---
p <- ggplot(waterfall_data, aes(x = Step)) +
  
  # Connecting lines between bars
  geom_segment(
    data = waterfall_data %>% filter(Type != "end"),
    aes(x = as.numeric(Step) + 0.4, xend = as.numeric(Step) + 0.6,
        y = End, yend = End),
    color = "#888888", linewidth = 0.5, linetype = "dashed"
  ) +
  
  # Bars
  geom_rect(
    aes(xmin = as.numeric(Step) - 0.4, xmax = as.numeric(Step) + 0.4,
        ymin = Start, ymax = End, fill = Fill),
    color = "#2E3440", linewidth = 0.4
  ) +
  
  # Value labels
  geom_text(
    aes(y = (Start + End) / 2,
        label = sprintf("%+.1f", ifelse(Type == "change", Value, End - Start))),
    size = 4, fontface = "bold", family = "Arial", color = "white"
  ) +
  
  # Zero line
  geom_hline(yintercept = 0, color = "#2E3440", linewidth = 0.6) +
  
  # Scales
  scale_y_continuous(
    limits = c(-3, 7),
    breaks = seq(-2, 6, 2),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  scale_fill_manual(
    values = c(
      "Initial" = "#2E3440",
      "Final" = "#2E3440",
      "Decrease" = "#D64545",  # Red for loss
      "Increase" = "#45A845"   # Green for gain (but we subtract it, so it's a loss to gap)
    ),
    guide = "none"
  ) +
  
  labs(
    title = "Decomposition: Why does the gap narrow under EWI?",
    subtitle = "Hyperforin loses slightly; Quercetin gains—but the lead persists",
    x = NULL,
    y = "Contribution to Z-score gap",
    caption = str_wrap(
      "[MECHANISTIC INSIGHT] The RWR advantage (+5.8) decomposes under expression weighting: Hyperforin's high-expression targets (CYP2C9, CYP3A4) provide less relative advantage when the entire network is expression-weighted. Quercetin's CFB (1115 TPM) boosts its signal. Net effect: gap narrows to +3.5 but never closes. GTEx v8 liver expression (TPM ≥1). Data: STRING v12.0 (≥900).",
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
    axis.text.x = element_text(face = "bold", size = 9, lineheight = 0.9),
    axis.title.y = element_text(face = "bold", size = 11, margin = margin(r = 8)),
    
    # Titles
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#4A4A4A", margin = margin(b = 15)),
    plot.caption = element_text(size = 8.5, hjust = 0, color = "#5A5A5A", 
                                lineheight = 1.3, margin = margin(t = 12)),
    
    # Grid
    panel.grid.major.y = element_line(color = "#F0F0F0", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    
    # Margins
    plot.margin = margin(15, 20, 15, 15)
  )

# --- Display & Save ---
print(p)

ggsave(
  filename = here("figures", "main", "fig8_waterfall.pdf"),
  plot = p,
  width = 180, height = 130, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = here("figures", "main", "fig8_waterfall.tiff"),
  plot = p,
  width = 180, height = 130, units = "mm",
  dpi = 300, compression = "lzw"
)

message("✓ Figure 8 (Waterfall Chart - Gap Decomposition) saved")
