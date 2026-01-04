# ============================================================================
# fig5_bootstrap.R - Figure 5: Bootstrap Sensitivity Analysis
# Defensive figure: Is Hyperforin's advantage just a target-count artefact?
# Distribution plot + observed reference (Option A: Density + Vertical Line)
# ============================================================================

source("R/00_setup_pub.R")
source("R/01_load_data.R")

# --- Load Bootstrap Data ---
bootstrap_iter <- read_csv(here("results", "bootstrap_sensitivity.csv"))
bootstrap_summary <- read_csv(here("results", "tables", "bootstrap_summary.csv"))

# --- Extract Key Values ---
hyperforin_observed <- bootstrap_summary$observed_influence[1]
ci_lower <- bootstrap_summary$ci_95_lower[1]
ci_upper <- bootstrap_summary$ci_95_upper[1]
bootstrap_mean <- bootstrap_summary$bootstrap_mean[1]

# --- Create Plot ---
p <- ggplot(bootstrap_iter, aes(x = quercetin_sampled_influence)) +
  
  # 95% CI shaded region

  annotate(
    "rect",
    xmin = ci_lower, xmax = ci_upper,
    ymin = -Inf, ymax = Inf,
    fill = "#E8E8E8", alpha = 0.6
  ) +
  
  # Density curve
  geom_density(
    fill = "#B0B0B0", color = "#4A4A4A",
    alpha = 0.7, linewidth = 0.6
  ) +
  
  # Bootstrap mean (dashed line)
  geom_vline(
    xintercept = bootstrap_mean,
    linetype = "dashed", color = "#666666", linewidth = 0.6
  ) +
  
  # Hyperforin observed (solid bold line)
  geom_vline(
    xintercept = hyperforin_observed,
    linetype = "solid", color = "#2E3440", linewidth = 1.2
  ) +
  
  # Label for Hyperforin
  annotate(
    "text",
    x = hyperforin_observed, y = Inf,
    label = sprintf("Hyperforin\nobserved = %.3f", hyperforin_observed),
    hjust = 1.1, vjust = 1.5,
    size = 3.5, fontface = "bold", family = "Arial", color = "#2E3440",
    lineheight = 0.9
  ) +
  
  # Label for 95% CI (at bottom of shaded region)
  annotate(
    "label",
    x = ci_upper - 0.01, y = 2,
    label = "95% CI",
    hjust = 0.5, vjust = 0,
    size = 3, fontface = "bold", family = "Arial", color = "#505050",
    fill = "white", label.size = 0
  ) +
  
  # Scales
  scale_x_continuous(
    limits = c(0, 0.13),
    breaks = seq(0, 0.12, 0.02),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1))
  ) +
  
  labs(
    title = "Bootstrap sensitivity analysis excludes target-count confounding",
    subtitle = "Distribution of influence scores from random 10-target samples",
    x = "Influence score (RWR)",
    y = "Density",
    caption = str_wrap(
      "[ROBUSTNESS CONTROL] Bootstrap sensitivity: 100 random 10-target samples from Quercetin's pool, scored by random walk with restart (RWI). Shaded = 95% CI. Hyperforin (solid line) exceeds entire distribution. Data: STRING v12.0 (≥900).",
      width = 90
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
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#4A4A4A", margin = margin(b = 15)),
    plot.caption = element_text(size = 8.5, hjust = 0, color = "#5A5A5A", 
                                lineheight = 1.3, margin = margin(t = 12)),
    
    # Grid
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Margins
    plot.margin = margin(15, 15, 15, 25)
  )

# --- Display & Save ---
print(p)

ggsave(
  filename = here("figures", "main", "fig5_bootstrap.pdf"),
  plot = p,
  width = 160, height = 120, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = here("figures", "main", "fig5_bootstrap.tiff"),
  plot = p,
  width = 160, height = 120, units = "mm",
  dpi = 300, compression = "lzw"
)

message("✓ Figure 5 (Bootstrap Sensitivity) saved - Defensive hypothesis killer")
