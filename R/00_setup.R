# ============================================================================
# 00_setup.R - Lancet-Quality Visualization Setup
# H. perforatum Network Toxicology Analysis
# ============================================================================

# --- Package Installation ---
required_packages <- c(
  # Core visualization
  "ggplot2", "ggpubr", "patchwork", "scales",
  
  # Network visualization  
  "igraph", "ggraph", "tidygraph",
  
  # Publication themes
  "ggsci", "ggthemes", "extrafont",
  
  # Statistical annotations
  "ggrepel", "ggsignif", "ggdist",
  
  # Data wrangling
  "tidyverse", "readr", "here",
  
  # Tables
  "gt", "kableExtra"
)

# Install missing packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

invisible(lapply(required_packages, install_if_missing))

# --- Lancet Theme ---
theme_lancet <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +  # Use default font
    theme(
      # Clean background
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "#EBEBEB", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      
      # Professional axes
      axis.line = element_line(color = "#333333", linewidth = 0.5),
      axis.ticks = element_line(color = "#333333", linewidth = 0.3),
      axis.text = element_text(color = "#333333", size = 10),
      axis.title = element_text(color = "#000000", size = 11, face = "bold"),
      
      # Clean legend
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      
      # Title styling
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 11, color = "#666666"),
      plot.caption = element_text(size = 9, color = "#999999", hjust = 1),
      
      # Margins
      plot.margin = margin(15, 15, 15, 15)
    )
}

# --- Lancet Color Palette ---
lancet_colors <- c(
  "Hyperforin" = "#00468B",   # Deep blue
  "Quercetin" = "#ED0000",    # Lancet red
  "DILI" = "#42B540",         # Green
  "Null" = "#ADB6B6",         # Gray
  "Significant" = "#0099B4",  # Teal
  "NS" = "#925E9F"            # Purple
)

# Scale functions
scale_fill_lancet <- function(...) {
  scale_fill_manual(values = lancet_colors, ...)
}

scale_color_lancet <- function(...) {
  scale_color_manual(values = lancet_colors, ...)
}

# --- Export Settings ---
fig_width <- 180  # mm (journal column width)
fig_height <- 120 # mm
fig_dpi <- 300    # print quality

save_figure <- function(plot, filename, width = fig_width, height = fig_height) {
  ggsave(
    filename = here::here("figures", "main", paste0(filename, ".pdf")),
    plot = plot,
    width = width, height = height, units = "mm", dpi = fig_dpi
  )
  ggsave(
    filename = here::here("figures", "main", paste0(filename, ".png")),
    plot = plot,
    width = width, height = height, units = "mm", dpi = 600
  )
  message(paste("Saved:", filename))
}

# --- Project Paths ---
data_dir <- here::here("results", "tables")
fig_dir <- here::here("figures", "main")

message("✓ Setup complete. Run 01_load_data.R next.")
