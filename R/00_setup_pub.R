# ============================================================================
# 00_setup_pub.R - Strict Lancet/Q1 Visualization Setup
# H. perforatum Network Toxicology Analysis
# ============================================================================

# --- Package Management ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  # Core visualization
  ggplot2, ggpubr, ggsci, scales,
  
  # Text and labels
  ggrepel, extrafont,
  
  # Layout and composition
  grid, gridExtra, cowplot, patchwork,
  
  # Data wrangling
  tidyverse, here, readr,
  
  # Tables
  gt, kableExtra
)

# --- Font Setup ---
# Lancet uses Arial/Helvetica. Register system fonts on Windows.
tryCatch({
  loadfonts(device = "win", quiet = TRUE)
}, error = function(e) {
  message("Note: Font registration skipped (non-Windows or fonts already loaded)")
})
main_font <- "Arial"

# --- Lancet Color Palette ---
# True Lancet colors from ggsci
lancet_cols <- pal_lancet("lanonc")(9)
cols <- c(
  "Hyperforin" = lancet_cols[1],  # Deep Blue (#00468B)
  "Quercetin"  = lancet_cols[2],  # Deep Red (#ED0000)
  "DILI"       = lancet_cols[3],  # Green (#42B540)
  "Neutral"    = lancet_cols[9]   # Grey
)

# Scale functions for easy use
scale_fill_lancet <- function(...) scale_fill_manual(values = cols, ...)
scale_color_lancet <- function(...) scale_color_manual(values = cols, ...)

# --- Publication Theme ---
theme_lancet_pub <- function(base_size = 11, base_family = main_font) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # Plot Margins and Background
      plot.margin = margin(15, 15, 15, 15),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      
      # Axis Lines & Ticks (Thicker for print)
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      
      # Text Hierarchy
      plot.title = element_text(face = "bold", size = 14, hjust = 0, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 12, color = "#404040", margin = margin(b = 15)),
      plot.caption = element_text(size = 9, color = "#606060", hjust = 0, margin = margin(t = 15)),
      
      # Axis Text
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 10, color = "black"),
      
      # Legend (Clean and minimal)
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 10),
      legend.background = element_blank(),
      legend.key = element_blank(),
      
      # Grid (Minimal or horizontal only)
      panel.grid.major.y = element_line(color = "#E0E0E0", linewidth = 0.4),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Facets
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 11, hjust = 0)
    )
}

# --- Export Function (Publication Standard) ---
save_pub_plot <- function(plot, filename, w = 180, h = 150) {
  # Ensure output directory exists
  dir.create(here("figures", "main"), showWarnings = FALSE, recursive = TRUE)
  
  # 1. PDF (Vector high quality)
  ggsave(
    filename = here("figures", "main", paste0(filename, ".pdf")),
    plot = plot,
    width = w, height = h, units = "mm",
    device = cairo_pdf
  )
  
  # 2. TIFF (300 DPI Raster for submission)
  ggsave(
    filename = here("figures", "main", paste0(filename, ".tiff")),
    plot = plot,
    width = w, height = h, units = "mm",
    dpi = 300, compression = "lzw"
  )
  
  message(paste("✓ Saved", filename, "(PDF + 300dpi TIFF)"))
}

# --- Project Paths ---
data_dir <- here("results", "tables")
fig_dir <- here("figures", "main")

message("✓ Setup complete (00_setup_pub.R)")
