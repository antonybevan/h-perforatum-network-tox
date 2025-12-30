# ============================================================================
# 00_setup_pub.R - Strict Lancet/Q1 Visualization Theme
# ============================================================================

# --- Packages ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  ggplot2, ggpubr, ggsci, extrafont, 
  grid, gridExtra, cowplot, 
  tidyverse, here, scales
)

# --- Fonts ---
# Lancet uses Arial/Helvetica. On Windows, we need to register it.
loadfonts(device = "win", quiet = TRUE)
main_font <- "Arial"

# --- Lancet Color Palette (ggsci) ---
# True Lancet colors: Blue (#00468B), Red (#ED0000), Green (#42B540)
lancet_cols <- pal_lancet("lanonc")(9)
cols <- c(
  "Hyperforin" = lancet_cols[1],  # Deep Blue
  "Quercetin"  = lancet_cols[2],  # Deep Red
  "DILI"       = lancet_cols[3],  # Green
  "Neutral"    = lancet_cols[9]   # Grey
)

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
save_pub_plot <- function(plot, filename, w=180, h=150) {
  # 1. PDF (Vector high quality)
  ggsave(
    filename = here("figures", "main", paste0(filename, ".pdf")),
    plot = plot,
    width = w, height = h, units = "mm",
    device = cairo_pdf  # Superior font handling
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
