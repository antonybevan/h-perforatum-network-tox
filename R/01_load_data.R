# ============================================================================
# 01_load_data.R - Load and Validate All Result Tables
# H. perforatum Network Toxicology Analysis
# ============================================================================

# --- Dependencies ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, readr)

# --- Define Expected Files ---
required_files <- c(
  "standard_rwr_lcc_permutation_results.csv",
  "expression_weighted_rwr_permutation_results.csv",
  "shortest_path_permutation_results.csv",
  "bootstrap_summary.csv",
  "chemical_similarity_summary.csv",
  "consolidated_results.csv"
)

# --- Data Integrity Check ---
message("\n--- Data Integrity Check ---")
missing_files <- c()
for (f in required_files) {
  path <- here("results", "tables", f)
  if (!file.exists(path)) {
    missing_files <- c(missing_files, f)
    message(paste("✗ Missing:", f))
  } else {
    message(paste("✓ Found:", f))
  }
}

if (length(missing_files) > 0) {
  stop(paste("\n❌ Missing required files:", paste(missing_files, collapse = ", "),
             "\nRun the analysis pipeline first."))
}

# --- Load Result Tables ---
message("\n--- Loading Data ---")
rwr_results <- read_csv(here("results", "tables", "standard_rwr_lcc_permutation_results.csv"), show_col_types = FALSE)
ewr_results <- read_csv(here("results", "tables", "expression_weighted_rwr_permutation_results.csv"), show_col_types = FALSE)
sp_results <- read_csv(here("results", "tables", "shortest_path_permutation_results.csv"), show_col_types = FALSE)
bootstrap_results <- read_csv(here("results", "tables", "bootstrap_summary.csv"), show_col_types = FALSE)
bootstrap_iter <- read_csv(here("results", "bootstrap_sensitivity.csv"), show_col_types = FALSE)
chemsim_results <- read_csv(here("results", "tables", "chemical_similarity_summary.csv"), show_col_types = FALSE)
consolidated <- read_csv(here("results", "tables", "consolidated_results.csv"), show_col_types = FALSE)

# --- Validate Data Structure ---
message("\n--- Validating Data Structure ---")

# Check RWR results
stopifnot("compound" %in% names(rwr_results))
stopifnot("observed_influence" %in% names(rwr_results))
stopifnot("z_score" %in% names(rwr_results))
stopifnot(nrow(rwr_results) > 0)
message("✓ RWR results validated")

# Check EWR results
stopifnot("compound" %in% names(ewr_results))
stopifnot("observed_influence" %in% names(ewr_results))
stopifnot("z_score" %in% names(ewr_results))
stopifnot(nrow(ewr_results) > 0)
message("✓ EWR results validated")

# Check Shortest Path results
stopifnot("compound" %in% names(sp_results))
stopifnot("z_score" %in% names(sp_results))
stopifnot(nrow(sp_results) > 0)
message("✓ Shortest path results validated")

# Check Bootstrap results
stopifnot("observed_influence" %in% names(bootstrap_results))
stopifnot("bootstrap_mean" %in% names(bootstrap_results))
stopifnot(nrow(bootstrap_iter) >= 100)
message("✓ Bootstrap results validated")

# Check Chemical Similarity results
stopifnot("compound" %in% names(chemsim_results))
stopifnot("max_sim_DILI_positive" %in% names(chemsim_results))
stopifnot(nrow(chemsim_results) == 2)
message("✓ Chemical similarity results validated")

# --- Filter to STRING ≥900 (primary analysis) ---
rwr_900 <- rwr_results %>% filter(network_threshold == 900)
ewr_900 <- ewr_results %>% filter(network_threshold == 900)
sp_900 <- sp_results %>% filter(network_threshold == 900)

stopifnot(nrow(rwr_900) == 2)  # Hyperforin + Quercetin
stopifnot(nrow(ewr_900) == 2)
stopifnot(nrow(sp_900) == 2)
message("✓ STRING ≥900 filter applied (2 compounds each)")

# --- Create Master Summary Table ---
master_summary <- tibble(
  Compound = c("Hyperforin", "Quercetin"),
  Targets = c(10, 62),
  RWI_Z = c(
    rwr_900 %>% filter(compound == "Hyperforin") %>% pull(z_score),
    rwr_900 %>% filter(compound == "Quercetin") %>% pull(z_score)
  ),
  EWI_Z = c(
    ewr_900 %>% filter(compound == "Hyperforin") %>% pull(z_score),
    ewr_900 %>% filter(compound == "Quercetin") %>% pull(z_score)
  ),
  SP_Z = c(
    sp_900 %>% filter(compound == "Hyperforin") %>% pull(z_score),
    sp_900 %>% filter(compound == "Quercetin") %>% pull(z_score)
  ),
  RWI_Influence = c(
    rwr_900 %>% filter(compound == "Hyperforin") %>% pull(observed_influence),
    rwr_900 %>% filter(compound == "Quercetin") %>% pull(observed_influence)
  ),
  EWI_Influence = c(
    ewr_900 %>% filter(compound == "Hyperforin") %>% pull(observed_influence),
    ewr_900 %>% filter(compound == "Quercetin") %>% pull(observed_influence)
  ),
  Efficiency_RWR = RWI_Influence / Targets,
  Efficiency_EWI = EWI_Influence / Targets
)

# --- Validate Key Results ---
message("\n--- Validating Key Results ---")
stopifnot(master_summary$RWI_Z[1] > master_summary$RWI_Z[2])  # Hyperforin > Quercetin
stopifnot(master_summary$EWI_Z[1] > master_summary$EWI_Z[2])
stopifnot(master_summary$Efficiency_RWR[1] > master_summary$Efficiency_RWR[2])
message("✓ Hyperforin dominates across all metrics")

# --- Print Summary ---
message("\n--- Master Summary (STRING ≥900) ---")
print(master_summary)

message("\n✓ Data loaded and validated successfully")
message("✓ Ready for figure generation")
