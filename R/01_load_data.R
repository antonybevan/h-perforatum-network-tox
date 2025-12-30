# ============================================================================
# 01_load_data.R - Load All Result Tables
# ============================================================================

source(here::here("R", "00_setup.R"))

# --- Load Result Tables ---
rwr_results <- read_csv(here::here("results", "tables", "standard_rwr_lcc_permutation_results.csv"))
ewr_results <- read_csv(here::here("results", "tables", "expression_weighted_rwr_permutation_results.csv"))
sp_results <- read_csv(here::here("results", "tables", "shortest_path_permutation_results.csv"))
bootstrap_results <- read_csv(here::here("results", "tables", "bootstrap_summary.csv"))
chemsim_results <- read_csv(here::here("results", "tables", "chemical_similarity_summary.csv"))
consolidated <- read_csv(here::here("results", "tables", "consolidated_results.csv"))

# --- Filter to STRING ≥900 (primary analysis) ---
rwr_900 <- rwr_results %>% filter(network_threshold == 900)
ewr_900 <- ewr_results %>% filter(network_threshold == 900)
sp_900 <- sp_results %>% filter(network_threshold == 900)

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
  PTNI_RWI = c(
    rwr_900 %>% filter(compound == "Hyperforin") %>% mutate(ptni = observed_influence/n_targets) %>% pull(ptni),
    rwr_900 %>% filter(compound == "Quercetin") %>% mutate(ptni = observed_influence/n_targets) %>% pull(ptni)
  )
)

# --- Print Summary ---
message("\n✓ Data loaded successfully")
message("\n--- Master Summary (STRING ≥900) ---")
print(master_summary)

message("\n✓ Ready for figure generation. Run 02_fig_violin.R next.")
