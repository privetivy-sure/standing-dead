#
# ==============================================================================
# BAYESIAN NON-LINEAR MIXED EFFECTS MODEL for lying deadwood
# ==============================================================================
# ==============================================================================


write.table(trees_fagus_dc1, 
            file = "trees_fagus_dc1.txt", 
            sep = "\t",          
            dec = ".",           
            row.names = FALSE, 
            quote = FALSE,       
            na = "",             
            fileEncoding = "UTF-8")



# FILTER ERROR-FREE and no dead RECORDS
# ==============================================================================

trees_clear_L <- trees_bayes_L %>%
  filter(
    errorlife_tree == FALSE & 
      discard_nodead == FALSE &
      (status3 != "A" | is.na(status3)) # Removes live stems (A)
  )



# ==============================================================================
# COMPUTING VOLUME
# ==============================================================================

# ------------------------------------------------------------------------------
# Methodological correction for LFRLP-FAWF: Move diameter_1 to diameter_middle
# ------------------------------------------------------------------------------

trees_clear_L <- trees_clear_L %>%
  mutate(
    # Create diameter_middle column and move diameter_1 into it for LFRLP-FAWF
    diameter_middle = case_when(
      institute == "LFRLP-FAWF" ~ diameter_1,
      TRUE ~ NA_real_
    ),
    
    # Clear diameter_1 for LFRLP-FAWF
    diameter_1 = case_when(
      institute == "LFRLP-FAWF" ~ NA_real_,
      TRUE ~ diameter_1
    )
  )


# ------------------------------------------------------------------------------
# Diameter type classification
# ------------------------------------------------------------------------------

trees_clear_L <- trees_clear_L %>%
  mutate(
    diam_type = case_when(
      # diameter_middle is present (highest priority, e.g., for LFRLP-FAWF)
      !is.na(diameter_middle) ~ "d_middle",
      
      # Both diameter_1 and diameter_2 are present (diameter_130 is ignored)
      !is.na(diameter_1) & !is.na(diameter_2) ~ "d_1_2",
      
      # Only diameter_1 is present (diameter_2 is missing, diameter_130 is ignored)
      !is.na(diameter_1) & is.na(diameter_2)  ~ "d_1_x",
      
      # diameter_1 is missing, but diameter_130 is present (regardless of diameter_2)
      is.na(diameter_1) & !is.na(diameter_130) ~ "d_x_x_130",
      
      # Both diameter_1 and diameter_130 are missing, but diameter_2 is present
      is.na(diameter_1) & is.na(diameter_130) & !is.na(diameter_2) ~ "d_x_2",
      
      # Catch-all for any other edge cases (e.g., all diameters missing)
      TRUE ~ NA_character_
    )
  )

library(dplyr)
library(tidyr)





# FILL MISSING VALUES DIAMETER_1 AND DIAMETER_2
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Fill missing values from previous inventory years (VUK logic)
# ------------------------------------------------------------------------------

trees_bayes_L3 <- trees_bayes_L3 %>%
  arrange(tree_unique_id, inventory_year) %>%
  group_by(tree_unique_id) %>%
  mutate(
    prev_diameter_1 = lag(diameter_1),
    prev_diameter_2 = lag(diameter_2)
  ) %>%
  ungroup()

# ------------------------------------------------------------------------------
# 2. Impute new diam_1 and diam_2 (Diameters in mm, length in meters)
# ------------------------------------------------------------------------------

trees_bayes_L3 <- trees_bayes_L3 %>%
  mutate(
    # Step A: Compute NEW column diam_1 [mm]
    diam_1b = case_when(
      # If original diameter_1 exists, keep it in the new column
      !is.na(diameter_1) ~ diameter_1,
      
      # Imputation rules when diameter_1 is missing:
      diam_type == "d_middle" ~ diameter_middle + (length * 5),
      diam_type == "d_x_x_130" & institute == "LFRLP-FAWF" ~ diameter_130 + 13,
      diam_type == "d_x_x_130" & institute == "VUK" & piece_id == "a" ~ diameter_130 + 13,
      diam_type == "d_x_x_130" & institute == "VUK" & piece_id != "a" ~ prev_diameter_1,
      
      TRUE ~ NA_real_
    ),
    
    # Step B: Compute NEW column diam_2 [mm]
    diam_2 = case_when(
      # If original diameter_2 exists, keep it in the new column
      !is.na(diameter_2) ~ diameter_2,
      
      # Imputation rules when diameter_2 is missing:
      diam_type == "d_middle" ~ diameter_middle - (length * 5),
      diam_type == "d_x_x_130" & institute == "VUK" & piece_id != "a" & !is.na(prev_diameter_2) ~ prev_diameter_2,
      !is.na(diam_1) ~ diam_1 - (length * 10),
      
      TRUE ~ NA_real_
    ),
    
    # Step C: Apply safety threshold (values <= 10 mm set to 10 mm default)
    diam_2 = if_else(!is.na(diam_2) & diam_2 <= 10, 10, diam_2)
  ) %>%
  # Remove temporary lookup columns
  select(-prev_diameter_1, -prev_diameter_2)


# ------------------------------------------------------------------------------
# Calculate log volume (volume in m3) using Smalian's formula
# Inputs: diam_1 [mm], diam_2 [mm], length [m]
# ------------------------------------------------------------------------------

trees_bayes_L3 <- trees_bayes_L3 %>%
  mutate(
    # Smalian's formula for truncated cone section volume (m3)
    volume = (pi / 8) * (((diam_1 / 1000)^2) + ((diam_2 / 1000)^2)) * length
  )

library(dplyr)

trees_bayes_L3 <- trees_bayes_L3 %>%
  mutate(
    vol_fraver = (pi * length / 8) * ((diam_1 / 1000)^2 + (diam_2 / 1000)^2)
  )


# ==============================================================================
# DECAY STAGES
# ==============================================================================

# SUMMARY OF DECAY PER INSTITUTE
# ------------------------------------------------------------------------------
df_long <- trees_clear_L %>%
  distinct(institute, decay) %>%
  arrange(institute, decay)

library(dplyr)
library(tidyr)


# DECAY CLASS RECLASIFICATION
# ------------------------------------------------------------------------------

trees_bayes_L <- trees_bayes_L %>%
  mutate(
    decay_new = case_when(
      decay %in% c(11, 101) ~ 1,
      decay == 201         ~ 2,
      decay %in% c(12, 301) ~ 3,
      decay %in% c(13, 401) ~ 4,
      TRUE                  ~ decay  # Copy other decay (1-5)
    )
  )

# Summary table
decay_summary <- trees_bayes_L %>%
  group_by(institute, decay_new) %>%
  summarise(n = n(), .groups = "drop")

library(dplyr)

# Summary table for NA in decay
decay_na_summary <- trees_bayes_L %>%
  group_by(institute, site_name) %>%
  summarise(
    total_records = n(),                            
    na_count = sum(is.na(decay)),                    
    pct_na = round((na_count / total_records) * 100, 2), 
    .groups = "drop"
  ) %>%
  arrange(institute, site_name)





# ==============================================================================
# WOOD DENSITY, CARBON FRACTION, CARBON STOCK
# ==============================================================================

library(dplyr)
library(stringr)

trees_bayes_L <- trees_bayes_L %>%
  mutate(
    # Extraction of first word from full_scientific
    first_word = word(full_scientific, 1),
    
    genus = case_when(
      first_word %in% c("Alnus", "Carpinus", "Fagus", "Fraxinus", 
                        "Quercus", "Abies", "Picea", "Pinus") ~ first_word,
      
      # Specific cases
      full_scientific %in% c("Coniferous", "Unknown conifer") ~ "Conifer",
      full_scientific == "Unknown species"                   ~ "Unknown",
      full_scientific == "changed"                           ~ "",
      
      # All other values to Broadleave
      TRUE                                                  ~ "Broadleave"
    )
  ) %>%
  # Odstranění pomocného sloupce
  select(-first_word)



# 1. COMPUTING AVERAGE DEADWOOD DENSITY FOR UNKNOWN SPECIES
# -----------------------------------------------------------------------------
# Aggregate raw data to ensure each Genus x Decay_class combination appears at most once
deadwood_density_clean <- deadwood_density %>%
  group_by(Genus, Decay_class) %>%
  summarise(
    Deadwood_density_g_cm_ = mean(Deadwood_density_g_cm_, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate mean deadwood density across all genera for "Unknown" species per Decay_class
unknown_species <- deadwood_density_clean %>%
  group_by(Decay_class) %>%
  summarise(
    Deadwood_density_g_cm_ = mean(Deadwood_density_g_cm_, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Genus = "Unknown") 

# Merge clean species density data with the Unknown category
# distinct() prevents duplicate rows if the script is re-run in the same session
deadwood_density_prep <- bind_rows(deadwood_density_clean, unknown_species) %>%
  distinct(Genus, Decay_class, .keep_all = TRUE)


# 2. JOIN DEADWOOD DENSITY TO TREES TABLE
# ------------------------------------------------------------------------------
# Store initial row count before joining to detect unwanted cartesian expansion
rows_before <- nrow(trees_bayes_L1)

# Join deadwood density values based on genus and decay class
trees_bayes_L2 <- trees_bayes_L1 %>%
  left_join(
    deadwood_density_prep %>% select(Genus, Decay_class, Deadwood_density_g_cm_),
    by = c("genus" = "Genus", "decay_new" = "Decay_class")
  ) %>%
  rename(density2 = Deadwood_density_g_cm_)

# Verify row count consistency after the join operation
rows_after <- nrow(trees_bayes_L2)
cat("Row count before join:", rows_before, "| After join:", rows_after, "\n")

if (rows_before != rows_after) {
  warning("WARNING: Row count changed during left_join! Check join keys in deadwood_density.")
}

# ASSIGN CARBON FRACTION TO TREES TABLE
# ------------------------------------------------------------------------------

library(dplyr)

trees_bayes_L3 <- trees_bayes_L2 %>%
  mutate(
    c_fraction = case_when(
      genus %in% c("Alnus", "Broadleave", "Carpinus", "Fagus", "Fraxinus", "Quercus") ~ 0.475,
      genus %in% c("Abies", "Conifer", "Picea", "Pinus") ~ 0.505,
      genus == "Unknown" ~ (0.475 + 0.505) / 2, # Hodnota 0.490
      TRUE ~ NA_real_ 
    )
  )

# COMPUTING CARBON QUANTITY
# ------------------------------------------------------------------------------

library(dplyr)

trees_bayes_L3 <- trees_bayes_L3 %>%
  mutate(
    # Computing carbon stock in t
    carbon_t = volume * density2 * c_fraction,
    
    # Carbon stock in kg
    carbon_kg = volume * density2 * c_fraction * 1000
  )


# DYNAMICS OF CARBON, VOLUME, DECAY (DECREASE OR INCREASE THROUGH TIME) IN %
# ------------------------------------------------------------------------------


# 1. Identify valid trees (exclude any tree_unique_id containing NA in carbon_kg, volume, or decay_new)
valid_trees <- trees_bayes_L3 %>%
  group_by(tree_unique_id) %>%
  filter(!any(is.na(carbon_kg)) & !any(is.na(volume)) & !any(is.na(decay_new))) %>%
  pull(tree_unique_id) %>%
  unique()

# 2. Sort records chronologically and calculate differences (diff) between consecutive inventory years
diff_data <- trees_bayes_L3 %>%
  # Retain only trees with complete time series across all evaluated variables
  filter(tree_unique_id %in% valid_trees) %>%
  # Sort chronologically by tree ID and inventory year
  arrange(tree_unique_id, inventory_year) %>%
  # Group by individual tree to calculate lagged differences
  group_by(tree_unique_id) %>%
  mutate(
    diff_carbon = carbon_kg - lag(carbon_kg),
    diff_volume = volume - lag(volume),
    diff_decay  = decay_new - lag(decay_new)
  ) %>%
  # Remove the first record of each tree (where lag produces NA)
  filter(!is.na(diff_carbon)) %>%
  ungroup()

# 3. Create a summary of "unnatural" transitions for carbon, volume, and decay class
summary_table2 <- diff_data %>%
  summarise(
    total_comparisons      = n(),
    
    # Unnatural Carbon dynamics (carbon increases despite decomposition)
    carbon_positive_count  = sum(diff_carbon > 0),
    carbon_positive_pct    = round(mean(diff_carbon > 0) * 100, 2),
    
    # Unnatural Volume dynamics (volume increases due to geometry measurement error)
    volume_positive_count  = sum(diff_volume > 0),
    volume_positive_pct    = round(mean(diff_volume > 0) * 100, 2),
    
    # Unnatural Decay dynamics (decay class decreases / log becomes "fresher")
    decay_negative_count   = sum(diff_decay < 0),
    decay_negative_pct     = round(mean(diff_decay < 0) * 100, 2),
    
    # Combined anomalies: Carbon increases while volume decreases or stays constant
    carbon_up_volume_down  = sum(diff_carbon > 0 & diff_volume <= 0)
  )



library(dplyr)
library(purrr)

# ==============================================================================
# FILTERING AND SORTING DATA FOR BAYES REGRESSION
# ==============================================================================
library(purrr)

trees_bayes_updated <- trees_bayes_L3 %>%
  mutate(
    species_group = case_when(
      # Extract Carpinus betulus into its own group
      tree_species == "Carpinus betulus" ~ "Carpinus",
      # Merge remaining broad_hard and broad_soft into broadleave
      species_group %in% c("broad_hard", "broad_soft") ~ "broadleave",
      # Keep all other species groups unchanged (e.g., Picea, Pinus, etc.)
      TRUE ~ species_group
    )
  )




# ==============================================================================
# COMPUTING BASIC ATTRIBUTES AND PARAMETERS FOR BAYES
# ==============================================================================

library(dplyr)



# COMPUTATION OF TIME, RELATIVE CARBON AND INITIAL DBH
# ==============================================================================


library(dplyr)
library(tidyr)


library(dplyr)
library(tidyr)


# IDENTIFY LAST INVENTORY YEAR PER SITE AND UNIQUE SITE-YEARS
# =========================================================================
site_max_years <- trees_bayes_updated %>%
  group_by(composed_site_id) %>%
  summarise(
    site_last_year = max(inventory_year, na.rm = TRUE),
    .groups = "drop"
  )

# Extract strictly unique combinations of site and inventory years
site_inventory_years <- trees_bayes_updated %>%
  distinct(composed_site_id, inventory_year)



# GENERATE ZERO-CARBON ROWS FOR VANISHED TREES (GHOST TREES)
# =========================================================================
missing_trees_zero_rows <- trees_bayes_updated %>%
  left_join(site_max_years, by = "composed_site_id") %>%
  group_by(tree_unique_id) %>%
  summarise(
    composed_site_id  = first(composed_site_id),
    species_group     = first(species_group),
    last_seen_year    = max(inventory_year),
    site_last_year    = first(site_last_year),
    first_decay_new = first(decay_new),
    first_diam_1      = first(diam_1),
    first_carbon_kg   = first(carbon_kg),
    .groups = "drop"
  ) %>%
  # Keep only trees that disappeared BEFORE the site's final inventory
  filter(last_seen_year < site_last_year) %>%
  
  # Join site inventory years with explicit many-to-many parameter to silence warning
  left_join(
    site_inventory_years, 
    by = "composed_site_id",
    relationship = "many-to-many"
  ) %>%
  filter(inventory_year > last_seen_year) %>%
  group_by(tree_unique_id) %>%
  slice_min(inventory_year, n = 1) %>% # Retain only the first year of absence
  ungroup() %>%
  
  # Construct the terminal zero-carbon record
  transmute(
    tree_unique_id,
    composed_site_id,
    species_group,
    inventory_year,
    decay_new = NA,              # Tree no longer exists as a log
    diam_1      = first_diam_1,    # Preserve initial dimension reference
    carbon_kg   = 0                # Critical zero value (100% loss)
  )



# MERGE DATASETS AND CALCULATE TREE-LEVEL BASELINES
# =========================================================================
trees_bayes <- bind_rows(trees_bayes_updated, missing_trees_zero_rows) %>%
  # Ensure strict chronological sequence per tree before calculating metrics
  arrange(tree_unique_id, inventory_year) %>%
  
  # Group records by individual tree to calculate longitudinal baseline variables
  group_by(tree_unique_id) %>%
  mutate(
    # Time elapsed (in years) since the first inventory visit (t = 0)
    decaying_time = inventory_year - min(inventory_year),
    
    # Initial carbon mass at t = 0 (serves as the 100% reference point)
    carbon_init = first(carbon_kg),
    
    # Proportion of carbon remaining (at t = 0, carbon_rel equals 1.0)
    carbon_rel = carbon_kg / carbon_init,
    
    # Baseline diameter recorded during the first visit
    diam_1_init = first(diam_1),
    
    # Initial decay class at t = 0
    decay_init = first(decay_new)
  ) %>%
  ungroup()




# ==============================================================================
# CLEANING
# ==============================================================================

library(dplyr)

# Summary of removed rows (carbon_kg NA)
removed_rows_summary <- trees_bayes_final %>%
  filter(is.na(carbon_kg)) %>%
  count(composed_site_id, name = "removed_rows_count")

# Cleaning main table:
#    - remove columns vol_fraver, diam_1b, carbon_t
#    - remove rows with NA in carbon_kg
trees_bayes_fin <- trees_bayes_final %>%
  select(-any_of(c("vol_fraver", "diam_1b", "carbon_t"))) %>%
  filter(!is.na(carbon_kg))

library(dplyr)


# IDENTIFY AND REMOVE TREES FIRST MEASURED IN THE LAST SITE INVENTORY
# =========================================================================

# 1. Identify trees whose FIRST observation occurred in the LAST inventory year of their site
trees_first_seen_in_last_inv <- trees_bayes_fin %>%
  # Determine site-level final inventory year
  group_by(composed_site_id) %>%
  mutate(site_last_year = max(inventory_year, na.rm = TRUE)) %>%
  
  # Determine tree-level first observation year
  group_by(tree_unique_id) %>%
  mutate(tree_first_year = min(inventory_year)) %>%
  
  # Filter trees where first observation matches the site's final census year
  filter(tree_first_year == site_last_year) %>%
  ungroup()

# Summary count of excluded trees
n_removed_trees <- n_distinct(trees_first_seen_in_last_inv$tree_unique_id)
message("ℹ️ Number of trees first recorded in the last site inventory: ", n_removed_trees)


# 2. Filter out these single-visit terminal trees from the main dataset
trees_bayes_clean <- trees_bayes_fin %>%
  group_by(composed_site_id) %>%
  mutate(site_last_year = max(inventory_year, na.rm = TRUE)) %>%
  group_by(tree_unique_id) %>%
  mutate(tree_first_year = min(inventory_year)) %>%
  
  # Keep ONLY trees that were first observed BEFORE the final site census
  filter(tree_first_year < site_last_year) %>%
  
  # Clean up temporary helper columns
  select(-site_last_year, -tree_first_year) %>%
  ungroup()


library(dplyr)


# ADD HELPER COLUMNS FOR FIRST INVENTORY & INITIAL DECAY CLASS
# -----------------------------------------------------------------------------
trees_bayes_clean <- trees_bayes_clean %>%
  group_by(tree_unique_id) %>%
  mutate(
    # Identify the exact calendar year of the tree's first visit
    first_inventory_year = min(inventory_year),
    
    # Flag indicating whether a row belongs to the first inventory visit
    is_first_inventory = (inventory_year == first_inventory_year)
    
  ) %>%
  ungroup()




library(ggplot2)
library(dplyr)

# ==============================================================================
# FIRST ROUND BAYESIAN MODEL - FALLEN TREES RECORDED WITH DC 1
# ==============================================================================


# EXPLORATORY DATA ANALYSIS: FAGUS CARBON DECAY TRAJECTORY
# ==============================================================================

# FILTER FALLEN TREES OF SELECTED DECAY CLASS (1) FOR SPECIES
# ------------------------------------------------------------------------------

# List of sites to exclude due to low log counts
sites_to_remove <- c(
  "LFRLP-FAWF__203__Hollaenderschlag__NA", 
  "LFRLP-FAWF__232__Wuesttal__NA", 
  "LFRLP-FAWF__240__Schwappelbruch__NA", 
  "LFRLP-FAWF__246__Adelsberg__1",
  "URK__5__Zarnowka - Babia Gora N.P.__NA",
  "VUK__22__Velka Ples__NA", 
  "VUK__4__Bila Opava__NA"
)

trees_fagus_dc1 <- trees_bayes_ready %>%
  filter(species_group == "Fagus") %>%
  # Filter out low-sample sites
  filter(!composed_site_id %in% sites_to_remove) %>%
  filter(
    !is.na(carbon_rel),
    !is.na(decaying_time),
    !is.na(diam_1_init),
    !is.na(composed_site_id),
    !is.na(decay_init)
  ) %>%
  # SELECT ONLY LOGS ENTERING AT DECAY STAGE 1 (freshly fallen logs)
  filter(decay_init == 1 | decay_init == "1") %>%
  mutate(
    # Create new site variable with merged Žofín sub-sites
    site_merged = case_when(
      composed_site_id %in% c("VUK__1__Zofin__a", "VUK__1__Zofin__c") ~ "VUK__1__Zofin__ac",
      composed_site_id %in% c("VUK__1__Zofin__b", "VUK__1__Zofin__d") ~ "VUK__1__Zofin__bd",
      TRUE ~ composed_site_id  # Keep all other sites unchanged
    ),
    
    decay_init = as.factor(decay_init),
    carbon_rel = ifelse(carbon_rel >= 1, 0.9999, carbon_rel),
    # Standardize diameter (z-score) specifically for the dc1 subset
    diam_z = as.numeric(scale(diam_1_init))
  )


# OVERALL TRAJECTORY PLOT WITH GAM CURVE
# ==============================================================================

# Create scatter plot with jittering to resolve point overplotting at y = 1 and y = 0
ggplot(trees_fagus_class1, aes(x = decaying_time, y = carbon_rel)) +
  # Add slightly transparent points with minor vertical/horizontal jitter
  geom_jitter(
    alpha = 0.35, 
    width = 0.25, 
    height = 0.01, 
    color = "#2c7fb8"
  ) +
  
  # Fit a non-parametric GAM smooth trend to visualize empirical trajectory
  geom_smooth(
    method = "gam", 
    formula = y ~ s(x, bs = "cs"), 
    color = "#111111", 
    fill = "#9ebcda",
    se = TRUE
  ) +
  
  # Format axes boundaries and labels
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Relative Carbon Loss Trajectory for Beech (Fagus)",
    subtitle = "Trees initially recorded in Decay Class 1; black curve shows non-parametric GAM trend",
    x = "Decaying Time [Years] (decaying_time)",
    y = "Relative Carbon Fraction (carbon_rel)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )



saveRDS(fit_fagus_fixed, file = "fit_fagus_exp_fixed.rds")
fit_fagus_fixed <- readRDS(file = "fit_fagus_exp_fixed.rds")

# ==============================================================================
# BAYESIAN ZERO-INFLATED BETA DECAY MODELS: TESTING SITE-EFFECT FRAMEWORKS
# ==============================================================================
# Comparing three spatial frameworks:
# 1. Partial Pooling (Hierarchical / Random Intercepts)
# 2. No-Pooling (Independent Fixed Site Effects)
# 3. No-Site / Global Model (Universal Population Trajectory)
# ==============================================================================

library(brms)

# Shared sampling configuration
MCMC_CONTROL <- list(
  chains  = 3,
  cores   = 3,
  iter    = 2000,
  warmup  = 1000,
  seed    = 123,
  init    = 0,
  control = list(adapt_delta = 0.90)
)


# ------------------------------------------------------------------------------
# MODEL 1: PARTIAL POOLING (HIERARCHICAL / RANDOM INTERCEPTS)
# ------------------------------------------------------------------------------
# Shrinks site-level estimates toward a common population mean.
# Useful for stabilizing estimates across sites with varying sample sizes.

fit_fagus_partial_pooling <- brm(
  formula = bf(
    carbon_rel ~ exp(- exp(k) * decaying_time),
    
    # Hierarchical random intercept for site
    k ~ 1 + diam_z + (1 | site_merged),
    
    # Zero-inflation increases over time with site variation
    zi ~ 1 + decaying_time + (1 | site_merged), 
    
    nl = TRUE
  ),
  data   = trees_fagus_dc1,
  family = zero_inflated_beta(),
  prior  = c(
    prior(normal(-1.44, 0.5), nlpar = "k", class = "b", coef = "Intercept"),
    prior(normal(0, 0.5), nlpar = "k", class = "b", coef = "diam_z"),
    prior(normal(-3.5, 1), class = "Intercept", dpar = "zi"),
    prior(normal(0.08, 0.04), class = "b", coef = "decaying_time", dpar = "zi")
  ),
  chains  = MCMC_CONTROL$chains,
  cores   = MCMC_CONTROL$cores,
  iter    = MCMC_CONTROL$iter,
  warmup  = MCMC_CONTROL$warmup,
  seed    = MCMC_CONTROL$seed,
  init    = MCMC_CONTROL$init,
  control = MCMC_CONTROL$control
)


# ------------------------------------------------------------------------------
# MODEL 2: NO-POOLING (INDEPENDENT FIXED SITE EFFECTS)
# ------------------------------------------------------------------------------
# Models each site completely independently (`0 + site_merged`).
# Prevents shrinkage bias caused by high-sample-size sites (e.g., Žofín).

fit_fagus_no_pooling <- brm(
  formula = bf(
    carbon_rel ~ exp(- exp(k) * decaying_time),
    
    # Independent fixed intercepts per site (no global intercept)
    k ~ 0 + site_merged + diam_z,
    
    # Independent zero-inflation intercepts per site
    zi ~ 0 + site_merged + decaying_time, 
    
    nl = TRUE
  ),
  data   = trees_fagus_dc1,
  family = zero_inflated_beta(),
  prior  = c(
    # Independent priors for each site's log-k parameter
    prior(normal(-1.44, 1.0), nlpar = "k", class = "b"),
    
    # Prior for diameter moderation effect
    prior(normal(0, 0.5), nlpar = "k", class = "b", coef = "diam_z"),
    
    # Priors for site-specific initial zero-inflation and slope over time
    prior(normal(-3.5, 1.5), class = "b", dpar = "zi"),
    prior(normal(0.08, 0.04), class = "b", coef = "decaying_time", dpar = "zi")
  ),
  chains  = MCMC_CONTROL$chains,
  cores   = MCMC_CONTROL$cores,
  iter    = MCMC_CONTROL$iter,
  warmup  = MCMC_CONTROL$warmup,
  seed    = MCMC_CONTROL$seed,
  init    = MCMC_CONTROL$init,
  control = MCMC_CONTROL$control
)


# ------------------------------------------------------------------------------
# MODEL 3: NO-SITE / GLOBAL MODEL (UNIVERSAL POPULATION TRAJECTORY)
# ------------------------------------------------------------------------------
# Omits the site variable entirely.
# Generates a single universal equation for external plot-level modeling,
# dominated by long-term high-power dataset sources.

fit_fagus_no_site <- brm(
  formula = bf(
    carbon_rel ~ exp(- exp(k) * decaying_time),
    
    # Global decomposition rate as a function of diameter only
    k ~ 1 + diam_z,
    
    # Global zero-inflation probability increasing over time
    zi ~ 1 + decaying_time, 
    
    nl = TRUE
  ),
  data   = trees_fagus_dc1,
  family = zero_inflated_beta(),
  prior  = c(
    prior(normal(-1.44, 0.5), nlpar = "k", class = "b", coef = "Intercept"),
    prior(normal(0, 0.5), nlpar = "k", class = "b", coef = "diam_z"),
    prior(normal(-3.5, 1), class = "Intercept", dpar = "zi"),
    prior(normal(0.08, 0.04), class = "b", coef = "decaying_time", dpar = "zi")
  ),
  chains  = MCMC_CONTROL$chains,
  cores   = MCMC_CONTROL$cores,
  iter    = MCMC_CONTROL$iter,
  warmup  = MCMC_CONTROL$warmup,
  seed    = MCMC_CONTROL$seed,
  init    = MCMC_CONTROL$init,
  control = MCMC_CONTROL$control
)


# ==============================================================================
# MODEL COMPARISON (LOO-CV)
# ==============================================================================
# Compare predictive performance using Leave-One-Out Cross-Validation

loo_partial  <- loo(fit_fagus_partial_pooling)
loo_nopool   <- loo(fit_fagus_no_pooling)
loo_nosite   <- loo(fit_fagus_no_site)

# Compare expected log predictive density (ELPD)
loo_compare(loo_partial, loo_nopool, loo_nosite)




# Display model summary in console
summary(fit_fagus_fixed)

# Comparison of sites between themselves
ranef(fit_fagus_fixed)
# k_intercept for sites
coef(fit_fagus_exp_diam)$site_merged

# ==============================================================================
# COMPUTING DECOMPOSITION TIME FROM NEGATIVE-EXPONENTIAL KOEF. k
# ==============================================================================

library(dplyr)
library(tibble)

# 1. Extract k_Intercept coefficients for each individual site
k_table <- coef(fit_fagus_fixed)$site_merged[, , "k_Intercept"] %>%
  as.data.frame() %>%
  rownames_to_column(var = "site_merged") %>%
  mutate(
    # Real decomposition rate parameter k (for mean log diameter, diam_z = 0)
    k_real = exp(Estimate),
    k_real_lower = exp(Q2.5),
    k_real_upper = exp(Q97.5),
    
    # Time to 50% carbon loss (Half-life t0.5) in years
    t0.5_years = log(2) / k_real,
    t0.5_lower = log(2) / k_real_upper, # Higher k = shorter decay time
    t0.5_upper = log(2) / k_real_lower,
    
    # Time to 95% carbon loss (t0.05) in years
    t0.05_years = -log(0.05) / k_real,
    t0.05_lower = -log(0.05) / k_real_upper,
    t0.05_upper = -log(0.05) / k_real_lower
  ) %>%
  # Round numeric values for cleaner display
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  # Sort sites by decomposition speed (from fastest to slowest)
  arrange(t0.5_years)

# 2. Print table to console
print(k_table)

# 3. Export table to CSV for manuscript / Excel
write.csv(k_table, "fagus_dc1_decay_by_site2.csv", row.names = FALSE)



# 1. Extraction of parameters from summary/fixef
# (Alternatively, you can manually set: k_int <- -2.31; zi_int <- -3.50)
fix_eff <- fixef(fit_fagus_fixed)

k_int  <- fix_eff["k_Intercept", "Estimate"]
zi_int <- fix_eff["zi_Intercept", "Estimate"]



# OPTION 1: PURE WOOD DECAY RATE (Calculated solely from k_Intercept)
# Recommended for standard ecological publication
# ==============================================================================

# Transform log-k to real decomposition rate parameter
k_real <- exp(k_int)

# Half-life (50% remaining carbon)
t0.5_pure <- log(2) / k_real

# Time to 95% carbon loss (5% remaining carbon)
t0.05_pure <- -log(0.05) / k_real   # or log(20) / k_real



# OPTION 2: COMBINED POPULATION DECAY RATE (Including initial zero-inflation zi)
# Account for initial zero probability at t = 0
# ==============================================================================

# Convert zi_Intercept from logit scale to initial zero-inflation probability
p_zi_init <- 1 / (1 + exp(-zi_int))   # or plogis(zi_int)

# Half-life including initial zi (50% remaining population carbon)
t0.5_comb <- log(2 * (1 - p_zi_init)) / k_real

# Time to 95% carbon loss including initial zi
t0.05_comb <- log(20 * (1 - p_zi_init)) / k_real


# ==============================================================================
# SECOND ROUND BAYESIAN MODEL - FALLEN TREES RECORDED IN ALL DC
# ==============================================================================



library(ggplot2)
library(dplyr)
library(mgcv)

# List of sites to exclude due to low log counts
sites_to_remove <- c(
  "LFRLP-FAWF__203__Hollaenderschlag__NA", 
  "LFRLP-FAWF__232__Wuesttal__NA", 
  "LFRLP-FAWF__240__Schwappelbruch__NA", 
  "LFRLP-FAWF__246__Adelsberg__1",
  "URK__5__Zarnowka - Babia Gora N.P.__NA",
  "VUK__22__Velka Ples__NA", 
  "VUK__4__Bila Opava__NA"
)


# DATA PREPARATION & EXCLUSIONS
# ==============================================================================

trees_fagus <- trees_bayes_ready %>%
  filter(species_group == "Fagus") %>%
  # Filter out low-sample sites
  filter(!composed_site_id %in% sites_to_remove) %>%
  filter(
    !is.na(carbon_rel),
    !is.na(decaying_time),
    !is.na(diam_1_init),
    !is.na(composed_site_id),
    !is.na(decay_init)
  ) %>%
  mutate(
    # Create new site variable with merged Žofín sub-sites
    site_merged = case_when(
      composed_site_id %in% c("VUK__1__Zofin__a", "VUK__1__Zofin__c") ~ "VUK__1__Zofin__ac",
      composed_site_id %in% c("VUK__1__Zofin__b", "VUK__1__Zofin__d") ~ "VUK__1__Zofin__bd",
      TRUE ~ composed_site_id  # Keep all other sites unchanged
    ),
    
    decay_init = as.factor(decay_init),
    carbon_rel = ifelse(carbon_rel >= 1, 0.9999, carbon_rel),
    diam_z = as.numeric(scale(diam_1_init))
  )


# ==============================================================================
# OVERALL TRAJECTORY PLOT WITH GROUP-SPECIFIC & OVERALL GAMM CURVES
# ==============================================================================

ggplot(trees_fagus, aes(x = decaying_time, y = carbon_rel)) +
  # Jittered scatter points colored by initial decay stage
  geom_jitter(
    aes(color = decay_init, shape = decay_init), 
    alpha = 0.25, 
    width = 0.25, 
    height = 0.01, 
    size = 1.4
  ) +
  
  # 1. GAM smooth curves ONLY for decay_init 1, 3, and 4
  geom_smooth(
    data = filter(trees_fagus, decay_init %in% c("1", "3", "4")),
    aes(color = decay_init, group = decay_init),
    method = "gam",
    formula = y ~ s(x, bs = "cs", k = 5),
    se = FALSE, 
    linewidth = 0.9,
    linetype = "solid"
  ) +
  
  # 2. OVERALL GAM smooth curve for ALL trees combined (thick black line + green ribbon)
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    color = "#000000",        # Výrazná černá čára
    fill = "#a8dda8",         # Světle zelený pás spolehlivosti
    alpha = 0.45,
    linewidth = 1.3,
    se = TRUE
  ) +
  
  # Formatting scales and Dark2 palette
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.2)) +
  scale_color_brewer(palette = "Dark2", name = "Initial Decay Stage") +
  scale_shape_manual(values = c(16, 17, 15, 18, 8), name = "Initial Decay Stage") +
  
  # Labels and publication theme
  labs(
    title = "Relative Carbon Loss Trajectory for Beech (Fagus sylvatica)",
    subtitle = "Thick Black Line = Overall GAM trend; Colored Lines = GAM trends by Initial Decay Stage (decay_init)",
    x = "Decaying Time [Years] (decaying_time)",
    y = "Relative Carbon Fraction (carbon_rel)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

library(tidybayes)

# 1. Prepare grid for prediction
grid_data <- data.frame(
  decaying_time = seq(0, 50, length.out = 100),
  diam_z = 0
)

# 2. Extract ONLY the mu component (ignoring zero-inflation zi effect on average line)
preds_mu <- grid_data %>%
  add_epred_draws(fit_fagus_fixed, dpar = "mu", re_formula = NA, allow_new_levels = TRUE, ndraws = 200)

# 3. Plot trajectory starting near 1.0
ggplot() +
  geom_point(
    data = trees_fagus_dc1, 
    aes(x = decaying_time, y = carbon_rel), 
    alpha = 0.15, color = "darkgreen"
  ) +
  stat_lineribbon(
    data = preds_mu, 
    aes(x = decaying_time, y = .epred), 
    .width = c(.95, .80), 
    fill = "#a8dda8", 
    color = "#1b6738"
  ) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Beech Wood Density/Carbon Loss (Mu Parameter Trajectory)",
    subtitle = "Excludes zi-inflation pull; curve correctly anchors near t=0, C_rel=1.0",
    x = "Decaying Time [Years]",
    y = "Relative Carbon Fraction (mu)"
  ) +
  theme_bw()


# FACETED PLOT BY SITE (INSPECTING SITE VARIABILITY)
# ==============================================================================

# Overview of decomposition trajectories broken down by individual sites
ggplot(trees_fagus, aes(x = decaying_time, y = carbon_rel)) +
  geom_jitter(
    aes(color = decay_init), 
    alpha = 0.4, 
    width = 0.2, 
    height = 0.01, 
    size = 1.2
  ) +
  
  # Local GAM smooth line per site
  geom_smooth(
    method = "gam", 
    formula = y ~ s(x, bs = "cs", k = 4), # Adjusted k for smaller site sample sizes
    se = FALSE, 
    color = "#1b6738", 
    linewidth = 0.8
  ) +
  
  # Separate panel for each site
  facet_wrap(~ site_merged, scales = "free_x") + 
  
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.2)) +
  scale_color_brewer(palette = "Dark2", name = "Initial Decay Stage") +
  labs(
    title = "Beech Wood Carbon Decomposition by Site (composed_site_id)",
    subtitle = "Highlights variability in decomposition rates across climate regions",
    x = "Decaying Time [Years]",
    y = "Relative Carbon Fraction"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Abies alba - fragment, not ready yet

# B) Subset for Abies alba (Fir) starting in Decay Class 1
trees_abies_dc3 <- trees_bayes_ready %>%
  filter(
    species_group == "Abies alba",
    decay_init == 3
  )

# Non-linear Weibull model to account for the lag phase and site variability (Abies)
priors_weibull <- c(
  prior(normal(-2.5, 1), nlpar = "k"),
  prior(gamma(2, 1), nlpar = "alpha", lb = 0) # Prior forcing alpha > 0
)

fit_abies_weibull <- brm(
  formula = bf(
    carbon_rel ~ exp(- (exp(k) * decaying_time)^alpha),
    k ~ 1 + diam_1_init + (1 | composed_site_id) + (1 | tree_unique_id),
    alpha ~ 1,
    nl = TRUE
  ),
  data = trees_abies_class1,
  family = zero_inflated_beta(),
  prior = priors_weibull,
  chains = 4, cores = 4, iter = 3000
)