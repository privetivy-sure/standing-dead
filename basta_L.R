# ==============================================================================
# BaSTA round 1
# Different adjustments dealing with BaSTA
# ==============================================================================

# Currently applies for BaSTA 1.9.4.

write.table(BaSTA_L_quercus, 
            file = "BaSTA_L_quercus.txt", 
            sep = "\t",          
            dec = ".",           
            row.names = FALSE, 
            quote = FALSE,       
            na = "",             
            fileEncoding = "UTF-8")


# FILTER ERROR-FREE and no dead RECORDS
# ==============================================================================

trees_clear_L <- trees_consistent_L %>%
  filter(
    errorlife_tree == FALSE & 
      discard_nodead == FALSE 
  )

# CREATE TEMPORAL TYPE CLASSIFICATION BASED ON BIRTH AND DEATH RECORDS
# ==============================================================================

trees_clear_L <- trees_clear_L %>%
  mutate(
    temporal_type = case_when(
      birth > 0 & death > 0  ~ "observed",
      birth == 0 & death > 0 ~ "left-truncated",
      birth > 0 & death == 0 ~ "right-censored",
      birth == 0 & death == 0 ~ "unobserved",
      TRUE                    ~ NA_character_ # Fallback for any unexpected values (e.g., negative numbers or NAs)
    )
  )



# REPAIRING ERRORS IN THE FINAL DATASET
# ==============================================================================
# --- IDENTIFY TREES WITH NA's IN LIFE, dbh_first, birth, death ---
# ------------------------------------------------------------------------------

library(dplyr)
library(readr)

# 1. IDENTIFY ALL UNIQUE TREE IDs WITH ANY NA VALUES IN THE TARGET COLUMNS
ids_with_nas <- trees_clear_L %>%
  # Filter rows where at least one critical variable is missing
  filter(is.na(status3) | is.na(diam_first) | is.na(birth) | is.na(death)) %>%
  distinct(tree_unique_id) %>%
  pull(tree_unique_id)

# 2. EXTRACT THE COMPLETE HISTORY FOR THESE TREES (FOR QUALITY CHECK)
trees_with_nas_history <- trees_clear_L %>%
  filter(tree_unique_id %in% ids_with_nas)

# 3. EXPORT THE EXTRACTED HISTORY TO A CSV FILE
write_csv(trees_with_nas_history, "trees_with_critical_nas_history.csv")

# 4. REMOVE THESE TREES FROM THE MAIN DATAFRAME
trees_clear_L <- trees_clear_L %>%
  filter(!tree_unique_id %in% ids_with_nas)




# (PREPARE covMat FOR BaSTA ANALYSIS) - OLD
# ==============================================================================
# We extract a unique list of trees with their corresponding categorical covariate.
# The first column is strictly named 'ID'.

# ------------------------------------------------------------------------------
# CATEGORICAL COVARITE: FOREST TYPE
# ------------------------------------------------------------------------------
# Aggregate the data to have exactly one row per unique tree
covMat <- trees_surv_50 %>%
  group_by(tree_unique_id) %>%
  # Extract the first available record of fortype, dbh_first and tree_species for each tree
  summarize(
    fortype   = first(fortype),
    dbh_first = first(dbh_first),
    genus = first(genus),
    .groups   = "drop" # Automatically ungroup after summarizing
  ) %>%
  # Rename the unique ID column to 'ID' as strictly required by BaSTA
  rename(ID = tree_unique_id) %>%
  rename(type_ = fortype) %>%
  rename(eea_ = EEA_fortype)
#CRITICAL CONDITION: Keep ONLY rows belonging to selected genus
  filter(genus == "Fagus")

# Generate the covariate matrix. 
covMatSharp <- MakeCovMat(x = ~ genus_ + dbh_first, data = covMat)

# Convert the resulting matrix to a dataframe so it can be safely joined via dplyr
covMatSharp_df <- as.data.frame(covMatSharp)

# MODIFY DUMMY COVARIATE COLUMNS DIRECTLY IN covMatSharp_df
# ------------------------------------------------------------------------------
# 1. Remove the column "type_Alpine.coniferous.forest" (excluding the category)

if ("type_Alpine.coniferous.forest" %in% colnames(covMatSharp_df)) {
  covMatSharp_df <- covMatSharp_df %>% select(-type_Alpine.coniferous.forest)
}

if ("type_other" %in% colnames(covMatSharp_df)) {
  covMatSharp_df <- covMatSharp_df %>% select(-type_other)
}

if ("type_Oak.beech" %in% colnames(covMatSharp_df)) {
  covMatSharp_df <- covMatSharp_df %>% select(-type_Oak.beech)
}

# 2. Rename the remaining forest type columns as requested
covMatSharp_df <- covMatSharp_df %>%
  rename(
    type_Beech             = any_of("type_Beech.forest"),
    type_Mountainous.Beech = any_of("type_Mountainous.beech.forest"),
    # type_Oak.beech         = any_of("type_Oak.mixed.forests"),
    # type_Alpine            = any_of("type_Alpine.coniferous.forest")
  )



# ==============================================================================
# covMat - CATEGORICAL COVARITE: GENUS / GENUS GROUPS
# ==============================================================================
library(dplyr)


# ==============================================================================
# COVARIATE DATA PREPARATION FOR ONE SPECIES
# ==============================================================================
# Filter dataset for Abies alba, extract first diameter per tree stem,
# calculate z-score (diam_std) and filter missing values.

covMat_quercus <- trees_clear_L %>%
  filter(species_group == "Quercus") %>%
  group_by(tree_unique_id) %>%
  summarize(
    diam_first = first(diam_first),
    .groups    = "drop"
  ) %>%
  rename(ID = tree_unique_id) %>%
  filter(!is.na(diam_first)) %>%
  mutate(
    # Standardize diam_first (Mean = 0, SD = 1)
    diam_std = as.vector(scale(diam_first))
  )


# Z-SCORE STANDARDIZED DIAMETER MATRIX FOR SINGLE SPECIES
# ==============================================================================
# 1. Create covariate matrix using only standardized diameter (no species factor)
covMatSharp_std <- MakeCovMat(x = ~ diam_std, data = covMat_quercus)

# 2. Convert to data frame for BaSTA input
covMatSharp_std_df <- as.data.frame(covMatSharp_std)

# Optional: Store mean and SD for future predictions/back-transformation
abies_dbh_mean <- mean(covMat_abies$diam_first)
abies_dbh_sd   <- sd(covMat_abies$diam_first)


# ACCOUNTING FOR TEMPERATURE AND RAINFALL
# ==============================================================================
# Compute z-score, mean and SD from all unique stems per dataset (one species)
# e.g. the same way as for dbh - each stem has its own temp. and rainf.




# ==============================================================================
# PREPARE birthDeath MATRIX FOR BaSTA ANALYSIS
# ==============================================================================
# We extract a unique list of trees with their respective birth and death years.
# Columns are named strictly: 'ID', 'BIRTH', and 'DEATH'.

birthDeath <- trees_clear_L %>%
  group_by(tree_unique_id) %>%
  summarize(
    # Extract birth and death years for each unique tree
    BIRTH = first(birth),
    DEATH = first(death),
    .groups = "drop"
  ) %>%
  # Rename the tree identifier column to 'ID' as requested by BaSTA
  rename(ID = tree_unique_id)



# ==============================================================================
# PREPARE censusMat FOR BaSTA ANALYSIS
# ==============================================================================
library(dplyr)
library(tidyr)

# --- 1. DEFINE FULL STUDY TIME SPAN ---
all_years <- min(trees_clear_L$inventory_year, na.rm = TRUE):max(trees_clear_L$inventory_year, na.rm = TRUE)

# --- 2. EXTRACT OBSERVED MORTALITY RECORDS ---
census_dead <- trees_clear_L %>%
  filter(life == "D") %>%
  select(tree_unique_id, inventory_year) %>%
  mutate(recorded_dead = 1) %>%
  distinct()

# --- 3. ALIGN LONG-FORMAT DATA TO FULL TIME SERIES ---
census_complete <- census_dead %>%
  complete(tree_unique_id, inventory_year = all_years, fill = list(recorded_dead = 0))

# --- 4. RESHAPE TO CENSUS CAPTURE-HISTORY MATRIX ---
censusMat_1 <- census_complete %>%
  pivot_wider(
    names_from = inventory_year, 
    values_from = recorded_dead,
    names_sort = TRUE
  ) %>%
  rename(ID = tree_unique_id)

# --- 5. ENFORCE CHRONOLOGICAL MATRIX STRUCTURE ---
censusMat <- censusMat_1[, c("ID", as.character(all_years))]


# ==============================================================================
# COMPOSING BaSTA inputMat
# ==============================================================================
library(dplyr)

# Ensure all ID columns are of the same data type before joining
birthDeath$ID   <- as.character(birthDeath$ID)
censusMat$ID  <- as.character(censusMat$ID)
covMatSharp_std_df$ID   <- as.character(covMatSharp_std_df$ID)

# The inner_joins will now automatically filter birthDeath and censusMat 
# to include ONLY the IDs that remain in the selected genus filter in covMat
BaSTA_L_quercus <- birthDeath %>%
  inner_join(censusMat, by = "ID") %>%
  inner_join(covMatSharp_std_df, by = "ID")


str(BaSTA_L_fraxinus)



# ------------------------------------------------------------------------------
# --- ROUND BIRTH and DEATH ---
# ------------------------------------------------------------------------------
# Neccassary for dataCheck
BaSTA_L_quercus <- BaSTA_L_quercus %>%
  mutate(
    BIRTH = ifelse(BIRTH > 0, round(BIRTH), BIRTH),
    DEATH = ifelse(DEATH > 0, round(DEATH), DEATH)
  )

# Search for NA censuses
# ------------------------------------------------------------------------------
library(tidyr)
na_report <- BaSTA_genus_L_zs %>%
  select(ID, as.character(1972:2025)) %>% 
  pivot_longer(
    cols = -ID, 
    names_to = "inventory_year", 
    values_to = "value"
  ) %>%
  filter(is.na(value)) 

# Display
View(na_report)


# ==============================================================================
# DATA CHECK
# ==============================================================================


verify <- DataCheck(
  object = BaSTA_L_quercus, 
  studyStart = 1972, 
  studyEnd = 2025, 
  autofix = rep(1, 7), 
  silent = FALSE
)



# ALTERNATIVE ERROR DETECTION, REMOVAL OF WRONG ROWS
# ------------------------------------------------------------------------------
# They are removed now - repare them later.


# 1. Row indices reported directly by DataCheck()
problem_rows <- c(
  88
)

# 2. Extract IDs of these trees from the BaSTA matrix
error_ids <- BaSTA_L_quercus$ID[problem_rows]

# 3. View these 21 problematic trees in your ORIGINAL source dataset
trees_errors_source <- trees_clear_L %>%
  filter(tree_unique_id %in% error_ids)

View(trees_errors_source)


# 2. Remove flawed rows from BaSTA
BaSTA_L_quercus <- BaSTA_L_quercus[-problem_rows, ]


# BaSTA for testing
# ------------------------------------------------------------------------------

genus_L <- basta(
  object = BaSTA_genus_L_zs_clean, 
  studyStart = 1972, 
  studyEnd = 2025, 
  covarsStruct = "fused", 
  model = "GO", 
  shape = "simple",
  nsim = 1,              
  parallel = FALSE,
  niter = 1000,          
  burnin = 100,
  thinning = 5,
  updateJumps = FALSE    # Turn off the MCMC chains
)
summary(genus_L)

# FULL BASTA MODEL
# ==============================================================================


# Multi BaSTA
# ------------------------------------------------------------------------------
quercus_multi_fused <- multibasta(
  object = BaSTA_L_quercus, 
  studyStart = 1972, 
  studyEnd = 2025,
  models = c("GO", "WE", "LO"),            # Gompertz, Weibull, Logistic
  shapes = c("simple", "Makeham"),        
  covarsStruct = "fused",
  nsim = 3,
  ncpus = 3,
  parallel = TRUE,
  niter = 10000,                           
  burnin = 2000, 
  thinning = 10
)

summary(quercus_multi_fused)

# Regular BaSTA
# ------------------------------------------------------------------------------
quercus_WE_fused <- basta(
  object = BaSTA_L_quercus, 
  studyStart = 1972, 
  studyEnd = 2025, 
  covarsStruct = "fused", 
  model = "WE", 
  shape = "simple",
  nsim = 3,
  ncpus = 3,
  parallel = TRUE,
  niter = 50000,       
  burnin = 10000,       
  thinning = 50,     
  updateJumps = TRUE
)


# Zobrazení porovnání podle DIC
print(abies_multi$DIC)

summary(ulmus_multi_fused)

# Complete BASTA plot function:
BaSTA:::plot.basta

# Save RDS
# -------------------------------------------------------------------------------
saveRDS(quercus_WE_fused, file = "L_quercus_WE_fused.rds")
acer_GO_fused <- readRDS(file = "L_acer_GO_fused.rds")

plot(quercus_WE_fused,plot.trace=FALSE,xlim=c(0, 100), noCI = FALSE)

# ORDER OF CATEGORICAL VARIABLES (alphabetical)
names(genus_200_799_full$survQuant)
summary(all_full)



# =========================================================================
# Interpolation function to calculate exact age (X-axis) for given S(x) (Y-axis)
# =========================================================================
library(purrr)
library(dplyr)

# Extract survival time quantiles (mean and 95% C.I.) for specified target survival levels S(x)
get_survival_quantiles_ci <- function(basta_obj, target_levels = c(0.5, 0.05)) {
  
  # Extract list of survival quantile matrices from the BaSTA model output
  surv_list <- basta_obj$survQuant
  
  # Clean group labels by removing unwanted prefixes
  group_names <- sub("^[_x.]*(genus_group)?", "", names(surv_list))
  
  # Loop over each covariate group matrix
  results <- map2_df(surv_list, group_names, function(surv_mat, g_name) {
    # X-axis: Age/time vector extracted from column names
    ages <- as.numeric(colnames(surv_mat))
    
    # Y-axis vectors from matrix rows: 1 = Mean, 3 = Lower 95% C.I. (2.5%), 4 = Upper 95% C.I. (97.5%)
    s_mean  <- surv_mat[1, ]
    s_lower <- surv_mat[3, ]
    s_upper <- surv_mat[4, ]
    
    # Perform linear interpolation via approx() for each target survival threshold S(x)
    res_list <- map(target_levels, function(lvl) {
      age_mean  <- approx(x = s_mean,  y = ages, xout = lvl)$y
      age_lower <- approx(x = s_lower, y = ages, xout = lvl)$y
      age_upper <- approx(x = s_upper, y = ages, xout = lvl)$y
      
      out <- c(
        mean  = round(age_mean, 2),
        lower = round(age_lower, 2),
        upper = round(age_upper, 2)
      )
      
      # Name columns according to the specified S(x) level
      names(out) <- paste0("age_S", lvl, "_", c("mean", "lower_ci", "upper_ci"))
      return(out)
    })
    
    # Flatten list into a single vector row and attach the group name
    combined_res <- unlist(res_list)
    c(group = g_name, combined_res)
  })
  
  # Convert text-formatted numeric columns back to double/numeric
  results <- results %>%
    mutate(across(-group, as.numeric))
  
  return(results)
}

# Run extraction on BaSTA model object
surv_summary <- get_survival_quantiles_ci(fraxinus_GO, target_levels = c(0.5, 0.05))
print(surv_summary)



# =========================================================================
# Calculation function for user-input DBH including 95% C.I.
# =========================================================================

# Calculate persistence times (t0.50 and t0.05) with 95% C.I. for any input diameter (mm)
calculate_persistence_fraxinus_ci <- function(dbh_mm) {
  # 1. Baseline persistence values (mean and 95% C.I. in years) for an average trunk (474 mm)
  t05_mean  <- 23.4; t05_lower  <- 21.2; t05_upper  <- 25.8
  t005_mean <- 39.5; t005_lower <- 35.8; t005_upper <- 43.9
  
  # 2. Population diameter parameters for tree species
  mean_d <- 474.0
  sd_d   <- 286.0
  
  # 3. Diameter effect coefficient from BaSTA (prop.haz)
  gamma  <- -0.22927
  
  # 4. Standardize the user-input diameter (z-score)
  dbh_std <- (dbh_mm - mean_d) / sd_d
  
  # 5. Recalculate mean and 95% C.I. persistence times based on proportional hazards
  exp_factor <- exp(-gamma * dbh_std)
  
  t_050_mean  <- t05_mean  * exp_factor
  t_050_lower <- t05_lower * exp_factor
  t_050_upper <- t05_upper * exp_factor
  
  t_005_mean  <- t005_mean  * exp_factor
  t_005_lower <- t005_lower * exp_factor
  t_005_upper <- t005_upper * exp_factor
  
  return(data.frame(
    DBH_mm            = dbh_mm,
    t_0.50_mean       = round(t_050_mean, 1),
    t_0.50_CI_lower   = round(t_050_lower, 1),
    t_0.50_CI_upper   = round(t_050_upper, 1),
    t_0.05_mean       = round(t_005_mean, 1),
    t_0.05_CI_lower   = round(t_005_lower, 1),
    t_0.05_CI_upper   = round(t_005_upper, 1)
  ))
}

# Test the function with example diameters (200 mm, average 474 mm, 800 mm)
calculate_persistence_fraxinus_ci(dbh_mm = c(200, 474, 800))
# ---------------------------------------------------------

