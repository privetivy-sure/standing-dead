# ==============================================================================
# BaSTA round 1
# Different adjustments dealing with problematic version of BaSTA
# ==============================================================================

# Currently applies for BaSTA 1.9.4.

write.table(trees_surv, 
            file = "trees_surv.txt", 
            sep = "\t",          
            dec = ".",           
            row.names = FALSE, 
            quote = FALSE,       
            na = "",             
            fileEncoding = "UTF-8")



# ==============================================================================
# REPAIRING ERRORS IN THE FINAL DATASET
# ==============================================================================
# --- IDENTIFY TREES WITH NA's IN LIFE, dbh_first, birth, death ---
# ------------------------------------------------------------------------------

library(dplyr)
library(readr)

# 1. IDENTIFY ALL UNIQUE TREE IDs WITH ANY NA VALUES IN THE TARGET COLUMNS
ids_with_nas <- trees_surv_50 %>%
  # Filter rows where at least one critical variable is missing
  filter(is.na(life) | is.na(dbh_first) | is.na(birth) | is.na(death)) %>%
  distinct(tree_unique_id) %>%
  pull(tree_unique_id)

# 2. EXTRACT THE COMPLETE HISTORY FOR THESE TREES (FOR QUALITY CHECK)
trees_with_nas_history <- trees_surv_50 %>%
  filter(tree_unique_id %in% ids_with_nas)

# 3. EXPORT THE EXTRACTED HISTORY TO A CSV FILE
write_csv(trees_with_nas_history, "trees_with_critical_nas_history.csv")

# 4. REMOVE THESE TREES FROM THE MAIN DATAFRAME
trees_dead9 <- trees_dead9 %>%
  filter(!tree_unique_id %in% ids_with_nas)




# ==============================================================================
# CREATE TEMPORAL TYPE CLASSIFICATION BASED ON BIRTH AND DEATH RECORDS
# ==============================================================================

trees_surv <- trees_surv %>%
  mutate(
    temporal_type = case_when(
      birth > 0 & death > 0  ~ "observed",
      birth == 0 & death > 0 ~ "left-truncated",
      birth > 0 & death == 0 ~ "right-censored",
      birth == 0 & death == 0 ~ "unobserved",
      TRUE                    ~ NA_character_ # Fallback for any unexpected values (e.g., negative numbers or NAs)
    )
  )


# STEP 3: FILTER ERROR-FREE RECORDS
# ==============================================================================
# We subset the 'Fagus' dataset to keep only rows where both tree-level 
# quality control flags (errorlife_tree and errorstatus_tree) are equal to 0.

trees_surv_clean <- trees_surv_50 %>%
  filter(errorlife_tree == 0 & errorstatus_tree == 0)




# STEP 7: PREPARE covMat FOR BaSTA ANALYSIS
# ==============================================================================
# We extract a unique list of trees with their corresponding categorical covariate.
# The first column is strictly named 'ID'.

# Aggregate the data to have exactly one row per unique tree
covMat <- trees_surv_50 %>%
  group_by(tree_unique_id) %>%
  # Extract the first available record of fortype, dbh_first and tree_species for each tree
  summarize(
    fortype   = first(fortype),
    dbh_first = first(dbh_first),
    species_group = first(species_group),
    .groups   = "drop" # Automatically ungroup after summarizing
  ) %>%
  # Rename the unique ID column to 'ID' as strictly required by BaSTA
  rename(ID = tree_unique_id) %>%
  rename(type_ = fortype) %>%
  # CRITICAL CONDITION: Keep ONLY rows belonging to "FagusGroup"
  filter(genus == "Picea")

# Generate the covariate matrix. 
covMatSharp <- MakeCovMat(x = ~ type_ + dbh_first, data = covMat)

# Convert the resulting matrix to a dataframe so it can be safely joined via dplyr
covMatSharp_df <- as.data.frame(covMatSharp)

# MODIFY DUMMY COVARIATE COLUMNS DIRECTLY IN covMatSharp_df
# ------------------------------------------------------------------------------
# 1. Remove the column "type_Alpine.coniferous.forest" (excluding the category)
#if ("type_Alpine.coniferous.forest" %in% colnames(covMatSharp_df)) {
#  covMatSharp_df <- covMatSharp_df %>% select(-type_Alpine.coniferous.forest)
#}

#if ("type_other" %in% colnames(covMatSharp_df)) {
#  covMatSharp_df <- covMatSharp_df %>% select(-type_other)
#}

# 2. Rename the remaining forest type columns as requested
covMatSharp_df <- covMatSharp_df %>%
  rename(
    type_Beech             = any_of("type_Beech.forest"),
    type_Mountainous.Beech = any_of("type_Mountainous.beech.forest"),
    type_Oak.beech         = any_of("type_Oak.mixed.forests"),
    type_Alpine            = any_of("type_Alpine.coniferous.forest")
  )

# STEP 8: PREPARE birthDeath MATRIX FOR BaSTA ANALYSIS
# ==============================================================================
# We extract a unique list of trees with their respective birth and death years.
# Columns are named strictly: 'ID', 'BIRTH', and 'DEATH'.

birthDeath <- trees_surv_50 %>%
  group_by(tree_unique_id) %>%
  summarize(
    # Extract birth and death years for each unique tree
    BIRTH = first(birth),
    DEATH = first(death),
    .groups = "drop"
  ) %>%
  # Rename the tree identifier column to 'ID' as requested by BaSTA
  rename(ID = tree_unique_id)


# ------------------------------------------------------------------------------
# --- ROUND BIRTH and DEATH ---
# ------------------------------------------------------------------------------
# Important to round, otherwise DataCheck within BaSTA doesn't compute well
birthDeath <- birthDeath %>%
  mutate(
    BIRTH = ifelse(BIRTH > 0, round(BIRTH), BIRTH),
    DEATH = ifelse(DEATH > 0, round(DEATH), DEATH)
  ) 


# ==============================================================================
# STEP 9: PREPARE censusMat FOR BaSTA ANALYSIS
# ==============================================================================
library(dplyr)
library(tidyr)

# --- 1. DEFINE FULL STUDY TIME SPAN ---
all_years <- min(trees_surv_50$inventory_year, na.rm = TRUE):max(trees_surv_50$inventory_year, na.rm = TRUE)

# --- 2. EXTRACT OBSERVED MORTALITY RECORDS ---
census_dead <- trees_surv_50 %>%
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
covMatSharp_df$ID   <- as.character(covMatSharp_df$ID)

# The inner_joins will now automatically filter birthDeath and censusMat 
# to include ONLY the IDs that remain in the selected genus filter in covMat
BaSTA_fagus <- birthDeath %>%
  inner_join(censusMat, by = "ID") %>%
  inner_join(covMatSharp_df, by = "ID")


str(BaSTA_fagus)


# Search for NA censuses
# ------------------------------------------------------------------------------
library(tidyr)
na_report <- BaSTA_fagus %>%
  select(ID, as.character(1936:2026)) %>% 
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


verify <- DataCheck(BaSTA_fagus, studyStart = 1936, studyEnd = 2026, autofix = rep(1, 7), silent = FALSE)



# ALTERNATIVE ERROR DETECTION
# ------------------------------------------------------------------------------

# 1. Identify recapture year columns (numeric columns)
all_cols <- colnames(BaSTA_fagus)
year_cols <- all_cols[!all_cols %in% c("ID", "BIRTH", "DEATH", "fortype", "dbh_first", "weight")]

# 2. Find row indices where recapture is 1 in the birth year
true_error_indices <- which(sapply(1:nrow(BaSTA_fagus), function(i) {
  b_year <- as.character(BaSTA_fagus$BIRTH[i])
  if (b_year %in% year_cols) {
    return(BaSTA_fagus[i, b_year] == 1)
  }
  return(FALSE)
}))

# 3. Extract the ACTUAL wrong rows using the correct indices
if (length(true_error_indices) > 0) {
  basta_TRUE_errors <- BaSTA_fagus[true_error_indices, ]
  
  cat("\n--- TRUE ERROR DETECTION COMPLETED ---\n")
  cat("Found", length(true_error_indices), "actual conflicting rows in BaSTA_fagus.\n")
  cat("----------------------------------------\n")
  
  # View the REAL problematic trees
  View(basta_TRUE_errors)
  
  # Create a truly clean matrix
  BaSTA_fagus_clean <- BaSTA_fagus[-true_error_indices, ]
} else {
  cat("\nNo actual birth-recapture conflicts found directly in BaSTA_fagus!\n")
}


# SIMPLE BaSTA FOR TESTING
# ================================================================================================

fagus_test <- basta(
  object = BaSTA_fagus, 
  studyStart = 1936, 
  studyEnd = 2026, 
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
# ORDER OF CATEGORICAL VARIABLES (alphabetical)
names(picea_test2$survQuant)
summary(fagus_test)


# FULL BASTA MODEL - still takes too long to compute
# ==============================================================================
SM <- basta(
  object = BaSTA_input, 
  studyStart = 1972, 
  studyEnd = 2024, 
  covarsStruct = "fused", 
  model = "GO", 
  shape = "simple",
  nsim = 5,
  ncpus = 5,
  parallel = TRUE,
  niter = 15000,       
  burnin = 2000,       
  thinning = 20,     
  updateJumps = TRUE
)

# Complete BASTA plot function:
BaSTA:::plot.basta

# Save RDS
# -------------------------------------------------------------------------------
saveRDS(full_test, file = "full_test.rds")
SM <- readRDS(file = "basta_sm1.rds")

plot(full_test,plot.trace=FALSE,xlim=c(0, 70))

# IMPROVED GRAPH
# ===============================================================================

my_basta_plot(full_test, 
              plot.trace = FALSE, 
              xlim = c(0, 70),
              names.legend = c("", "", ""))

par(mfg = c(1, 1))

legend("bottomleft", 
       legend = c(
         "Beech Mountainous",    # eea_bmtn (red)
         "Beech Submountainous", # eea_bsubmtn (blue)
         "Spruce Mountainous"    # eea_obm (green)
       ), 
       col = c("#E41A1C", "#377EB8", "#4DAF4A"), # Základní paleta BaSTA barev
       lty = 1,
       lwd = 5,
       bty = "n",
       cex = 1.1)


