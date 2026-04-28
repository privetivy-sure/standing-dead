# ==============================================================================
# DATA PREPARATION: UNIQUE IDs AND INVENTORY FILTERING
# Project: Standing Deadwood Analysis (Cox Model Preparation)
# ==============================================================================

library(tidyr)
library(dplyr)
library(stringr)
library(stringi)
library(DBI)

trees <- dbGetQuery(con, "
  SELECT * FROM trees 
  WHERE position = 'S'") 

plots <- dbGetQuery(con, "
  SELECT * FROM plots
")

library(dplyr)

# --- JOIN TABLES ---
# Join selected columns from 'plots' into 'trees' based on 'plot_record_id'
# left_join ensures we keep all records in the 'trees' table
trees <- trees %>%
  left_join(
    plots %>% select(plot_record_id, plot_id, inventory_type, plot_sampled),
    by = "plot_record_id"
  )

# --- JOIN VALIDATION ---
# Check if any trees failed to find a match in the plots table
na_count <- sum(is.na(trees$plot_id))

if (na_count > 0) {
  warning(paste("⚠️ Data inconsistency:", na_count, "rows in 'trees' have no matching 'plot_record_id' in 'plots'."))
} else {
  message("✅ Join successful: All trees matched with plot data.")
}

# --- GENERATE PLOT UNIQUE IDs ---
# Extracting hierarchical levels and keeping them as separate columns
# ------------------------------------------------------------------------------

trees <- trees %>%
  mutate(
    # A. Split composed_site_id by double underscores
    comp_split = stri_split_fixed(composed_site_id, "__", simplify = TRUE),
    
    # B. Create standalone columns (source: composed_site_id)
    institute = if_else(comp_split[, 1] == "" | is.na(comp_split[, 1]), "NA", as.character(comp_split[, 1])),
    site_id   = if_else(comp_split[, 2] == "" | is.na(comp_split[, 2]), "NA", as.character(comp_split[, 2])),
    site_name = if_else(comp_split[, 3] == "" | is.na(comp_split[, 3]), "NA", as.character(comp_split[, 3])),
    
    # C. Handle sub_id with specific logic (NA string is common here)
    sub_id = case_when(
      is.na(comp_split[, 4]) | comp_split[, 4] == "" | comp_split[, 4] == "NA" | comp_split[, 4] == "\\N" ~ "NA",
      TRUE ~ as.character(comp_split[, 4])
    ),
    
    # D. Extract dataset from inventory_id (2nd part)
    dataset = stri_split_fixed(inventory_id, "__", simplify = TRUE)[, 2],
    dataset = if_else(
      is.na(dataset) | dataset == "" | dataset == "\\N" | dataset == "NA", 
      "NA", as.character(dataset)
    ),
    
    # E. Create the final plot-level unique identifier
    # Using the standardized columns we just created above
    plot_unique_id = paste(
      institute, 
      site_id, 
      sub_id, 
      site_name, 
      plot_id,          # Joined from 'plots' table in step 2
      inventory_type,   # Joined from 'plots' table in step 2
      dataset, 
      sep = "_"
    )
  ) %>%
  # Remove only the temporary matrix and the raw dataset helper
  select(-comp_split)

# --- CREATE TREE_UNIQUE_ID ----
# Handling stem_id errors and generating unique stem-level identifiers
# ------------------------------------------------------------------------------

# Search for problematic stem_ids (NA, 0, or empty strings)
stem_errors <- trees %>%
  filter(is.na(stem_id) | stem_id %in% c("NA", "", 0, "0"))

if (nrow(stem_errors) > 0) {
  cat("⚠️ Found", nrow(stem_errors), "suspicious stem_id records. These will be standardized to 'NA'.\n")
}

# Generate the unique stem identifier
trees <- trees %>%
  mutate(
    # Ensure tree_id and stem_id are treated as character to avoid paste errors
    stem_id_fix = if_else(is.na(stem_id) | stem_id == "" | stem_id == "0", "NA", as.character(stem_id)),
    tree_id_char = as.character(tree_id),
    
    # Combine plot base with tree and stem identifiers
    tree_unique_id = paste(plot_unique_id, tree_id_char, stem_id_fix, sep = "_")
  ) %>%
  # Remove temporary fix columns
  select(-stem_id_fix, -tree_id_char) %>%
  # Final sorting by ID and Year to ensure chronological history per stem
  arrange(tree_unique_id, inventory_year)

# Final Verification
cat("✅ TREE_UNIQUE_ID generation complete.\n")
cat("Example ID:", trees$tree_unique_id[1], "\n")
cat("Total rows in 'trees' table:", nrow(trees), "\n")


# Check for any potential NAs in the final IDs
total_nas <- sum(is.na(trees$tree_unique_id))
cat("\nNumber of NA values in tree_unique_id:", total_nas, "\n")


# ==============================================================================
# INVENTORY SUMMARY & FILTERING (>= 3 INVENTORIES) ----
# Identifying plots with enough temporal depth for survival analysis
# ==============================================================================

# Calculate number of inventories directly in the trees table
trees <- trees %>%
  group_by(plot_unique_id) %>%
  mutate(
    no_inventories = n_distinct(inventory_year),
    is_3inv = if_else(no_inventories >= 3, "Y", "N")
  ) %>%
  ungroup()

# --- Create the summary table for reporting ---
# Including plot metadata for easier filtering and identification
summary_inventory <- trees %>%
  group_by(plot_unique_id) %>%
  summarise(
    # Core inventory metrics
    no_inventories = first(no_inventories),
    is_3inv        = first(is_3inv),
    years_sampled  = paste(sort(unique(inventory_year)), collapse = ", "),
    
    # Metadata columns (kept for easier reporting)
    institute      = first(institute),
    site_id        = first(site_id),
    site_name      = first(site_name),
    sub_id         = first(sub_id),
    dataset        = first(dataset),
    
    .groups = "drop"
  )

# Quick console check
cat("\n--- INVENTORY DEPTH SUMMARY ---\n")
plots_3plus <- sum(summary_inventory$no_inventories >= 3)
cat("Number of plots with >= 3 inventories:", plots_3plus, "\n")
print(head(summary_inventory))


# FILTERING TREES FOR ANALYSIS ----
# ------------------------------------------------------------------------------

# Create the final filtered table for analysis
trees_3inv <- trees %>%
  filter(is_3inv == "Y") %>%
  # Sort by plot, tree, stem, and year
  # Using the newly created tree_unique_id for perfectly organized records
  arrange(tree_unique_id, inventory_year)

cat("\n--- DATA VOLUME CHECK ---\n")
cat("Total tree records in original table: ", nrow(trees), "\n")
cat("Tree records in filtered table (3inv):", nrow(trees_3inv), "\n")
cat("Percentage of data retained:          ", 
    round((nrow(trees_3inv) / nrow(trees)) * 100, 2), "%\n")

# Verification of the first few IDs in the final table
cat("\nPreview of filtered table structure:\n")
trees_3inv %>% 
  select(tree_unique_id, inventory_year, no_inventories) %>% 
  head(10) %>% 
  print()


## --- INVENTORY RANKING ----
# Assigning chronological rank to each inventory per plot directly in the main table
# ------------------------------------------------------------------------------

trees_3inv <- trees_3inv %>%
  group_by(plot_unique_id) %>%
  # dense_rank sorts years chronologically and assigns 1, 2, 3...
  mutate(rank_inventory = dense_rank(inventory_year)) %>%
  ungroup() %>%
  # Keep the data organized
  arrange(tree_unique_id, inventory_year)

# --- QUICK VERIFICATION ---
cat("\n--- RANKING CHECK ---\n")
trees_3inv %>% 
  select(plot_unique_id, inventory_year, rank_inventory) %>% 
  distinct() %>% 
  head(10) %>% 
  print()

# ==============================================================================
# CONSISTENT_ID CORRECTION
# ==============================================================================
# --- ADDING CONSISTENT_SWITCH COLUMN ---
# 1 = Tree flipped between TRUE and FALSE
# 0 = Tree stayed stable (Only TRUE, only FALSE, or only NA)
# Switchers can appear - tree was re-identified in second inventory, but not
# re-identified in third inventory. etc.
# This column is not important for further filtering
# It tells only if consistent_id was changed, but it includes FALSE throughout
# ------------------------------------------------------------------------------

trees_3inv <- trees_3inv %>%
  group_by(tree_unique_id) %>%
  mutate(
    # n_distinct(..., na.rm = TRUE) counts how many different non-NA values exist
    consistent_switch = if_else(n_distinct(consistent_id, na.rm = TRUE) > 1, 1, 0)
  ) %>%
  ungroup()

# --- VERIFICATION ---
# Check the distribution of the new column
switch_summary <- trees_3inv %>%
  group_by(consistent_switch) %>%
  summarise(
    total_records = n(),
    unique_trees = n_distinct(tree_unique_id),
    .groups = "drop"
  )

cat("\n--- CONSISTENT_SWITCH DISTRIBUTION ---\n")
print(switch_summary)

# --- SHOWCASE EXAMPLES ---
# Display a few trees where the switch occurred
if (any(trees_3inv$consistent_switch == 1)) {
  cat("\n--- EXAMPLES OF SWITCHER TREES (consistent_switch == 1) ---\n")
  trees_3inv %>%
    filter(consistent_switch == 1) %>%
    select(tree_unique_id, rank_inventory, inventory_year, consistent_id, consistent_switch) %>%
    arrange(tree_unique_id, rank_inventory) %>%
    head(15) %>%
    print()
}

# HARMONIZATION OF CONSISTENT_ID (overwriting NA)
library(dplyr)
library(tidyr)

# --- CREATE 'consistent_tree' COLUMN (Tree history logic) ---
# Pravidlo: Pokud má strom v kterémkoliv roce consistent_id == FALSE, 
# celá jeho historie (všechny řádky daného tree_unique_id) je označena jako FALSE.
trees_3inv <- trees_3inv %>%
  group_by(tree_unique_id) %>%
  mutate(
    consistent_tree = if_else(any(consistent_id == FALSE, na.rm = TRUE), FALSE, TRUE)
  ) %>%
  ungroup()


# --- FINAL FILTRATION - CREATING trees_consistent TABLE ---
# Do finální tabulky jdou pouze stromy, které jsou konzistentní napříč všemi roky.
trees_consistent <- trees_3inv %>%
  filter(consistent_tree == TRUE)


# --- CONTROL OUTPUTS ---
# 1. Row-level calculations
total_rows <- nrow(trees_3inv)
final_rows <- nrow(trees_consistent)
excluded_rows <- total_rows - final_rows
pct_rows_kept <- round((final_rows / total_rows) * 100, 2)

# 2. Tree-level calculations (Unique tree_unique_id)
total_trees <- n_distinct(trees_3inv$tree_unique_id)
final_trees <- n_distinct(trees_consistent$tree_unique_id)
excluded_trees <- total_trees - final_trees
pct_trees_kept <- round((final_trees / total_trees) * 100, 2)

# 3. Printing the report
cat("\n--- DATA CONSISTENCY REPORT ---\n")
cat("RECORD LEVEL (Rows):\n")
cat("  Original row count:     ", total_rows, "\n")
cat("  Final row count:        ", final_rows, "\n")
cat("  Excluded rows:          ", excluded_rows, " (", 100 - pct_rows_kept, "%)\n")
cat("  Data kept:              ", pct_rows_kept, "%\n\n")

cat("INDIVIDUAL STEM LEVEL (Unique IDs):\n")
cat("  Original unique trees:  ", total_trees, "\n")
cat("  Final unique trees:     ", final_trees, "\n")
cat("  Excluded trees:         ", excluded_trees, " (", 100 - pct_trees_kept, "%)\n")
cat("  Stems kept:             ", pct_trees_kept, "%\n")
cat("-------------------------------\n")



# ==============================================================================
# FILTERING TREES WITH AT LEAST ONE "DEAD" RECORD
# ==============================================================================
library(dplyr)

# 1. Filter trees that have at least one "D" (dead) record in their history
# We use the property where first_death_rank is Inf for trees that stayed alive
trees_dead <- trees_consistent %>%
  filter(first_death_rank != Inf) %>%
  # Sort by ID and time for better overview
  arrange(tree_unique_id, rank_inventory)

# 2. Control output
cat("Original row count:", nrow(trees_consistent), "\n")
cat("New row count (history of dead trees):", nrow(trees_dead), "\n")
cat("Number of unique dead trees:", n_distinct(trees_dead$tree_unique_id), "\n")



# LIFE-CYCLE LOGIC & MORTALITY RANKING
# Analyzing tree history to identify death events and "zombie" errors
# ==============================================================================

trees_consistent <- trees_consistent %>%
  group_by(tree_unique_id) %>%
  # Crucial: Ensure chronological order within each tree's history
  arrange(rank_inventory, .by_group = TRUE) %>%
  mutate(
    # 1. Identify the EARLIEST rank where the tree was marked as Dead (D)
    # If the tree never died, min() returns Inf
    first_death_rank = suppressWarnings(min(rank_inventory[life == "D"], na.rm = TRUE)),
    
    # 2. Flag 'was_dead_in_past': 
    # TRUE if current rank is greater than the rank of first recorded death
    was_dead_in_past = if_else(!is.infinite(first_death_rank) & 
                                 rank_inventory > first_death_rank, TRUE, FALSE),
    
    # 3. Detect 'error_life' (Zombie Error):
    # Flag as 1 if the tree was dead in the past but is recorded as Alive (A) now
    error_life = if_else(was_dead_in_past & life == "A", 1, 0),
    
    # 4. Global life status: 
    # Mark as "D" if the tree recorded at least one death event in its history
    life_total = if_else(!is.infinite(first_death_rank), "D", "A")
  ) %>%
  # 5. Mark the entire tree history as invalid if a zombie error occurred anywhere
  mutate(has_lifecycle_error = any(error_life == 1)) %>%
  ungroup()

# --- CLEANING & REPORTING ---

zombie_counts <- trees_consistent %>% 
  filter(has_lifecycle_error) %>% 
  summarise(n = n_distinct(tree_unique_id))

cat("Number of trees excluded due to life-cycle errors (zombies):", zombie_counts$n, "\n")

# Final removal of problematic trees
trees_consistent <- trees_consistent %>%
  filter(has_lifecycle_error == FALSE) %>%
  select(-has_lifecycle_error)





# --- CREATE 'status' ATTRIBUTE - Stump classification by height ---
# ==============================================================================

trees_dead <- trees_dead %>%
  mutate(
    # Temporary handling of NA values for string concatenation
    L = if_else(is.na(life), "", as.character(life)),
    I = if_else(is.na(integrity), "", as.character(integrity)),
    
    # Basic status (combining two letters, e.g., "DC", "DF")
    status = paste0(L, I),
    
    # Special rule for DP (Dead Stump)
    # Applied if life is "D", integrity is "F" (or NA), 
    # and height is between 0.1 and 1.2 m
    status = if_else(
      life == "D" & 
        (integrity == "F" | is.na(integrity)) & 
        !is.na(height) & height >= 0.1 & height <= 1.2,
      "DP",
      status
    ),
    
    # If the status is an empty string (both were NA), return NA
    status = if_else(status == "", NA_character_, status)
  ) %>%
  # Remove temporary helper columns
  select(-L, -I)

# --- STATUS SUMMARY REPORT ---

status_summary <- trees_dead %>%
  count(status) %>%
  mutate(pct = round(n / sum(n) * 100, 2))

print("--- Overview of created statuses ---")
print(status_summary)


# ---- EXTENDED ATTRIBUTES EXTRACTION ----
# ==============================================================================

library(dplyr)
library(tidyr)
library(stringr)


# EXTRACTION OF STATUS2 FROM EXTENDED_ATTRIBUTES
# ------------------------------------------------------------------------------

# 1. Extract status2 and create a summary report of combinations
institute_status2_report <- trees_dead %>%
  # Extract status2 from extended_attributes using regex
  mutate(
    status2_extracted = str_extract(extended_attributes, "status2:\\s*([^,}]+)") %>%
      str_remove("status2:\\s*") %>%
      trimws()
  ) %>%
  # Group by institute and the newly extracted status2
  group_by(institute, status2_extracted) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  # Sort to keep combinations from the same institute together
  arrange(institute, desc(count))

# 2. Rename column for cleaner output
colnames(institute_status2_report)[2] <- "status2"

# 3. Export the resulting summary table
# Note: Using a relative path is better for GitHub portability
write.table(institute_status2_report, 
            file = "institute_status2_summary.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = TRUE)

# 4. Print report to console
print(institute_status2_report)


# APPLY EXTRACTION TO THE MAIN DATASET
# ------------------------------------------------------------------------------
# We use the regex: "status2:\\s*([^,} ]+)"
# It finds "status2:", skips spaces, and captures everything until the next comma, 
# space, or closing brace.

trees_dead <- trees_dead %>%
  mutate(
    status2 = str_match(extended_attributes, "status2:\\s*([^,} ]+)")[, 2]
  )


# JOIN STATUS DESCRIPTIONS ----
# Joining the status table which contains descriptions of status2 values from metadata
# ------------------------------------------------------------------------------

trees_dead <- trees_dead %>%
  left_join(
    status %>% select(institute, status2, value_description), 
    by = c("institute" = "institute", "status2" = "status2")
  )


# 1. PREVIEW UNIQUE COMBINATIONS ----
# Checking if descriptions were correctly mapped
check_join <- trees_dead %>%
  filter(!is.na(status2)) %>%
  distinct(institute, status2, value_description) %>%
  arrange(institute, status2)

print("--- Overview of mapped codes and their descriptions ---")
print(check_join)


# 2. IDENTIFY MISSING DESCRIPTIONS ----
# If there are codes in the data that are missing from the lookup table, list them:
# (Note: Assumes match_stats or similar check was performed)

unmatched_list <- trees_dead %>%
  filter(!is.na(status2) & is.na(value_description)) %>%
  distinct(institute, status2)

if (nrow(unmatched_list) > 0) {
  print("--- Warning: The following codes were not found in the status lookup table (by institute) ---")
  print(unmatched_list)
} else {
  message("Success: All status2 codes were successfully matched with descriptions.")
}