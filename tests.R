# ==============================================================================
# CONSISTENCY TESTS - FILTERING CONSISTENT_ID YES
# vyřešit stromy, které začínají consistent_id NULL a tím skončí
# ==============================================================================


library(dplyr)
library(tidyr)

# --- DISTRIBUTION BY INVENTORY RANK ---
# Checking if rank_inventory == 1 is always NA and how TRUE/FALSE develops
consistency_by_rank <- trees_3inv %>%
  group_by(rank_inventory) %>%
  summarise(
    count_NA = sum(is.na(consistent_id)),
    count_TRUE = sum(consistent_id == TRUE, na.rm = TRUE),
    count_FALSE = sum(consistent_id == FALSE, na.rm = TRUE),
    total_records = n(),
    .groups = "drop"
  ) %>%
  mutate(
    pct_NA = round(100 * count_NA / total_records, 1),
    pct_TRUE = round(100 * count_TRUE / total_records, 1),
    pct_FALSE = round(100 * count_FALSE / total_records, 1)
  )

cat("\n--- DISTRIBUTION OF CONSISTENT_ID BY INVENTORY RANK ---\n")
print(consistency_by_rank)

# --- CLASSIFYING ALL TREES BY CONSISTENCY HISTORY ---
# We analyze the entire timeline for every tree_unique_id in trees_3inv
tree_consistency_report <- trees_3inv %>%
  group_by(tree_unique_id) %>%
  # Chronological sort within group
  arrange(rank_inventory, .by_group = TRUE) %>%
  summarise(
    # Get the sequence of consistent_id (skipping the first NA)
    timeline = paste(na.omit(consistent_id), collapse = " -> "),
    # Count unique states (excluding NA)
    unique_states = n_distinct(consistent_id, na.rm = TRUE),
    # Check if FALSE ever appeared
    ever_false = any(consistent_id == FALSE, na.rm = TRUE),
    # Check if TRUE ever appeared
    ever_true = any(consistent_id == TRUE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    category = case_when(
      timeline == ""               ~ "Only one inventory (No consistency check)",
      unique_states > 1            ~ "Switchers (TRUE <-> FALSE)",
      ever_false & !ever_true      ~ "Always FALSE",
      ever_true & !ever_false      ~ "Always TRUE",
      TRUE                         ~ "Other / Error"
    )
  )

# --- SUMMARY STATISTICS ---
# Calculating total counts and percentages for the entire dataset
final_summary <- tree_consistency_report %>%
  group_by(category) %>%
  summarise(
    number_of_trees = n(),
    .groups = "drop"
  ) %>%
  mutate(
    percentage = round(100 * number_of_trees / sum(number_of_trees), 2)
  ) %>%
  arrange(desc(number_of_trees))

cat("\n--- FINAL CONSISTENCY REPORT (ALL TREES IN trees_3inv) ---\n")
print(final_summary)

# --- DETAILED LOOK AT TRANSITIONS ---
# For those that are not simple "Always TRUE" or "Always FALSE"
cat("\n--- DETAILED TRANSITION PATTERNS ---")
transition_detail <- tree_consistency_report %>%
  filter(category %in% c("Switchers (TRUE <-> FALSE)", "Other / Error")) %>%
  count(timeline, sort = TRUE)

print(transition_detail)

# --- IDENTIFY SWITCHER IDs ---
# We use the previous logic to find trees that flip states
switcher_data <- trees_3inv %>%
  group_by(tree_unique_id) %>%
  filter(n_distinct(consistent_id, na.rm = TRUE) > 1) %>%
  # Keep only one row per tree to get the plot info
  slice(1) %>%
  ungroup()

# --- CREATE THE SUMMARY TABLE ---
# Listing the specific locations where these inconsistent trees are found
switchers_location_list <- switcher_data %>%
  select(institute, site_id, site_name, plot_id, tree_unique_id) %>%
  group_by(institute, site_id, site_name, plot_id) %>%
  summarise(
    number_of_switcher_trees = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(number_of_switcher_trees))

# --- DISPLAY THE TABLE ---
cat("\n--- PLOTS CONTAINING INCONSISTENT 'SWITCHER' TREES ---\n")
if (nrow(switchers_location_list) > 0) {
  print(as.data.frame(switchers_location_list))
} else {
  cat("No switcher trees found in the dataset.\n")
}

# Optional: Export this list to CSV if you need to investigate in Excel
write.csv(switchers_location_list, "switcher_plots_report.csv", row.names = FALSE)



# --- 1. IDENTIFY RECRUIT CONSISTENCY ERRORS ---
# A tree has an error if consistent_id is NA even after its first appearance
# It reveals also duplicated tree_ids
# ------------------------------------------------------------------------------

trees_consistent <- trees_consistent %>%
  group_by(tree_unique_id) %>%
  arrange(rank_inventory, .by_group = TRUE) %>%
  mutate(
    # Check if NA exists in 2nd, 3rd... occurrence of the tree
    # (The first occurrence of a tree always has a 'tree-level' rank of 1)
    tree_rank = row_number(),
    has_bad_na = any(tree_rank > 1 & is.na(consistent_id)),
    cons_recruits_error = if_else(has_bad_na, 1, 0)
  ) %>%
  select(-tree_rank, -has_bad_na) %>% # Cleanup temp columns
  ungroup()

# --- 2. GLOBAL SUMMARY OF ERRORS ---
error_total_summary <- trees_consistent %>%
  summarise(
    total_trees = n_distinct(tree_unique_id),
    trees_with_error = n_distinct(tree_unique_id[cons_recruits_error == 1]),
    pct_error = round(100 * trees_with_error / total_trees, 2)
  )

cat("\n--- GLOBAL RECRUIT ERROR SUMMARY ---\n")
print(error_total_summary)

# --- 3. SUMMARY BY INSTITUTE AND SITE ---
error_site_summary <- trees_consistent %>%
  filter(cons_recruits_error == 1) %>%
  group_by(institute, site_id, site_name) %>%
  summarise(
    n_error_trees = n_distinct(tree_unique_id),
    .groups = "drop"
  ) %>%
  arrange(desc(n_error_trees))

cat("\n--- RECRUIT ERRORS BY LOCATION ---\n")
if (nrow(error_site_summary) > 0) {
  print(as.data.frame(error_site_summary))
} else {
  cat("No recruit consistency errors found. Database logic is correct.\n")
}


# DUPLICATE CHECK ----
# Hlavně duplicity count 2, jestli jich nejsou stovky na ploše (nebo Magda)
# ------------------------------------------------------------------------------

# --- CHECKING FOR DUPLICATE RECORDS WITHIN THE SAME YEAR ---
# Identifying if any tree_unique_id appears more than once in a single inventory
# ------------------------------------------------------------------------------

# 1. Identify combinations of ID and Year that appear multiple times
yearly_duplicates <- trees %>%
  group_by(tree_unique_id, inventory_year) %>%
  summarise(duplicity_count = n(), .groups = 'drop') %>%
  filter(duplicity_count > 1)

# 2. Join this info back to the main table and clean up the count column
trees <- trees %>%
  left_join(
    yearly_duplicates, 
    by = c("tree_unique_id", "inventory_year")
  ) %>%
  # If no duplicate was found, duplicity_count is 1 (the record itself)
  mutate(duplicity_count = if_else(is.na(duplicity_count), 1, as.numeric(duplicity_count)))

# --- REPORTING ---
cat("\n--- DUPLICITY REPORT ---\n")
if (nrow(yearly_duplicates) > 0) {
  total_dup_rows <- sum(yearly_duplicates$duplicity_count)
  unique_dup_trees <- n_distinct(yearly_duplicates$tree_unique_id)
  
  cat("WARNING: Duplicates found!\n")
  cat("  - Affected ID/Year combinations: ", nrow(yearly_duplicates), "\n")
  cat("  - Total rows involved:           ", total_dup_rows, "\n")
  cat("  - Unique trees affected:         ", unique_dup_trees, "\n")
} else {
  cat("SUCCESS: Each stem has exactly one record per inventory year.\n")
}

# --- OPTIONAL: INVESTIGATE DUPLICATES BY LOCATION ---
if (nrow(yearly_duplicates) > 0) {
  dup_location_summary <- trees %>%
    filter(duplicity_count > 1) %>%
    group_by(institute, site_name) %>%
    summarise(n_duplicates = n(), .groups = "drop") %>%
    arrange(desc(n_duplicates))
  
  cat("\n--- DUPLICATES BY LOCATION ---\n")
  print(as.data.frame(dup_location_summary))
}
