write.table(plot_timeline, 
            file = "plot_timeline.txt", 
            sep = "\t",          
            dec = ".",           
            row.names = FALSE, 
            quote = FALSE,       
            na = "",             
            fileEncoding = "UTF-8")



# ==============================================================================
# DATA PREPARATION: UNIQUE IDs AND INVENTORY FILTERING
# Project: Standing Deadwood Analysis (Cox Model Preparation)
# ==============================================================================

library(tidyr)
library(dplyr)
library(stringr)
library(stringi)
library(DBI)

con <- dbConnect(
  RPostgres::Postgres(),
  host = Sys.getenv("DB_HOST"),
  port = Sys.getenv("DB_PORT"),
  dbname = Sys.getenv("DB_NAME"),
  user = Sys.getenv("DB_USER"),
  password = Sys.getenv("DB_PASS")
)
dbListTables(con)

library(RPostgres)


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
# CONSISTENT_ID CORRECTION and FILTRATION
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
# Rule: If a tree has consistent_id == FALSE in any year,
# its entire history (all rows of the given tree_unique_id) is marked as FALSE.
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

# Export trees_consistent
write.table(trees_consistent, 
            file = "trees_consistent.txt", 
            sep = "\t",           
            row.names = FALSE,     
            quote = FALSE,         
            fileEncoding = "UTF-8") 




library(dplyr)

# EXPORTS AND IMPORTS
# ==============================================================================

  # Export trees_dead
  write.table(trees_dead, 
              file = "trees_dead.txt", 
              sep = "\t",           
              row.names = FALSE,     
              quote = FALSE,         
              fileEncoding = "UTF-8") 

# Import trees_dead
library(readr)

file_path <- "D:/DEADWOOD_D/WildCard/Residence/standing-dead-git/trees_dead.txt"

trees_dead <- read_delim(file_path, 
                               delim = "\t", 
                               escape_double = FALSE, 
                               trim_ws = TRUE,
                               locale = locale(encoding = "UTF-8"))

# Loading the metadata lookup table
status_meta <- read_delim("status.txt", delim = "\t", locale = locale(encoding = "UTF-8"))  

# ==============================================================================
# FILTER TREES WITH AT LEAST ONE "DEAD" RECORD
# ==============================================================================
# DETECT BIOLOGICAL ERRORS (ZOMBIES)
# ------------------------------------------------------------------------------
library(dplyr)

trees_dead <- trees_consistent %>%
  group_by(tree_unique_id) %>%
  filter(any(life == "D")) %>%
  # Crucial: Check if 'rank_inventory' is the correct column name!
  arrange(rank_inventory, .by_group = TRUE) %>%
  mutate(
    # Find the first record of death
    # We use min() because we are already grouped by tree_unique_id
    first_death_rank = min(rank_inventory[life == "D"]),
    
    # Identify zombies (Alive after being Dead)
    errorlife_record = if_else(life == "A" & rank_inventory > first_death_rank, 1, 0),
    
    # Mark the whole history of that tree if it ever failed
    errorlife_tree = any(errorlife_record == 1)
  ) %>%
  ungroup()

# --- REPORTING ---

# Count unique trees that became "zombies"
zombie_trees_count <- trees_dead %>% 
  filter(errorlife_tree) %>% 
  summarise(n = n_distinct(tree_unique_id)) %>% 
  pull(n)

cat("Original row count:", nrow(trees_consistent), "\n")
cat("Dead trees history row count:", nrow(trees_dead), "\n")
cat("Number of unique dead trees:", n_distinct(trees_dead$tree_unique_id), "\n")
cat("Number of trees with biological errors (D to A):", zombie_trees_count, "\n")


# ==============================================================================
# --- CREATE 'status' ATTRIBUTE ---
# Improved logic to handle NULL integrity and identify DP (Stumps)
# ==============================================================================
library(stringr)

trees_dead <- trees_dead %>%
  mutate(
    status = case_when(
      # 1. Special rule for DP (Dead Stump): 
      # Dead, (Fragmented OR NULL integrity), and height <= 1.2m
      life == "D" & 
        (integrity == "F" | is.na(integrity)) & 
        between(height, 0.1, 1.2) ~ "DP",
      
      # 2. If integrity is NULL, status is just the life value (A or D)
      is.na(integrity) ~ as.character(life),
      
      # 3. Standard case: Combine life and integrity (e.g., AC, DC, DF)
      TRUE ~ str_c(life, integrity)
    )
  )

# --- STATUS SUMMARY REPORT ---

status_summary <- trees_dead %>%
  group_by(status) %>%
  summarise(
    count = n(),
    avg_height = round(mean(height, na.rm = TRUE), 2)
  ) %>%
  arrange(desc(count))

print("--- Overview of created statuses ---")
print(status_summary)



library(dplyr)
library(tidyr)
library(stringr)


# ==============================================================================
# EXTRACTION OF STATUS2 FROM EXTENDED_ATTRIBUTES
# ==============================================================================
trees_dead <- trees_dead %>%
# Extract 'status2' value from the 'extended_attributes' JSON-like string
  mutate(
    # Extracting status2 value into the new column
    status2 = str_match(extended_attributes, "status2[\"']?:\\s*[\"']?([^\"',} ]+)")[, 2],
    
    # Cleaning the extracted strings
    status2 = str_remove_all(status2, "[\"']"),
    status2 = trimws(status2)
  )

# LIST OF VALUES ACCORDING TO INSTITUTE
# ------------------------------------------------------------------------------
# This creates the list you wanted to review before harmonization
institute_status2_summary <- trees_dead %>%
  filter(!is.na(status2)) %>%
  group_by(institute, status2) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  # Join descriptions to see what each institute means by their code
  left_join(
    status_meta %>% select(institute, status2, value_description) %>% distinct(),
    by = c("institute", "status2")
  ) %>%
  arrange(institute, status2)

# 4. EXPORT AND CONTROL
# -------------------------------------------------------------------------------
# Export of the summary for your review
write.table(institute_status2_summary, "status2_raw_summary.txt", sep = "\t", row.names = FALSE)


library(dplyr)

library(dplyr)

# HARMONIZATION AND UNIFICATION OF STATUS2 VALUES
# ------------------------------------------------------------------------------
# This script maps various institute-specific codes to a standardized 
# "status" format (combination of life and integrity) with a few exceptions.

trees_dead <- trees_dead %>%
  mutate(
    # Create a unified status column based on extracted codes and their descriptions
    status2_unified = case_when(
      # --- ALIVE (A) ---
      # Mapping to AC (Alive Complete)
      status2 %in% c("ASI", "AI", "LSv", "LS(v)") ~ "AC",  
      
      # Mapping to AF (Alive Fragmented)
      status2 %in% c("ASB", "AB", "LSs", "DAB", "DALB") ~ "AF", 
      
      # Mapping to AFC (Alive Fragmented Cracked)
      status2 %in% c("ASC", "AC") ~ "AFC",                
      
      # Mapping to AL (Alive Lying/Uprooted)
      status2 %in% c("AU", "AL", "ALM") ~ "AL",           
      
      # --- DEAD (D) ---
      # Mapping to DC (Dead Complete)
      status2 %in% c("DSI", "DI", "Di", "TSv", "TS(v)") ~ "DC", 
      
      # Mapping to DF (Dead Fragmented)
      status2 %in% c("DSB", "DB", "TSs", "TSsk") ~ "DF",  
      
      # Mapping to DP (Dead Stump)
      status2 %in% c("DP") ~ "DP",                        
      
      # Mapping to DL (Dead Lying)
      status2 %in% c("DU") ~ "DL",                        
      
      # --- SPECIAL / OTHER CATEGORIES ---
      # Dead standing but unspecified integrity
      status2 %in% c("TS") ~ "D",                         
      
      # Mixed Dead/Alive status (specific to NWFVA)
      status2 %in% c("ASv", "ASs", "AS(v)") ~ "A/DC",     
      
      # Missing or Dismissed objects
      status2 %in% c("MIS", "MIS_ABS", "MIS_CWD") ~ "MIS", 
      
      # Dismissed due to ultimate decay
      status2 %in% c("MIS_DEC") ~ "MIS_DEC",              
      
      # Keep original value if the code is not recognized in the list above
      TRUE ~ status2 
    )
  )

# VALIDATION OF RESULTS
# ------------------------------------------------------------------------------

# Create a summary table to check the success of the unification process
unification_check <- trees_dead %>%
  group_by(institute, status2, status2_unified) %>%
  summarise(
    count = n(), 
    .groups = "drop"
  ) %>%
  arrange(institute, status2_unified)

# Print the validation table to console
print("--- Summary of status2 unification by institute ---")
print(unification_check)

# EXPORT THE UNIFIED SUMMARY
# -----------------------------------------------------------------------------
# Exporting the summary for manual verification of the mapping logic
write.table(unification_check, 
            file = "status2_unified_check.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

library(dplyr)
library(readr)
library(stringr)


# CHECK OF DUPLICATES (FOR COMPOSED_SITE_ID and INVENTORY_YEAR)
# ==============================================================================

library(dplyr)

# Spuštění testu na duplicity
duplicity_test <- trees_dead %>%
  # Filtrujeme pouze na zájmové lokality, kde k problému dochází
  filter(composed_site_id %in% c(
    "VUK__1__Zofin__a", "VUK__1__Zofin__b", 
    "VUK__1__Zofin__c", "VUK__1__Zofin__d", 
    "VUK__15__Ranspurk__NA"
  )) %>%
  # Seskupíme podle roku, lokality a ID stromu/kmene
  # (Pozn.: Pokud se sloupec pro ID stromu jmenuje jinak než tree_id, upravte název)
  group_by(composed_site_id, inventory_year, tree_id, stem_id) %>%
  # Spočítáme, kolikrát tam tato kombinace je
  summarise(pocet_zaznamu = n(), .groups = "drop") %>%
  # Zajímají nás pouze ty, které jsou tam 2x a více
  filter(pocet_zaznamu > 1)

# Výpis výsledku do konzole
print(paste("Počet nalezených duplicitních kmenů:", nrow(duplicity_test)))
head(duplicity_test)


# ==============================================================================
# REMOVING TREES MEASURED TWICE IN ONE YEAR
# ==============================================================================

library(dplyr)

# 1. CATEGORIZE BY DBH AND CLEAN OVERLAPS *ONLY* WITHIN TARGET SITES
# ==============================================================================

trees_dead2 <- trees_dead %>%
  # Create DBH groups for all trees
  mutate(
    dbh_group = case_when(
      dbh >= 10 & dbh <= 100 ~ "10-100 mm",
      dbh > 100              ~ "> 100 mm",
      TRUE                   ~ "Undefined / NA"
    )
  ) %>%
  # Create a temporary grouping key: 
  # For target sites, it creates a unique biological key per year.
  # For ALL OTHER sites, it creates a unique key per row (using row_number()),
  # which guarantees that slice(1) will never delete them.
  group_by(
    group_key = if_else(
      composed_site_id %in% c("VUK__1__Zofin__a", "VUK__1__Zofin__b", "VUK__1__Zofin__c", "VUK__1__Zofin__d", "VUK__15__Ranspurk__NA"),
      paste(composed_site_id, inventory_year, tree_id, stem_id, sep = "_"),
      paste0("keep_all_", row_number())
    )
  ) %>%
  
  # Sort by plot_id to keep Plot 1 as priority for the target sites
  arrange(plot_id, .by_group = TRUE) %>%
  
  # This now removes duplicates ONLY in Žofín and Ranšpurk.
  # Other sites have a unique group_key per row, so slice(1) keeps them all.
  slice(1) %>%
  ungroup() %>%
  
  # Remove the temporary helper column
  select(-group_key)
  


# 2. GENERATING DATA CLEANING AND VERIFICATION REPORTS
# ==============================================================================

# Identify and count discarded duplicate records (rows from Plot 2 that were removed)
report_discarded <- trees_dead %>%
  mutate(
    dbh_group = case_when(
      dbh >= 10 & dbh <= 100 ~ "10-100 mm",
      dbh > 100              ~ "> 100 mm",
      TRUE                   ~ "Undefined / NA"
    )
  ) %>%
  # Find which rows from the original data are missing in trees_dead2
  filter(!paste(composed_site_id, inventory_year, tree_id, stem_id, plot_id, sep = "_") %in% 
           paste(trees_dead2$composed_site_id, trees_dead2$inventory_year, trees_dead2$tree_id, trees_dead2$stem_id, trees_dead2$plot_id, sep = "_")) %>%
  group_by(composed_site_id, inventory_year, dbh_group) %>%
  summarise(discarded_duplicates = n(), .groups = "drop")


# Count preserved trees in the new dataset for the target sites only (for verification)
report_preserved_target_sites <- trees_dead2 %>%
  filter(composed_site_id %in% c("VUK__1__Zofin__a", "VUK__1__Zofin__b", "VUK__1__Zofin__c", "VUK__1__Zofin__d", "VUK__15__Ranspurk__NA")) %>%
  group_by(composed_site_id, inventory_year, dbh_group) %>%
  summarise(preserved_trees = n(), .groups = "drop")


# 3. CONSOLE OUTPUT REPORT
# ==============================================================================

cat("VERIFICATION REPORT: PLOT OVERLAP CLEANING VIA BIOLOGICAL KEY\n")
cat("======================================================================\n")
cat("Total rows in ORIGINAL table (ALL sites): ", nrow(trees_dead), "\n")
cat("Total rows in NEW table (ALL sites):      ", nrow(trees_dead2), "\n")
cat("Total duplicate rows removed:             ", nrow(trees_dead) - nrow(trees_dead2), " rows\n\n")

cat("----------------------------------------------------------------------\n")
cat("1. DISCARDED DUPLICATES REPORT (Should ONLY show target sites):\n")
cat("----------------------------------------------------------------------\n")
if(nrow(report_discarded) == 0) {
  cat("No duplicate rows were removed.\n")
} else {
  print(report_discarded)
}

cat("2. TARGET SITES PRESERVED TREES (Verification that small trees are safe):\n")
cat("----------------------------------------------------------------------\n")
print(report_preserved_target_sites)


library(dplyr)


# ==============================================================================
# ==============================================================================
# CREATING BIRTH AND DEATH
# ==============================================================================
# ==============================================================================
#
library(dplyr)

target_sites <- c("VUK__1__Zofin__a", "VUK__1__Zofin__b", "VUK__1__Zofin__c", "VUK__1__Zofin__d", "VUK__15__Ranspurk__NA")

trees_dead2 <- trees_dead2 %>%
  mutate(plot_id2 = if_else(composed_site_id %in% target_sites, "1", as.character(plot_id)))


# 1. TIMELINE WITH RE-ALIGNED PLOT BOUNDARIES FOR ZOFIN AND RANSPURK (USING plot_id2)
# ==============================================================================
# Use the updated timeline with merged plots for the 5 target VUK sites
plot_timeline_clean <- plot_timeline %>%
  mutate(plot_id2 = if_else(composed_site_id %in% target_sites, "1", as.character(plot_id))) %>%
  distinct(composed_site_id, plot_id2, inventory_year, rank_inventory) %>%
  group_by(composed_site_id, plot_id2) %>%
  arrange(rank_inventory, .by_group = TRUE) %>%
  mutate(
    previous_inv_year = lag(inventory_year),
    next_inv_year     = lead(inventory_year),
    max_plot_rank     = max(rank_inventory, na.rm = TRUE)
  ) %>%
  ungroup()


# 2. SURVIVAL METRICS AGGREGATION (STREAMLINED SITE-PLOT-TREE HIERARCHY)
# ==============================================================================
# Extract unified survival intervals and censoring parameters per stem
stem_survival_metrics <- trees_dead2 %>%
  left_join(
    plot_timeline_clean, 
    by = c("composed_site_id", "plot_id2", "inventory_year", "rank_inventory")
  ) %>%
  group_by(composed_site_id, plot_id2, tree_id, stem_id) %>%
  summarise(
    has_standing_dead   = any(status %in% c("DC", "DF", "D"), na.rm = TRUE),
    year_entry_dead     = if_else(has_standing_dead, min(inventory_year[status %in% c("DC", "DF", "D") & rank_inventory == first_death_rank], na.rm = TRUE), NA_real_),
    year_last_standing  = if_else(has_standing_dead, max(inventory_year[status %in% c("DC", "DF", "D")], na.rm = TRUE), NA_real_),
    rank_last_standing  = if_else(has_standing_dead, max(rank_inventory[status %in% c("DC", "DF", "D")], na.rm = TRUE), NA_real_),
    
    # Secure cohort parameters from full study timeline
    first_death_rank    = min(first_death_rank, na.rm = TRUE),
    year_alive_prior    = if_else(has_standing_dead, min(previous_inv_year[rank_inventory == first_death_rank], na.rm = TRUE), NA_real_),
    next_inv_post_dead  = if_else(has_standing_dead, min(next_inv_year[rank_inventory == max(rank_inventory[status %in% c("DC", "DF", "D")], na.rm = TRUE)], na.rm = TRUE), NA_real_),
    
    first_obs_is_stump  = any(status == "DP" & rank_inventory == first_death_rank, na.rm = TRUE),
    has_turned_stump    = any(status == "DP" & inventory_year > year_last_standing, na.rm = TRUE),
    year_stump_observed = if_else(has_turned_stump, min(inventory_year[status == "DP" & inventory_year > year_last_standing], na.rm = TRUE), NA_real_),
    max_plot_rank       = max_plot_rank[1],
    
    # Exclude rows with data validation errors from risk set
    is_invalid_history  = any(errorlife_tree == TRUE | errorlife_record == 1 | is.na(life) | is.na(errorlife_tree), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(year_entry_dead, year_last_standing, year_stump_observed, year_alive_prior, next_inv_post_dead), 
                ~ if_else(is.infinite(.x), NA_real_, .x)))


# 3. ESTIMATION OF INITIALIZATION (BIRTH_2) AND CENSORING/EVENT TIMES (DEATH_2)
# ==============================================================================
# Pure mathematical midpoint estimations over the unified plot space
stem_numeric_intervals <- stem_survival_metrics %>%
  mutate(
    birth = case_when(
      is_invalid_history | first_obs_is_stump | !has_standing_dead ~ NA_real_,
      first_death_rank == 1 ~ 0, # Left-truncated entry
      !is.na(year_alive_prior) & !is.na(year_entry_dead) ~ (year_alive_prior + year_entry_dead) / 2, # Interval entry
      TRUE ~ NA_real_
    ),
    death = case_when(
      is_invalid_history | !has_standing_dead ~ NA_real_,
      rank_last_standing == max_plot_rank & !has_turned_stump ~ 0, # Right-censored at study termination
      has_turned_stump ~ (year_last_standing + year_stump_observed) / 2, # Event time: Transition to stump
      !is.na(year_last_standing) & !is.na(next_inv_post_dead) ~ (year_last_standing + next_inv_post_dead) / 2, # Event time: Missing row interpolation
      TRUE ~ NA_real_
    )
  )


# 4. FINAL ASSEMBLY AND TREE-LEVEL COHORT & CENSORING CLASSIFICATION
# ==============================================================================
trees_dead3 <- trees_dead2 %>%
  # Link estimated interval endpoints back to the master dataset
  left_join(
    stem_numeric_intervals %>% select(
      composed_site_id, plot_id2, tree_id, stem_id, birth, death, 
      is_invalid_history, has_standing_dead, first_obs_is_stump, has_turned_stump
    ), 
    by = c("composed_site_id", "plot_id2", "tree_id", "stem_id")
  ) %>%
  
  # Grouping strictly by tree history to guarantee uniform flags across records
  group_by(composed_site_id, plot_id2, tree_id, stem_id) %>%
  mutate(
    birth_null = if_else(all(is.na(birth)), 1, 0),
    death_null = if_else(all(is.na(death)), 1, 0),
    
    trajectory_note = case_when(
      is_invalid_history & any(is.na(life) | is.na(errorlife_tree)) ~ "missing",
      is_invalid_history ~ "error",
      birth == 0 & death == 0 ~ "left-t & right-c",
      !has_standing_dead & first_obs_is_stump ~ "from alive to stump",
      !has_standing_dead ~ "no standing dead record",
      birth == 0 ~ "left-truncated",
      is.na(birth) & first_obs_is_stump ~ "from alive to stump",
      first_death_rank > 1 & !((first_death_rank - 1) %in% rank_inventory) ~ "dead recruit",
      death == 0 ~ "right-censored",
      has_turned_stump ~ "censored (to stump)",
      TRUE ~ "censored (to decomposed)"
    )
  ) %>%
  ungroup() %>%
  
  # Strip data processing structural helpers to deliver clean final tables
  select(-is_invalid_history, -has_standing_dead, -first_obs_is_stump, -has_turned_stump)


# ==============================================================================
# SEQUENCE REVISION OF STATUS TRAJECTORIES
# ==============================================================================

trees_dead3 <- trees_dead3 %>%
  # Group by individual stem to track its precise historical development
  group_by(composed_site_id, plot_id, tree_id, stem_id) %>%
  # Crucial: Sort by rank_inventory to ensure timeline flows from past to future
  arrange(rank_inventory, .by_group = TRUE) %>%
  mutate(
    # Get the status from the PREVIOUS inventory year
    previous_status = lag(status),
    
    # Construct the status_check trajectory
    status_check = case_when(
      # If there is no previous record, it's the beginning of the tree's history
      is.na(previous_status) ~ paste0("START_", status),
      
      # Standard transition case: Combine previous and current status (e.g., "AC_DF")
      TRUE ~ paste0(previous_status, "_", status)
    )
  ) %>%
  ungroup() %>%
  # Clean up the temporary helper column
  select(-previous_status)

library(dplyr)

# ==============================================================================
# VALIDATION OF STATUS TRAJECTORIES (CREATING status_error)
# ==============================================================================

library(dplyr)


trees_dead3 <- trees_dead3 %>%
  mutate(
    errorstatus_record = case_when(
      # Special accepted cases (Minor measurement inconsistencies to be ignored)
      status_check %in% c("DF_DC", "AF_AC", "AF_DC") ~ 0,
      
      # Stump errors: If DP is at the beginning of the pair and followed by anything else
      grepl("^DP_", status_check) & !grepl("_DP$", status_check) ~ 1,
      
      # Zombie errors: Dead trees becoming alive again (DC or DF transitioning to AC, AF, or A)
      status_check %in% c("DC_AC", "DC_AF", "DC_A", 
                          "DF_AC", "DF_AF", "DF_A") ~ 1,
      
      # Broad zombie check safety net (Any D/DC/DF to A/AC/AF transition)
      grepl("^(DC_|DF_|D_)", status_check) & grepl("(_AC|_AF|_A)$", status_check) ~ 1,
      
      # All other cases are correct
      TRUE ~ 0
    )
  )

# 2. SAFE AGGREGATION TO errorstatus_tree VIA HELPER LISTS
# ==============================================================================

# A. Error list for STANDARD sites (grouped by plot_id)
bad_stems_standard <- trees_dead3 %>%
  filter(!composed_site_id %in% target_sites) %>%
  group_by(composed_site_id, plot_id, tree_id, stem_id) %>%
  summarise(has_error_standard = if_else(any(errorstatus_record == 1), 1, 0), .groups = "drop") %>%
  filter(has_error_standard == 1)

# B. Error list for VUK target sites (globally pooled, ignoring plot_id)
bad_stems_vuk <- trees_dead3 %>%
  filter(composed_site_id %in% target_sites) %>%
  group_by(composed_site_id, tree_id, stem_id) %>%
  summarise(has_error_vuk = if_else(any(errorstatus_record == 1), 1, 0), .groups = "drop") %>%
  filter(has_error_vuk == 1)

# 3. DATASET RE-ASSEMBLY & FINAL FORMATION OF errorstatus_tree
# ==============================================================================

trees_dead3 <- trees_dead3 %>%
  # Map error indicators from standard plots
  left_join(bad_stems_standard, by = c("composed_site_id", "plot_id", "tree_id", "stem_id")) %>%
  # Map error indicators from pooled VUK sites
  left_join(bad_stems_vuk, by = c("composed_site_id", "tree_id", "stem_id")) %>%
  
  # Assign final tree-level error based on where the error flags were matched
  mutate(
    errorstatus_tree = case_when(
      composed_site_id %in% target_sites & has_error_vuk == 1 ~ 1,
      !composed_site_id %in% target_sites & has_error_standard == 1 ~ 1,
      TRUE ~ 0
    )
  ) %>%
  # Clean up temporary columns from joins
  select(-has_error_standard, -has_error_vuk)


# KONTROLA BIRTH AND DEATH PODLE TAB. TIMELINE


# TASK - THRESHOLDS (e.g. exclude records)
# 

