write.table(trees_dead6, 
            file = "trees_dead6.txt", 
            sep = "\t",          
            dec = ".",           
            row.names = FALSE, 
            quote = FALSE,       
            na = "",             
            fileEncoding = "UTF-8")

write.table(fagus_1, 
            file = "fagus_1.txt", 
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

# First merge plot_id in Zofin and Ranspurk - plot_id2
# Define the 5 specific VUK locations where plots are globally pooled
target_sites <- c(
  "VUK__1__Zofin__a",
  "VUK__1__Zofin__b",
  "VUK__1__Zofin__c",
  "VUK__1__Zofin__d",
  "VUK__15__Ranspurk__NA"
)

# Create plot_id2 while preserving the original plot_id for design mapping
trees <- trees %>%
  mutate(
    plot_id2 = if_else(
      composed_site_id %in% target_sites, 
      "1", 
      as.character(plot_id)
    )
  )


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
      plot_id2,          # Joined from 'plots' table in step 2
      inventory_type,   # Joined from 'plots' table in step 2
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
    # Handle missing or zero stem IDs: convert to literal "NA" string
    stem_id_clean = if_else(is.na(stem_id) | stem_id == "" | stem_id == "0", "NA", as.character(stem_id)),
    
    # Ensure tree_id is treated as character to avoid format mismatch during concatenation
    tree_id_clean = as.character(tree_id),
    
    # Concatenate the plot base with tree and stem identifiers using an underscore separator
    tree_unique_id = paste(plot_unique_id, tree_id_clean, stem_id_clean, sep = "_")
  ) %>%
  # Remove temporary cleaning columns to keep the dataframe pristine
  select(-stem_id_clean, -tree_id_clean) %>%
  # Sort final data chronologically to maintain historical census records per stem
  arrange(tree_unique_id, inventory_year)

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
# IDENTIFY ZOMBIE TREES VIA CHRONOLOGICAL CUMULATIVE CHECKS
# ==============================================================================

trees_dead <- trees_consistent %>%
  group_by(tree_unique_id) %>%
  # Filter to keep only trees that have at least one dead (snag) record
  filter(any(life == "D")) %>%
  
  # Ensure the history is strictly ordered by real calendar years
  arrange(inventory_year, .by_group = TRUE) %>%
  
  mutate(
    # Was this specific tree already recorded as dead ("D") in any previous census?
    # cumany() turns TRUE on the first "D" and stays TRUE for all subsequent rows.
    was_dead_before = lag(cumany(life == "D"), default = FALSE),
    
    # A zombie record occurs if the tree is active ("A") but was already dead before
    errorlife_record = if_else(life == "A" & was_dead_before, 1, 0),
    
    # Flag the entire tree history if a zombie record is found anywhere
    errorlife_tree = as.integer(any(errorlife_record == 1))
  ) %>%
  ungroup() %>%
  # Clean up the temporary logical tracking column
  select(-was_dead_before)

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

# ==============================================================================
# 1. SEQUENCE REVISION OF STATUS TRAJECTORIES
# ==============================================================================

trees_dead2b <- trees_dead2 %>%
  # Group by plot_id2 space to track individual stem trajectories across unified plots
  group_by(composed_site_id, plot_id2, tree_id, stem_id) %>%
  # Sort by rank_inventory to ensure timeline flows chronologically
  arrange(rank_inventory, .by_group = TRUE) %>%
  mutate(
    # Fetch status from the previous inventory record
    previous_status = lag(status),
    
    # Construct sequential status transaction checkpoints
    status_check = case_when(
      is.na(previous_status) ~ paste0("START_", status),
      TRUE ~ paste0(previous_status, "_", status)
    )
  ) %>%
  ungroup() %>%
  # Drop the temporary lag helper column
  select(-previous_status)


# ==============================================================================
# 2. VALIDATION OF STATUS TRAJECTORIES (RECORD & TREE-LEVEL ERRORS)
# ==============================================================================

# Step A: Classify errors on individual records (rows)
trees_dead2b <- trees_dead2b %>%
  mutate(
    errorstatus_record = case_when(
      # Minor fieldwork inconsistencies allowed (to be ignored)
      status_check %in% c("DF_DC", "AF_AC", "AF_DC") ~ 0,
      
      # Stump sequence violations: DP at entry cannot transition to standing/alive
      grepl("^DP_", status_check) & !grepl("_DP$", status_check) ~ 1,
      
      # Zombie violations: Verified dead stems resurrecting as alive
      status_check %in% c("DC_AC", "DC_AF", "DC_A", 
                          "DF_AC", "DF_AF", "DF_A") ~ 1,
      
      # Comprehensive zombie safety net (any dead-to-alive progression)
      grepl("^(DC_|DF_|D_)", status_check) & grepl("(_AC|_AF|_A)$", status_check) ~ 1,
      
      # All other status transitions are valid
      TRUE ~ 0
    )
  ) %>%
  
  # Step B: Instantly propagate record errors to the entire tree history using plot_id2.
  # This completely replaces the previous complex split-and-join approach (A & B tables).
  group_by(composed_site_id, plot_id2, tree_id, stem_id) %>%
  mutate(
    errorstatus_tree = if_else(any(errorstatus_record == 1), 1, 0)
  ) %>%
  ungroup()

# ==============================================================================
# SEPARATING STUMPS (STATUS "DP") FROM THE MAIN DATASET
# ==============================================================================

# 1. Create a separate table containing ONLY stump records
trees_dead_stumps <- trees_dead2b %>%
  filter(status == "DP")

# 2. Keep only non-stump records in the main table for threshold mapping
trees_dead2c <- trees_dead2b %>%
  filter(status != "DP" | is.na(status))

# ==============================================================================
# AGAIN SEPARATING TREES WITH NO DEAD IN HISTORY FROM THE MAIN DATASET
# ==============================================================================
# GENERATING CHRONOLOGICAL STATUS SEQUENCES FOR STEM HISTORIES
# ==============================================================================

# 1. Define the 5 target sites that require plot_id unification
target_sites <- c("VUK__15__Ranspurk__NA", "VUK__16__Zofin__NA", "HUL__1__Zofin__NA", "HUL__2__Ranspurk__NA", "HUL__3__Boubin__NA")

# 2. Create the unified plot_id2 column and ensure chronological sorting
# We replace circle_no with plot_id2 in the sorting hierarchy
trees_dead2c <- trees_dead2c %>%
  mutate(
    plot_id2 = if_else(composed_site_id %in% target_sites, "1", as.character(plot_id))
  ) %>%
  arrange(composed_site_id, plot_id2, tree_id, stem_id, inventory_year)

# 3. Group by each unique stem context (using plot_id2) and collapse status history
stem_sequences <- trees_dead2c %>%
  group_by(composed_site_id, plot_id2, tree_id, stem_id) %>%
  summarise(
    sequence = paste(coalesce(status, "NA"), collapse = "_"),
    .groups = "drop"
  )

# 4. Map the sequence column back to the main dataset using the updated plot_id2 key
trees_dead2d <- trees_dead2c %>%
  left_join(
    stem_sequences,
    by = c("composed_site_id", "plot_id2", "tree_id", "stem_id")
  )

# FILTERING OF ISOLATED SEQUENCES AND COMPLETELY MISSING 'D' STATUSES
# ==============================================================================

# 1. Define the specific isolated sequences to be discarded
target_sequences <- c("A", "AC", "AF")

# 2. Apply both filters in a single pipeline to create trees_dead4b
trees_dead2e <- trees_dead2d_clean %>%
  # Filter A: Drop exact, isolated matches only (e.g., leaves "AC_DC" untouched)
  filter(!sequence %in% target_sequences) %>%
  
  # Filter B: Drop histories that never entered a standing dead phase (lack 'D')
  filter(stringr::str_detect(sequence, "D"))




# ==============================================================================
# DBH THRESHOLD MAPPING
# ==============================================================================

# Clean the design dataset, ensuring thresholds exist, and convert cm to mm.
# Include circle_no in the selection for strict spatial partitioning.
design_clean <- design %>%
  filter(!is.na(standing_alive_threshold), !is.na(standing_dead_threshold)) %>%
  mutate(
    a_threshold_mm = standing_alive_threshold * 10,
    d_threshold_mm  = standing_dead_threshold * 10
  ) %>%
  select(composed_site_id, inventory_year, inventory_type, inventory_id, circle_no, plots_list, a_threshold_mm, d_threshold_mm)

design_long <- design_clean %>%
  # 1. Clean the string to remove brackets and quotes: ["1", "2"] -> 1, 2
  mutate(plot_id = stringr::str_remove_all(plots_list, '[\"\\[\\]]')) %>%
  # 2. Split the comma-separated string into separate rows
  tidyr::separate_rows(plot_id, sep = ",\\s*") %>%
  # 3. Trim any potential accidental whitespace
  mutate(plot_id = stringr::str_trim(plot_id)) %>%
  # 4. Remove the helper column as it is no longer needed
  select(-plots_list)


# STANDARD SEAMLESS JOIN
# ------------------------------------------------------------------------------

trees_dead2e <- trees_dead2e %>%
  left_join(
    design_long,
    by = c("composed_site_id", "inventory_year", "inventory_id", "inventory_type", "circle_no", "plot_id")
  )

library(dplyr)

# CALCULATING MAXIMUM HISTORICAL DBH AND RE-EVALUATING THRESHOLDS
# ==============================================================================

trees_dead2e <- trees_dead2e %>%
  # 1. Group by the unique stem identifier
  group_by(composed_site_id, plot_id2, tree_id, stem_id) %>%
  
  # 2. Calculate max DBH safely by checking if any valid DBH exists first
  mutate(
    tree_dbh_max = if_else(
      condition = any(!is.na(dbh)),       # Is there at least one non-NA value in the group?
      true      = max(dbh, na.rm = TRUE),  # If YES, calculate the maximum
      false     = NA_real_                 # If NO (all are NA), assign NA straight away
    )
  ) %>%
  
  # 3. Ungroup to return to a standard flat table structure
  ungroup()

# ==============================================================================
# EVALUATION OF DBH AGAINST THRESHOLDS
# ==============================================================================

trees_dead3 <- trees_dead2e %>%
  mutate(
    # 1. Check DBH thresholds for ALIVE records (life == "A")
    a_thresh_check = case_when(
      life == "A" & tree_dbh_max >= a_threshold_mm ~ "within",
      life == "A" & tree_dbh_max <  a_threshold_mm ~ "under",
      TRUE ~ NA_character_  # Returns NULL (NA) for life == "D" or any other cases
    ),
    
    # 2. Check DBH thresholds for DEAD records (life == "D")
    d_thresh_check = case_when(
      life == "D" & tree_dbh_max >= d_threshold_mm ~ "within",
      life == "D" & tree_dbh_max <  d_threshold_mm ~ "under",
      TRUE ~ NA_character_  # Returns NULL (NA) for life == "A" or any other cases
    )
  )


# ==============================================================================
# STEM-LEVEL FILTERING FOR DISPARATE LIVE/DEAD THRESHOLDS
# ==============================================================================

# Step 1: Identify the UNIQUE KEYS of stems that failed the threshold 
# on ANY single row in their history (either live or dead check)
problematic_stem_keys <- trees_dead3 %>%
  filter(a_thresh_check == "under" | d_thresh_check == "under") %>%
  distinct(composed_site_id, plot_id2, tree_id, stem_id)


# Step 2: Extract ALL historical rows for these problematic stems.
# If a stem is under the limit as a dead tree, we pull its live history here too.
trees_dead_under <- trees_dead3 %>%
  semi_join(
    problematic_stem_keys, 
    by = c("composed_site_id", "plot_id2", "tree_id", "stem_id")
  )


# Step 3: Create the final clean analytical dataset (Strictly "within" always).
# We completely purge the main dataset of any stems found in the 'under' table.
trees_dead_thres <- trees_dead3 %>%
  anti_join(
    trees_dead_under, 
    by = c("composed_site_id", "plot_id2", "tree_id", "stem_id")
  )

# DECIDING ABOUT THRESHOLDS AND FOREST TYPES GROUPS
trees_dead4_thres2 <- trees_dead4_thres %>%
  left_join(
    wildcard_metadata,
    by = c("composed_site_id")
  )

trees_dead4_thres2 <- trees_dead4_thres %>%
  left_join(
    # Z tabulky wildcard_metadata vybereme klíč + jen ty sloupce, které chceme
    wildcard_metadata %>% select(composed_site_id, soil_water, soil_nutrient, fortype, EEA_fortype, ecoregion),
    by = c("composed_site_id")
  )


# V QGIS z rastru připojit teploty, srážky, nadmořskou výšku (centroid nebo polygon)
# Pak promýšlet environmntální proměnné (diskuze)
# U starých censů může být problém, že nemáme průměrné teploty a srážky z té doby
# - možno nadmořská výška?
# Pro finální výběr lokalit získat jejich souřadnice (centroid)


# ==============================================================================
# ==============================================================================
# CREATING BIRTH AND DEATH
# ==============================================================================
# ==============================================================================




# CONSTRUCTING PLOT TIMELINE FROM TREE-LEVEL DATASET
# CALCULATING CENSUS RANK AND BIRTH/DEATH (WITH 0 FOR UNKNOWN)
# ==============================================================================
# In trees_3inv table adjusted plot_id2 is necessary
plot_timeline <- trees_3inv %>% 
  # 1. Select relevant tracking columns directly from trees_3inv
  select(institute, composed_site_id, plot_id2, inventory_year) %>%
  
  # 2. Reduce tree-level data to unique plot-year combinations
  distinct(institute, composed_site_id, plot_id2, inventory_year) %>%
  filter(!is.na(inventory_year)) %>%
  
  # 3. Group and arrange chronologically per plot
  group_by(composed_site_id, plot_id2) %>%
  arrange(inventory_year, .by_group = TRUE) %>%
  
  # 4. Calculate all temporal variables, rank, and BaSTA-compliant birth/death
  mutate(
    no_inventories   = n_distinct(inventory_year),
    rank_inventory   = row_number(),
    
    # Lag and Lead functions to identify neighboring census years
    inventory_before = lag(inventory_year, order_by = inventory_year),
    inventory_next   = lead(inventory_year, order_by = inventory_year),
    
    # Birth mid-point (0 if it is the first inventory)
    birth = if_else(is.na(inventory_before), 0, (inventory_year + inventory_before) / 2),
    
    # Death mid-point (0 if it is the last inventory)
    death = if_else(is.na(inventory_next), 0, (inventory_year + inventory_next) / 2)
  ) %>%
  
  # 5. Clean up grouping and apply final sorting
  ungroup() %>%
  arrange(institute, composed_site_id, plot_id2, inventory_year) 


# In the next step, table was adjusted in table outside of R - for two overlapping plots
# in VUK Ranspurk and Zofin
# CALCULATE BIRTH AND DEATH MIDPOINTS FOR BIG AND SMALL PLOTS IN PLOT_TIMELINE
# Applies for Ranspurk and Zofin
# ==============================================================================

plot_timeline_updated <- plot_timeline_new %>%
  mutate(
    # --- 1. BIG PLOTS (Four-census variant) ---
    # Birth midpoint: average of current year and the previous valid census for BIG
    birth_big = if_else(
      is.na(inventory_before_big) | inventory_before_big == "", 
      0, 
      (inventory_year + as.numeric(inventory_before_big)) / 2
    ),
    
    # Death midpoint: average of current year and the next valid census for BIG
    death_big = if_else(
      is.na(inventory_next_big) | inventory_next_big == "", 
      0, 
      (inventory_year + as.numeric(inventory_next_big)) / 2
    ),
    
    # --- 2. SMALL PLOTS (Six-census variant) ---
    # Birth midpoint: average of current year and the previous valid census for SMALL
    birth_small = if_else(
      is.na(inventory_before_small) | inventory_before_small == "", 
      0, 
      (inventory_year + as.numeric(inventory_before_small)) / 2
    ),
    
    # Death midpoint: average of current year and the next valid census for SMALL
    death_small = if_else(
      is.na(inventory_next_small) | inventory_next_small == "", 
      0, 
      (inventory_year + as.numeric(inventory_next_small)) / 2
    )
  )

# ADJUSTING BIRTH AND DEATH
# ==============================================================================
# Step 1: Create the 'big_small' column on the tree level
# ------------------------------------------------------------------------------

trees_dead5b <- trees_dead5 %>%
  group_by(tree_unique_id) %>%
  mutate(
    # A. Basic site and year tracking
    is_ranspurk      = any(grepl("Ranspurk", tree_unique_id)),
    has_ranspurk_yrs = any(inventory_year %in% c(2015, 2025)),
    
    is_zofin         = any(grepl("Zofin", tree_unique_id)),
    has_zofin_yrs    = any(inventory_year %in% c(2012, 2017)),
    
    # B. DBH tracking for Ranšpurk edge-cases (trees only recorded in 2020)
    # Check if the tree has ONLY a 2020 record and belongs to the small DBH group
    is_only_2020_small_dbh = any(inventory_year == 2020) & 
      n_distinct(inventory_year) == 1 & 
      any(dbh_group == "10-99mm"), # Adjust column name if needed
    
    # C. Comprehensive classification logic
    big_small = case_when(
      # If it's Ranšpurk and is small DBH found only in 2020 -> must be "small"
      is_ranspurk & is_only_2020_small_dbh ~ "small",
      
      # Standard Ranšpurk rules
      is_ranspurk & has_ranspurk_yrs        ~ "small",
      is_ranspurk & !has_ranspurk_yrs       ~ "big",
      
      # Standard Žofín rules
      is_zofin & has_zofin_yrs              ~ "small",
      is_zofin & !has_zofin_yrs              ~ "big",
      
      # Default for all other standard sites
      TRUE                                  ~ "small"
    )
  ) %>%
  ungroup() %>%
  # Clean up temporary logical flags
  select(-is_ranspurk, -has_ranspurk_yrs, -is_zofin, -has_zofin_yrs, -is_only_2020_small_dbh)

# STEP 2: PREPARE TIMELINE LOOKUPS DIRECTLY VIA CALENDAR YEARS
# ==============================================================================
# We extract clean lookups from 'plot_timeline_updated' using years as keys.
# We ensure year columns are numeric and remove any row duplicates caused by sub-plots.

# 2A. Lookup table for BIG plots (fetching birth_big and death_big)
timeline_big <- plot_timeline_updated %>%
  filter(!is.na(inventory_year_big) & inventory_year_big != "") %>%
  mutate(inventory_year = as.numeric(inventory_year_big)) %>%
  select(composed_site_id, inventory_year, birth = birth_big, death = death_big) %>%
  distinct(composed_site_id, inventory_year, .keep_all = TRUE)

# 2B. Lookup table for SMALL plots (fetching birth_small and death_small)
timeline_small <- plot_timeline_updated %>%
  filter(!is.na(inventory_year_small) & inventory_year_small != "") %>%
  mutate(inventory_year = as.numeric(inventory_year_small)) %>%
  select(composed_site_id, inventory_year, birth = birth_small, death = death_small) %>%
  distinct(composed_site_id, inventory_year, .keep_all = TRUE)


# ==============================================================================
# STEP 3: SPLIT DATASETS, CALCULATE FIRST/LAST SNAG YEARS, AND JOIN
# ==============================================================================
# Data is separated by 'big_small'. For each tree, we find its absolute first and 
# last calendar year as a snag ("D"). These fixed years are used for the timeline join.

# --- 3A. PROCESS BIG PLOTS ---
trees_big_joined <- trees_dead5_prepped %>%
  filter(big_small == "big") %>%
  group_by(tree_unique_id) %>%
  mutate(
    # Find the absolute first and last calendar year where this tree was a snag ("D")
    first_d_year = min(inventory_year[life == "D"]),
    last_d_year  = max(inventory_year[life == "D"])
  ) %>%
  ungroup() %>%
  # Join BIRTH: Based on the FIRST year it became a snag
  left_join(
    timeline_big %>% select(composed_site_id, inventory_year, birth), 
    by = c("composed_site_id", "first_d_year" = "inventory_year")
  ) %>%
  # Join DEATH: Based on the LAST year it was seen standing
  left_join(
    timeline_big %>% select(composed_site_id, inventory_year, death), 
    by = c("composed_site_id", "last_d_year" = "inventory_year")
  )


# --- 3B. PROCESS SMALL PLOTS ---
trees_small_joined <- trees_dead5_prepped %>%
  filter(big_small == "small") %>%
  group_by(tree_unique_id) %>%
  mutate(
    # Find the absolute first and last calendar year where this tree was a snag ("D")
    first_d_year = min(inventory_year[life == "D"]),
    last_d_year  = max(inventory_year[life == "D"])
  ) %>%
  ungroup() %>%
  # Join BIRTH: Based on the FIRST year it became a snag
  left_join(
    timeline_small %>% select(composed_site_id, inventory_year, birth), 
    by = c("composed_site_id", "first_d_year" = "inventory_year")
  ) %>%
  # Join DEATH: Based on the LAST year it was seen standing
  left_join(
    timeline_small %>% select(composed_site_id, inventory_year, death), 
    by = c("composed_site_id", "last_d_year" = "inventory_year")
  )


# ==============================================================================
# STEP 4: RECOMBINE AND PRESERVE NA VALUES FOR QUALITY CONTROL
# ==============================================================================
# We merge the datasets back together and sort them chronologically.

trees_dead5_final <- bind_rows(trees_big_joined, trees_small_joined) %>%
  # Arrange chronologically to maintain historical sequence per tree stem
  arrange(tree_unique_id, inventory_year)
# Assign birth and death according to life change
# ==============================================================================


# TIME_TYPE (left-truncated, ...)
# ==============================================================================


# ==============================================================================
# EXTRACT INITIAL DBH AT FIRST SNAG DETECTION (DBH_FIRST)
# ==============================================================================
# Code is still not working, need to fix it. A lot of trees with numeric dbh
# receives dbh_first NA

trees_dead6 <- trees_dead5_final %>%
  group_by(tree_unique_id) %>%
  mutate(
    # 1. Exact DBH at the first snag year (returns NA if not found)
    exact_dbh = na.omit(if_else(inventory_year == first_d_year, dbh, NA_real_))[1],
    
    # 2. Closest DBH BEFORE first snag year (we take the LAST/LATEST one available)
    prev_dbh  = last(na.omit(if_else(inventory_year < first_d_year, dbh, NA_real_))),
    
    # 3. Closest DBH AFTER first snag year (we take the FIRST/EARLIEST one available)
    next_dbh  = first(na.omit(if_else(inventory_year > first_d_year, dbh, NA_real_)))
  ) %>%
  # 4. Combine into final dbh_first using the requested hierarchy
  mutate(
    dbh_first = case_when(
      is.na(first_d_year) ~ NA_real_,
      !is.na(exact_dbh)   ~ exact_dbh,
      !is.na(prev_dbh)    ~ prev_dbh,
      !is.na(next_dbh)    ~ next_dbh,
      TRUE                ~ NA_real_
    )
  ) %>%
  # Clean up temporary columns and ungroup
  select(-exact_dbh, -prev_dbh, -next_dbh) %>%
  ungroup()

# ==============================================================================
# IDENTIFY SPECIES STABILITY AND CLEAN UNKNOWN TAXONS PER STEM
# ==============================================================================

trees_dead6 <- trees_dead6 %>%
  group_by(tree_unique_id) %>%
  mutate(
    # 1. Gather all unique scientific names recorded for this specific stem history,
    # completely ignoring generic "Unknown" placeholders and NA values.
    valid_species_list = list(unique(full_scientific[!full_scientific %in% c("Unknown broadleaf", "Unknown species", "Unknown conifer", NA, "")])),
    
    # 2. Count how many distinct, valid species names exist in the stem's timeline
    distinct_real_species_count = length(valid_species_list[[1]]),
    
    # 3. Apply classification logic to handle taxonomy issues over time:
    tree_species = case_when(
      # If more than one valid species name is found, the taxonomy changed (error)
      distinct_real_species_count > 1 ~ "changed",
      
      # If exactly one valid species name is found, assign it to the entire stem history
      distinct_real_species_count == 1 ~ valid_species_list[[1]][1],
      
      # If the stem never had a specific name and was always unknown, flag as Unknown
      TRUE ~ "Unknown"
    )
  ) %>%
  ungroup() %>%
  # Remove temporary list and counter columns to keep the dataframe clean
  select(-valid_species_list, -distinct_real_species_count)


# Kontrola birth a death podle plot_timeline, jestli tam nejsou nějaké nesmysly 
# (záhadné birth a death)

# Test na přeskakující stromy


