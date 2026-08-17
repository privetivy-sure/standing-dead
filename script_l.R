# R SCRIPT FOR FILTERING AND ADJUSTING LYING STEMS FOR BASTA ANALYSIS
# ==============================================================================

# Basic functions
# ------------------------------------------------------------------------------


write.table(sites_clima, 
            file = "sites_clima2.txt", 
            sep = "\t",          
            dec = ".",           
            row.names = FALSE, 
            quote = FALSE,       
            na = "",             
            fileEncoding = "UTF-8")

write.table(trees_50_new3, 
            file = "trees_50_new.txt", 
            sep = "\t",          
            dec = ".",           
            row.names = FALSE, 
            quote = FALSE,       
            na = "",             
            fileEncoding = "UTF-8")

library(dplyr)
# Deletes columns
trees_50_new3 <- trees_50_new2 %>% select(-full_scientific, -drevina, -tree_species)

# Creates a list with numbers of combinations
kontingence_long <- trees_50_new2 %>%
  count(full_scientific)


# Removes rows with particular value in selected column
trees_50_surv3 <- trees_50_surv2 %>%
  filter(cons_recruits_error != 1)


# Creates a list
site_characteristics <- trees_50_new %>%
  select(composed_site_id, site_name, rain, temp) %>%
  distinct()

# ==============================================================================
# Project: BaSTA (Survival / Persistence) Analysis for Lying Deadwood
# DATA PREPARATION: UNIQUE IDs AND INVENTORY FILTERING
# 
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
dbListFields(con, "metadata")

library(RPostgres)


treesL <- dbGetQuery(con, "
  SELECT * FROM trees 
  WHERE position = 'L'") 

plots <- dbGetQuery(con, "
  SELECT * FROM plots
")

metadata <- dbGetQuery(con, "
  SELECT * 
  FROM metadata
  WHERE institute IN ('NWFVA', 'LFRLP-FAWF', 'DISAFA - UNITO')
")

metadata <- tbl(con, "metadata") %>%
  filter(institutes %in% c("NWFVA", "LFRLP-FAWF", "DISAFA - UNITO")) %>%
  collect()

library(dplyr)

# --- JOIN TABLES ---
# Join selected columns from 'plots' into 'trees' based on 'plot_record_id'
# left_join ensures we keep all records in the 'trees' table
treesL <- treesL %>%
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
# ==============================================================================

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
treesL <- treesL %>%
  mutate(
    plot_id2 = if_else(
      composed_site_id %in% target_sites, 
      "1", 
      as.character(plot_id)
    )
  )


treesL <- treesL %>%
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

# Verify missing values in the final unique identifier
total_nas <- sum(is.na(treesL$plot_unique_id))
cat("\nNumber of NA values in tree_unique_id:", total_nas, "\n")

# --- CREATE TREE_UNIQUE_ID ----
# Handling stem_id errors and generating unique stem-level identifiers
# ==============================================================================

library(dplyr)

# STANDARDIZE PIECE_ID TO LOWERCASE AND CHECK FOR SUSPICIOUS / INVALID IDs
# ------------------------------------------------------------------------------

# 1. Convert uppercase letters in piece_id to lowercase
treesL <- treesL %>%
  mutate(
    piece_id = tolower(as.character(piece_id))
  )

# 2. Check for suspicious or invalid IDs in stem_id and piece_id
stem_errors <- treesL %>%
  filter(is.na(stem_id) | stem_id %in% c("NA", "", 0, "0"))

piece_errors <- treesL %>%
  filter(is.na(piece_id) | piece_id %in% c("NA", "", 0, "0", "na"))

if (nrow(stem_errors) > 0) {
  cat("⚠️ Found", nrow(stem_errors), "suspicious stem_id records. These will be standardized to 'NA' text.\n")
}
if (nrow(piece_errors) > 0) {
  cat("⚠️ Found", nrow(piece_errors), "suspicious piece_id records. These will be standardized to 'NA' text.\n")
}
# ------------------------------------------------------------------------------
# STANDARDIZE IDs AND CONCATENATE WITH LITERAL "NA" STRINGS
# ------------------------------------------------------------------------------
treesL <- treesL %>%
  mutate(
    # Convert missing/zero stem IDs to the literal string "NA"
    stem_id_clean = if_else(
      is.na(stem_id) | stem_id %in% c("", "0", "NA"), 
      "NA", 
      as.character(stem_id)
    ),
    
    # Convert missing/zero piece IDs to the literal string "NA"
    piece_id_clean = if_else(
      is.na(piece_id) | piece_id %in% c("", "0", "NA"), 
      "NA", 
      as.character(piece_id)
    ),
    
    # Ensure tree_id is character to avoid concatenation issues
    tree_id_clean = as.character(tree_id),
    
    # Concatenate all components into plot_tree_stem_piece
    # Missing components will explicitly appear as "_NA" in the string
    tree_unique_id = paste(plot_unique_id, tree_id_clean, stem_id_clean, piece_id_clean, sep = "_")
  ) %>%
  # Remove temporary cleaning columns
  select(-stem_id_clean, -piece_id_clean, -tree_id_clean) %>%
  # Sort final dataset chronologically by tree identifier and census year
  arrange(tree_unique_id, inventory_year)

# ------------------------------------------------------------------------------
# VERIFY MISSING VALUES IN THE FINAL UNIQUE IDENTIFIER
# ------------------------------------------------------------------------------
total_nas <- sum(is.na(treesL$tree_unique_id))
cat("\nNumber of TRUE NA values in tree_unique_id:", total_nas, "\n")


# ==============================================================================
# INVENTORY SUMMARY & FILTERING (>= 3 INVENTORIES) ----
# Identifying plots with enough temporal depth for survival analysis
# ==============================================================================

# Calculate number of inventories directly in the trees table
treesL <- treesL %>%
  group_by(plot_unique_id) %>%
  mutate(
    no_inventories = n_distinct(inventory_year),
    is_3inv = if_else(no_inventories >= 3, "Y", "N")
  ) %>%
  ungroup()

# --- Create the summary table for reporting ---
# Including plot metadata for easier filtering and identification
summary_inventory <- treesL %>%
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
trees_3inv <- treesL %>%
  filter(is_3inv == "Y") %>%
  # Sort by plot, tree, stem, and year
  # Using the newly created tree_unique_id for perfectly organized records
  arrange(tree_unique_id, inventory_year)

cat("\n--- DATA VOLUME CHECK ---\n")
cat("Total tree records in original table: ", nrow(treesL), "\n")
cat("Tree records in filtered table (3inv):", nrow(trees_3inv), "\n")


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
# ------------------------------------------------------------------------------
library(dplyr)
library(tidyr)

# --- CREATE 'consistent_tree' COLUMN (Tree history logic) ---
# Rule: If a tree is flagged with consistent_id == FALSE in any census,
# its entire history (all rows for that tree_unique_id) is marked as FALSE.
# Single-census trees (with NA/NULL in consistent_id) are kept as TRUE.
trees_3inv <- trees_3inv %>%
  group_by(tree_unique_id) %>%
  mutate(
    consistent_tree = !any(consistent_id %in% FALSE)
  ) %>%
  ungroup()


# --- FINAL FILTRATION - CREATING trees_consistent TABLE ---
# Keep all trees EXCEPT those explicitly marked as inconsistent in any year
trees_consistent <- trees_3inv %>%
  filter(consistent_tree)


library(dplyr)


# ==============================================================================
# DETECT AND MARK DUPLICATE MEASUREMENTS WITHIN THE SAME YEAR
# ==============================================================================
# Applies for Zofin, Ranspurk and Mionsi where 2 plots are overlapping
# And in certain inventory_year trees are twice in the DB
# Removing trees from smaller plots
# Checking the whole dataset
library(dplyr)

trees_consistent <- trees_consistent %>%
  # 1. Group by tree unique ID and inventory year
  group_by(tree_unique_id, inventory_year) %>%
  
  # 2. Prioritize dataset removal logic within groups
  # Priority 2 = target datasets to mark as duplicate (sorted last)
  # Priority 1 = all other datasets (kept first)
  mutate(
    dataset_priority = case_when(
      site_name == "Mionsi" & dataset == "cm"                      ~ 2,
      site_name %in% c("Ranspurk", "Zofin") & dataset == "c2b"    ~ 2,
      TRUE                                                        ~ 1
    )
  ) %>%
  
  # 3. Sort so priority 1 comes first, followed by priority 2, then plot_id
  arrange(dataset_priority, plot_id, .by_group = TRUE) %>%
  
  # 4. Flag duplicate rows (1st record = FALSE, 2nd+ record = TRUE)
  mutate(
    duplicity = row_number() > 1
  ) %>%
  
  # Clean up temporary priority column and ungroup
  select(-dataset_priority) %>%
  ungroup()


# Summary of duplicate rows count for all composed_site_id with at least one duplicate
duplicity_by_site_report <- trees_consistent %>%
  filter(duplicity == TRUE) %>%
  count(institute, composed_site_id, name = "duplicate_rows_count", sort = TRUE)

# Print the report
print(duplicity_by_site_report, n = Inf)

# Filter out duplicate records if you want to remove them completely
trees_consistent <- trees_consistent %>% 
  filter(duplicity == FALSE)



# DETECT NON-STANDARD & MISSING 'life' VALUES (nalife)
# ==============================================================================

# Create a summary table of rows where 'life' is not purely "A" or "D"
non_standard_life <- trees_consistent %>%
  filter(is.na(life) | !life %in% c("A", "D")) %>%
  select(
    institute, 
    composed_site_id, 
    plot_id, 
    inventory_year, 
    tree_unique_id, 
    life
  ) %>%
  # Replace NA with explicit string "NA_value" for better table display/grouping
  mutate(life = if_else(is.na(life), "NA_value", as.character(life)))


# ==============================================================================
# FLAG NA VALUES AND CATEGORIZE DEVELOPMENT PATTERNS
# ==============================================================================
# Attribute nalife

trees_consistent2 <- trees_consistent %>%
  group_by(tree_unique_id) %>%
  # Ensure chronological ordering by census year
  arrange(inventory_year, .by_group = TRUE) %>%
  mutate(
    # TRUE for the specific census record if 'life' is missing or empty
    nalife_record = is.na(life) | life == "" | life == "NA",
    
    # TRUE for the entire tree history if 'life' is missing in ANY census
    nalife_tree = any(nalife_record),
    
    # --- DEVELOPMENT PATTERNS FOR NA TREES ---
    nalife_pattern = if_else(
      !nalife_tree, 
      "NO_NA",
      case_when(
        # 1. Single record tree which is NA
        n() == 1 ~ "ONLY_NA",
        
        # 2. Started with valid status (A/D) and developed into NA
        !nalife_record[1] & nalife_record[n()] ~ "VALID_TO_NA",
        
        # 3. Started with NA and developed into valid status (A/D)
        nalife_record[1] & !nalife_record[n()] ~ "NA_TO_VALID",
        
        # 4. NA is in the middle of history or both start/end are NA with valid data inside
        TRUE ~ "MIXED_NA"
      )
    )
  ) %>%
  ungroup()


# ==============================================================================
# IDENTIFY ZOMBIE TREES VIA CHRONOLOGICAL CUMULATIVE CHECKS (ALL TREES KEPT)
# ==============================================================================
# Attribute errorlife

trees_consistent2 <- trees_consistent2 %>%
  group_by(tree_unique_id) %>%
  # Ensure the history is strictly ordered by real calendar years
  arrange(inventory_year, .by_group = TRUE) %>%
  
  mutate(
    # Was this specific tree already recorded as dead ("D") in any previous census?
    # cumany() turns TRUE on the first "D" and stays TRUE for all subsequent rows.
    was_dead_before = lag(cumany(life == "D"), default = FALSE),
    
    # A zombie record occurs if the tree is active ("A") but was already dead before
    errorlife_record = life == "A" & was_dead_before,
    
    # Flag the entire tree history if a zombie record is found anywhere (logical TRUE/FALSE)
    errorlife_tree = any(errorlife_record)
  ) %>%
  ungroup() %>%
  # Clean up the temporary logical tracking column
  select(-was_dead_before)

library(dplyr)
library(tidyr)

# ==============================================================================
# SUMMARY OF TREE REPETITIONS AND INVENTORY TIMESPAN PER PLOT
# ==============================================================================
# Because some plots have all trees with only one repetition which is suspicious
# They need to be deleted

# 1. Calculate the number of census years per individual tree
tree_counts <- trees_consistent2 %>%
  group_by(plot_unique_id, tree_unique_id) %>%
  summarise(
    n_years = n_distinct(inventory_year),
    .groups = "drop"
  )

# 2. Pivot tree repetition counts into dynamic columns (1x, 2x, 3x, ...)
tree_repetitions_wide <- tree_counts %>%
  mutate(rep_category = paste0(n_years, "x")) %>%
  group_by(plot_unique_id, rep_category) %>%
  summarise(
    tree_count = n_distinct(tree_unique_id),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = rep_category,
    values_from = tree_count,
    values_fill = 0
  )

# 3. Calculate plot-level inventory statistics (no_inventories & inventories_span)
plot_inventory_stats <- trees_consistent2 %>%
  group_by(plot_unique_id) %>%
  summarise(
    # Taking the maximum value of no_inventories recorded for the plot
    no_inventories = max(no_inventories, na.rm = TRUE),
    
    # Timespan between the latest and earliest inventory year
    inventories_span = max(inventory_year, na.rm = TRUE) - min(inventory_year, na.rm = TRUE),
    
    # Total unique trees on the plot for validation
    total_unique_trees = n_distinct(tree_unique_id),
    .groups = "drop"
  )

# 4. Combine all metrics into a single final summary table
plot_tree_repetition_summary <- plot_inventory_stats %>%
  left_join(tree_repetitions_wide, by = "plot_unique_id") %>%
  # Arrange columns: plot_unique_id, plot stats, then repetition count columns (1x, 2x, ...)
  select(
    plot_unique_id, 
    no_inventories, 
    inventories_span, 
    total_unique_trees, 
    everything()
  )

# DISPLAY & EXPORT
# ==============================================================================

# Print preview of the result
print(head(plot_tree_repetition_summary))

# Export summary table to text file
write.table(
  plot_tree_repetition_summary, 
  "plot_tree_repetition_summary2.txt", 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)

# Removing all plots of BGD-NP, FVA-BW, INBO, NWFVA, TUZVO
# Removing trees from site_name "Mionsi" dataset "c1" (inventories only 2004, 2014)
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# because these plots have all dead lying trees with only one measurement (inventory)

library(dplyr)

# List of institutes to remove
excluded_institutes <- c("BGD-NP", "FVA-BW", "INBO", "NWFVA", "TUZVO")

# Filter out the specified institutes
trees_consistent3 <- trees_consistent2 %>%
  filter(!institute %in% excluded_institutes)

# Filter out part of Mionsi
remove_mionsi_c1 <- function(df) {
  df %>%
    filter(!(site_name == "Mionsi" & dataset == "c1"))
}
trees_consistent3 <- remove_mionsi_c1(trees_consistent3)

# Verification summary
cat("--- FILTERING SUMMARY ---\n")
cat("Original rows (trees_consistent2):", nrow(trees_consistent2), "\n")
cat("Filtered rows (trees_consistent3):", nrow(trees_consistent3), "\n")


# ==============================================================================
# EXTRACTION OF STATUS2 FROM EXTENDED_ATTRIBUTES
# ==============================================================================
library(stringr)

trees_consistent3 <- trees_consistent3 %>%
  mutate(
    # Extracts exact 'status2' value from JSON/JSON-like extended_attributes string
    status2 = str_match(
      extended_attributes, 
      "[\"']?status2[\"']?\\s*:\\s*[\"']?([^\"',}]+)[\"']?"
    )[, 2],
    
    # Clean up quotes and surrounding whitespace
    status2 = str_remove_all(status2, "[\"']"),
    status2 = trimws(status2)
  )

# CREATE SUMMARY BY INSTITUTE AND JOIN METADATA
# ------------------------------------------------------------------------------
institute_status2_summary <- trees_consistent3 %>%
  filter(!is.na(status2), status2 != "") %>%
  group_by(institute, status2) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  # Join with metadata table filtered specifically for 'status2' attribute
  left_join(
    metadata %>% 
      filter(attribute == "status2") %>% 
      select(institute, value_code, value_description) %>% 
      distinct(),
    # Map 'institute' -> 'institutes' AND 'status2' -> 'value_code'
    by = c("institute" = "institute", "status2" = "value_code")
  ) %>%
  arrange(institute, status2)

# PRINT SUMMARY & EXPORT TO TEXT FILE
# ------------------------------------------------------------------------------
cat("--- SUMMARY OF EXTRACTED STATUS2 VALUES WITH DESCRIPTIONS ---\n")
print(institute_status2_summary)

# Export summary table for review
write.table(
  institute_status2_summary, 
  "status2_raw_summary.txt", 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)

library(dplyr)
library(stringr)


# ==============================================================================
# CREATE SIMPLIFIED 'status3' COLUMN BASED ON 'status2' AND 'life' MAPPING
# ==============================================================================

library(dplyr)
library(stringr)

trees_consistent_L <- trees_consistent_L %>%
  mutate(
    # Clean whitespace and convert columns to character for robust matching
    status2_clean   = trimws(as.character(status2)),
    life_clean      = trimws(as.character(life)),
    institute_clean = trimws(as.character(institute)),
    
    status3 = case_when(
      # 1. Map to "NA" if life is empty/NA AND status2 is missing, empty, or contains \N
      (is.na(life_clean) | life_clean == "" | life_clean == "NA") & 
        (is.na(status2_clean) | status2_clean == "" | status2_clean == "NA" | str_detect(status2_clean, "^\\\\*N$")) ~ "NA",
      
      # 2a. Map DP to "O" ONLY for LFRLP-FAWF institute
      status2_clean == "DP" & institute_clean == "LFRLP-FAWF" ~ "O",
      
      # 2b. Map other stump and missing codes to "O"
      status2_clean %in% c("AP", "DH", "MIS", "MISABS", "MISCWD", "MISDEC") ~ "O",
      
      # 3. Map alive-related codes OR life == "A" to "A"
      status2_clean %in% c("AL", "ALU") | life_clean == "A" ~ "A",
      
      # 4. Map all other remaining values (including dead codes and DP from other institutes) to "D"
      TRUE ~ "D"
    )
  ) %>%
  # Remove temporary helper columns
  select(-status2_clean, -life_clean, -institute_clean)



library(dplyr)
library(readr)
library(stringr)



# ==============================================================================
# CREATE CHRONOLOGICAL STATUS3_SEQUENCE PER TREE HISTORY
# ==============================================================================

trees_consistent_L3 <- trees_consistent_L3 %>%
  # 1. Group by individual tree identity
  group_by(tree_unique_id) %>%
  
  # 2. Sort chronologically by measurement year
  arrange(inventory_year, .by_group = TRUE) %>%
  
  # 3. Construct sequence string from status3
  mutate(
    status3_sequence = {
      # Ensure status3 values are stored as string (replace actual NAs with "NA" text)
      clean_s3 <- if_else(is.na(status3), "NA", as.character(status3))
      
      # Check if at least one valid (non-NA) status3 value exists for this tree
      has_valid_val <- any(clean_s3 != "NA")
      
      # Concatenate values with underscore if valid record exists, else return NA
      if (has_valid_val) {
        paste(clean_s3, collapse = "_")
      } else {
        NA_character_
      }
    }
  ) %>%
  ungroup()


# ==============================================================================
# VALIDATION OF STATUS TRAJECTORIES (RECORD & TREE-LEVEL ERRORS)
# ==============================================================================
# Attribute statuserror

library(dplyr)
library(stringr)

# DETECT INVALID SEQUENCES WHERE 'O' COMES BEFORE 'D'
# ------------------------------------------------------------------------------

trees_consistent_L3 <- trees_consistent_L3 %>%
  # 1. Create statuserror_record for each individual row/measurement
  mutate(
    statuserror_record = if_else(
      !is.na(status3_sequence),
      str_detect(status3_sequence, "(^|_)O_.*D($|_)"),
      FALSE
    )
  ) %>%
  
  # 2. Propagate the error to statuserror_tree across the entire tree group
  group_by(tree_unique_id) %>%
  mutate(
    statuserror_tree = any(statuserror_record, na.rm = TRUE)
  ) %>%
  ungroup()


# ==============================================================================
# SEPARATING STUMPS, MIS AND ALIVE WHICH ARE THE WHOLE OBJECT TIMESPAN
# ==============================================================================

library(dplyr)

# Define vector of status3_sequence values targeted for removal
sequences_to_remove <- c(
  "A", "A_A", "A_A_A", "A_A_O", "A_O", 
  "NA_NA_O", "O", "O_D_O", "O_O", "O_O_O"
)

# Add 'remove' column based on exact sequence match
trees_consistent_L3 <- trees_consistent_L3 %>%
  mutate(
    discard_nodead = status3_sequence %in% sequences_to_remove
  )


# ==============================================================================
# SUMMARY OF DIAMETER_1 AND D_130
# ==============================================================================
library(dplyr)

# Summary of missing values in diameter_1 and diameter_130 for filtered data
attr_summary2 <- trees_consistent_L4 %>%
  # Apply data quality filters
  filter(
    discard_nodead == FALSE,
    statuserror_record == FALSE,
    errorlife_record == FALSE,
    nalife_record == FALSE
  ) %>%
  # Group by site and calculate count of NA values
  group_by(composed_site_id) %>%
  summarise(
    diam_unify_NA = sum(is.na(diam_unify)),
    diameter_1_NA = sum(is.na(diameter_1)),
    .groups = "drop"
  )


# ==============================================================================
# HARMONIZING DIAMETER_1 (DIAM_130)
# ==============================================================================

trees_consistent_L4 <- trees_consistent_L4 %>%
  # 1. Sort by tree ID and inventory year to ensure chronological order
  arrange(tree_unique_id, inventory_year) %>%
  
  # 2. Group by tree_unique_id to perform chronological fill for VUK
  group_by(tree_unique_id) %>%
  mutate(
    # Create a temporary column for VUK where NA is filled from previous years
    diam_vuk_filled = ifelse(institute == "VUK", diameter_1, NA_real_)
  ) %>%
  tidyr::fill(diam_vuk_filled, .direction = "down") %>%
  
  # 3. Apply final unification rules based on institute
  mutate(
    diam_unify = case_when(
      # --- LOGIC FOR VUK ---
      # A) VUK: Keep diameter_1 if present, or take from previous years (from fill)
      institute == "VUK" & !is.na(diam_vuk_filled) ~ diam_vuk_filled,
      
      # B) VUK: If still NA (no previous record), use diameter_130 ONLY for piece_id == "a"
      institute == "VUK" & is.na(diam_vuk_filled) & piece_id == "a" ~ diameter_130,
      
      # --- LOGIC FOR OTHER INSTITUTES ---
      # C) Others: Take diameter_1 if present
      institute != "VUK" & !is.na(diameter_1) ~ diameter_1,
      
      # D) Others: If diameter_1 is NA, take diameter_130 from the SAME ROW (regardless of history or piece_id)
      institute != "VUK" & is.na(diameter_1) ~ diameter_130,
      
      # Fallback for unhandled cases
      TRUE ~ NA_real_
    )
  ) %>%
  
  # Clean up temporary column and remove grouping
  select(-diam_vuk_filled) %>%
  ungroup()



library(dplyr)
library(tidyr)

# ==============================================================================
# CREATE DIAMETER_FIRST COLUMN AND BROADCAST ACROSS ALL ROWS OF EACH TREE
# ==============================================================================
trees_consistent_L4 <- trees_consistent_L4 %>%
  # 1. Sort chronologically by tree and inventory year
  arrange(tree_unique_id, inventory_year) %>%
  
  # 2. Group by tree_unique_id to evaluate each tree individually
  group_by(tree_unique_id) %>%
  
  # 3. Calculate diameter_first directly within each group
  mutate(
    diam_first = {
      # Find row indices where status3 is "D"
      dead_rows <- which(status3 == "D")
      
      if (length(dead_rows) > 0) {
        # If the tree died, take diam_unify from the FIRST year of death
        first_dead_idx <- dead_rows[1]
        diam_unify[first_dead_idx]
      } else {
        # If the tree NEVER died, return NA for all its rows
        NA_real_
      }
    }
  ) %>%
  
  # 4. Remove grouping
  ungroup()


# ADD first_NA COLUMN (TRUE = is NA, FALSE = has number)
trees_consistent_L4 <- trees_consistent_L4 %>%
  mutate(
    diamunify_NA = is.na(diam_unify)
  )

trees_consistent_L4 <- trees_consistent_L4 %>%
  mutate(
    diamfirst_NA = is.na(diam_first)
  )
str(trees_consistent_L4)


# ==============================================================================
# DIAMETER THRESHOLD MAPPING
# ==============================================================================

# Clean the design dataset, ensuring thresholds exist, and convert cm to mm.
# Include circle_no in the selection for strict spatial partitioning.
design_clean <- all_desgin %>%
  filter(!is.na(lying_alive_threshold), !is.na(lying_dead_threshold)) %>%
  mutate(
    a_threshold_mm = lying_alive_threshold * 10,
    d_threshold_mm  = lying_dead_threshold * 10
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

# Safe left_join function that unifies plot_id data types to character
safe_left_join_design <- function(x, y) {
  x %>%
    mutate(plot_id = as.character(plot_id)) %>%
    left_join(
      y %>% mutate(plot_id = as.character(plot_id)),
      by = c("composed_site_id", "inventory_year", "inventory_id", "inventory_type", "circle_no", "plot_id")
    )
}

trees_consistent_L4 <- safe_left_join_design(trees_consistent_L4, design_long)

library(dplyr)


# ==============================================================================
# EVALUATION OF DBH AGAINST THRESHOLDS
# ==============================================================================

trees_consistent_L4 <- trees_consistent_L4 %>%
  mutate(
 
    
    # 2. Check DBH thresholds for DEAD records (life == "D")
    d_thresh_check = case_when(
      status3 == "D" & diam_first >= d_threshold_mm ~ "within",
      status3 == "D" & diam_first <  d_threshold_mm ~ "under",
      TRUE ~ NA_character_  # Returns NULL (NA) for life == "A" or any other cases
    )
  )


# SOLUTION: The dataset also retains stems with diameters lower than threshold.
# Filtering to the desired diameter ranges will be done later.



# ==============================================================================
# JOING SITE-LEVEL VARIABLES TO LYING STEMS DATASET
# ==============================================================================

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




# ==============================================================================
# ==============================================================================
# CREATING BIRTH AND DEATH
# ==============================================================================
# ==============================================================================




# CONSTRUCTING PLOT TIMELINE FROM TREE-LEVEL DATASET
# CALCULATING CENSUS RANK AND BIRTH/DEATH (WITH 0 FOR UNKNOWN)
# ==============================================================================
# In trees_3inv table adjusted plot_id2 is necessary

# The following plot_time line is outdated because in some LWF sites
# 2 new inventories were added (2025, 2026)
# ------------------------------------------------------------------------------

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
    
  
  # 5. Clean up grouping and apply final sorting
  ungroup() %>%
  arrange(institute, composed_site_id, plot_id2, inventory_year)) 


# In the next step, table was adjusted in table outside of R - for two overlapping plots
# in VUK Ranspurk and Zofin
# CALCULATE BIRTH AND DEATH MIDPOINTS FOR BIG AND SMALL PLOTS IN PLOT_TIMELINE
# Applies for Ranspurk and Zofin
# ==============================================================================

plot_timeline_updated <- plot_timeline_updated %>%
  mutate(
    # --- 1. BIG PLOTS (Four-census variant) ---
    # Birth midpoint: average of current year and the previous valid census for BIG
    birth_big = if_else(
      is.na(inventory_before_big) | inventory_before_big == "", 
      0, 
      (inventory_year_big + as.numeric(inventory_before_big)) / 2
    ),
    
    # Death midpoint: average of current year and the next valid census for BIG
    death_big = if_else(
      is.na(inventory_next_big) | inventory_next_big == "", 
      0, 
      (inventory_year_big + as.numeric(inventory_next_big)) / 2
    ),
    
    # --- 2. SMALL PLOTS (Six-census variant) ---
    # Birth midpoint: average of current year and the previous valid census for SMALL
    birth_small = if_else(
      is.na(inventory_before_small) | inventory_before_small == "", 
      0, 
      (inventory_year_small + as.numeric(inventory_before_small)) / 2
    ),
    
    # Death midpoint: average of current year and the next valid census for SMALL
    death_small = if_else(
      is.na(inventory_next_small) | inventory_next_small == "", 
      0, 
      (inventory_year_small + as.numeric(inventory_next_small)) / 2
    )
  )

# ADJUSTING BIRTH AND DEATH
# ==============================================================================
# Step 1: Create the 'big_small' column on the tree level
# ------------------------------------------------------------------------------

trees_consistent_L5 <- trees_consistent_L4 %>%
  group_by(tree_unique_id) %>%
  mutate(
    # A. Site and year tracking
    is_ranspurk      = any(grepl("Ranspurk", tree_unique_id)),
    has_ranspurk_yrs = any(inventory_year %in% c(2015, 2025)),
    
    is_zofin         = any(grepl("Zofin", tree_unique_id)),
    has_zofin_yrs    = any(inventory_year %in% c(2012, 2017)),
    
    # B. Classification logic (without DBH threshold rules)
    big_small = case_when(
      # Ranšpurk: if recorded in 2015 or 2025 -> "small", otherwise "big"
      is_ranspurk & has_ranspurk_yrs  ~ "small",
      is_ranspurk & !has_ranspurk_yrs ~ "big",
      
      # Žofín: if recorded in 2012 or 2017 -> "small", otherwise "big"
      is_zofin & has_zofin_yrs        ~ "small",
      is_zofin & !has_zofin_yrs       ~ "big",
      
      # Default for all other sites -> "small"
      TRUE                            ~ "small"
    )
  ) %>%
  ungroup() %>%
  # Clean up temporary logical flags
  select(-is_ranspurk, -has_ranspurk_yrs, -is_zofin, -has_zofin_yrs)

library(tidyverse)

# ==============================================================================
# STEP 2: PREPARE SINGLE LONG TIMELINE LOOKUP FROM timeline_L
# ==============================================================================
# Reshape timeline_L into a long format containing: 
# composed_site_id, big_small, inventory_year, birth, and death

timeline_lookup <- timeline_L %>%
  pivot_longer(
    cols = c(ends_with("_big"), ends_with("_small")),
    names_to = c(".value", "big_small"),
    names_pattern = "(.*)_(big|small)"
  ) %>%
  filter(!is.na(inventory_year) & inventory_year != "") %>%
  mutate(inventory_year = as.numeric(inventory_year)) %>%
  select(composed_site_id, big_small, inventory_year, birth, death) %>%
  distinct(composed_site_id, big_small, inventory_year, .keep_all = TRUE)


# ==============================================================================
# STEP 3 & 4: CALCULATE DEAD STEM YEARS, JOIN TIMELINE & UPDATE trees_consistent_L5
# ==============================================================================

trees_consistent_L5 <- trees_consistent_L5 %>%
  # 1. Safely calculate first and last year with status3 == "D" at the tree level
  group_by(tree_unique_id) %>%
  mutate(
    first_d_year = if_else(
      any(status3 == "D", na.rm = TRUE),
      min(inventory_year[status3 == "D"], na.rm = TRUE),
      NA_real_
    ),
    last_d_year = if_else(
      any(status3 == "D", na.rm = TRUE),
      max(inventory_year[status3 == "D"], na.rm = TRUE),
      NA_real_
    )
  ) %>%
  ungroup() %>%
  
  # 2. Join BIRTH based on site, plot size group (BIG/SMALL), and first_d_year
  left_join(
    timeline_lookup %>% select(composed_site_id, big_small, inventory_year, birth),
    by = c("composed_site_id", "big_small", "first_d_year" = "inventory_year")
  ) %>%
  
  # 3. Join DEATH based on site, plot size group (BIG/SMALL), and last_d_year
  left_join(
    timeline_lookup %>% select(composed_site_id, big_small, inventory_year, death),
    by = c("composed_site_id", "big_small", "last_d_year" = "inventory_year")
  ) %>%
  
  # 4. Sort chronologically
  arrange(tree_unique_id, inventory_year)




# ==============================================================================
# IDENTIFY SPECIES STABILITY AND CLEAN UNKNOWN TAXONS PER STEM
# ==============================================================================

trees_consistent_L5 <- trees_consistent_L5 %>%
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


# ==============================================================================
# GROUPING TREE SPECIES, A GROUP MIN. OF 1000 INDIVIDUALS
# ==============================================================================

library(dplyr)

trees_consistent_L5 <- trees_consistent_L5 %>%
  mutate(
    species_group = case_when(
      tree_species %in% c("Abies alba", "Abies concolor", "Coniferous", "Pinus sylvestris") ~ "Abies alba",
      tree_species %in% c("Acer campestre", "Acer platanoides", "Acer pseudoplatanus") ~ "Acer",
      tree_species == "Fagus sylvatica" ~ "Fagus",
      tree_species %in% c("Fraxinus angustifolia", "Fraxinus excelsior") ~ "Fraxinus",
      tree_species == "Picea abies" ~ "Picea",
      tree_species %in% c("Quercus petraea", "Quercus robur") ~ "Quercus",
      tree_species %in% c("Ulmus", "Ulmus carpinifolia", "Ulmus glabra", "Ulmus laevis", "Ulmus minor") ~ "Ulmus",
      tree_species %in% c("Aesculus hippocastanum", "Carpinus betulus", "Juglans nigra", "Malus sylvestris", 
                          "Prunus avium", "Prunus spinosa", "Pyrus communis", "Pyrus pyraster", 
                          "Robinia pseudoacacia", "Sorbus aria", "Sorbus aucuparia", "Sorbus torminalis") ~ "broad_hard",
      tree_species %in% c("Alnus glutinosa", "Betula", "Betula pendula", "Betula spp.", 
                          "Corylus avellana", "Crataegus", "Crataegus monogyna", 
                          "Populus", "Populus alba", "Populus tremula", "Salix caprea", 
                          "Salix sp.", "Salix spp.", "Sambucus nigra", "Tilia cordata", "Tilia spp.") ~ "broad_soft",
      TRUE ~ "unknown"
    )
  )
