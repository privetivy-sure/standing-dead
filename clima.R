# ==============================================================================
# CLIMATIC VARIABLES
# ==============================================================================


# UPLOAD DESIGN, JOIN GEOM TO SWD ANALYZED SITES
# ==============================================================================
# Load required libraries
library(terra)
library(dplyr)
library(sf)


# FIX THE POSTGRES/PROJ CONFLICT
# Forces R to use its own spatial libraries instead of PostgreSQL local path
Sys.setenv("PROJ_LIB" = system.file("proj", package = "terra")[1])



# =========================================================================
# STEP A: JOIN AND REPAIR GEOMETRIES
# =========================================================================

# 3. JOIN: Bring the 'geom' column from 'design_all' into 'sites_cwd'
sites_cwd_joined <- sites_cwd %>%
  left_join(design_all %>% select(composed_site_id, geom), by = "composed_site_id")

# 4. PARSE THE GEOMETRY COLUMN (WKB hex string to SF object)
wkb_geometries <- structure(as.list(sites_cwd_joined$geom), class = "WKB")
geom_sfc <- st_as_sfc(wkb_geometries, EWKB = TRUE)

# Bind back to data and set coordinate system to WGS84 (EPSG:4326)
sites_sf <- st_sf(sites_cwd_joined %>% select(-geom), geometry = geom_sfc)
st_crs(sites_sf) <- 4326


# =========================================================================
# STEP B: CALCULATE CENTROIDS
# =========================================================================

# 5. Temporarily disable strict S2 geometry engine to allow GEOS tolerance
sf_use_s2(FALSE)

# Repair self-intersections or duplicate vertices (errors like Loop 0 is not valid)
sites_sf_repaired <- st_make_valid(sites_sf)

# Calculate centroids (works on points, polygons and geometry collections)
sites_centroids_sf <- suppressWarnings(st_centroid(sites_sf_repaired))

# Re-enable the S2 geometry engine
sf_use_s2(TRUE)

# 6. CONVERT TO 'TERRA' SPATVECTOR FOR RASTER EXTRACTION
# This allows the centroids to interact natively with WorldClim rasters
sites_centroids <- vect(sites_centroids_sf)



# =========================================================================
# SYSTEM SETUP & ENVIRONMENT
# =========================================================================

# Force R to use its own spatial libraries instead of PostgreSQL local path
Sys.setenv("PROJ_LIB" = system.file("proj", package = "terra")[1])

library(terra)
library(dplyr)
library(geodata)

# Define your local climate directory
climate_dir <- "D:/DEADWOOD_D/WildCard/Residence/standing-dead-git/temp_data/climate/geodat"
if(!dir.exists(climate_dir)) dir.create(climate_dir, recursive = TRUE)

# Define your exact European study area bounding box (Belgium to Carpathians, Bosnia to Poland)
study_extent <- ext(2.0, 26.5, 42.0, 55.5)


# =========================================================================
# PART 1: WORLDCLIM AUTOMATIC WEB DOWNLOAD & CROP (Runs only if files are missing)
# =========================================================================

# We download only the specific 30x30 degree WorldClim tile covering Central Europe.

tavg_path <- file.path(climate_dir, "cropped_wc2.1_30s_tavg.tif")
bio_path  <- file.path(climate_dir, "cropped_wc2.1_30s_bio.tif")
vapr_path <- file.path(climate_dir, "cropped_wc2.1_30s_vapr.tif")

if (!file.exists(tavg_path) || !file.exists(bio_path) || !file.exists(vapr_path)) {
  message("--- PART 1: Missing local files detected. Downloading from WorldClim web ---")
  
  # Coordinates of a point in the center of your study area (e.g., Czech Republic)
  # This ensures we download the exact single tile we need (Tile spanning Lon 0-30, Lat 30-60)
  target_lon <- 15.0
  target_lat <- 50.0
  
  # 1. Download only the specific tile to R temporary directory
  message("--> Downloading VAPR tile covering Central Europe...")
  vapr_tile <- geodata::worldclim_tile(var = "vapr", lon = target_lon, lat = target_lat, path = tempdir())
  
  message("--> Downloading TAVG tile covering Central Europe...")
  tavg_tile <- geodata::worldclim_tile(var = "tavg", lon = target_lon, lat = target_lat, path = tempdir())
  
  message("--> Downloading BIO tile covering Central Europe...")
  bio_tile  <- geodata::worldclim_tile(var = "bio", lon = target_lon, lat = target_lat, path = tempdir())
  
  
  # 2. Crop the tile to your exact study extent (excluding unnecessary areas like Mediterranean/North Africa)
  message("--> Cropping tile layers to exact European study area...")
  vapr_cropped <- crop(vapr_tile, study_extent)
  tavg_cropped <- crop(tavg_tile, study_extent)
  bio_cropped  <- crop(bio_tile, study_extent)
  
  # 3. Save the cropped versions locally to your disk
  message("--> Saving cropped layers to disk...")
  writeRaster(vapr_cropped, filename = vapr_path, overwrite = TRUE)
  writeRaster(tavg_cropped, filename = tavg_path, overwrite = TRUE)
  writeRaster(bio_cropped,  filename = bio_path,  overwrite = TRUE)
  
  message("--- PART 1: Tile download and cropping completed successfully! ---")
} else {
  message("--- PART 1: All local cropped files already exist on disk. Skipping download. ---")
}


# =========================================================================
# PART 2: LOAD LOCAL CROPPED RASTERS FROM DISK
# =========================================================================
# Fast loading of local lightweight files (no internet connection required).

message("--- PART 2: Loading cropped rasters from PC ---")

tavg_raster <- rast(tavg_path)
bio_raster  <- rast(bio_path)
vapr_raster <- rast(vapr_path)


# =========================================================================
# PART 3: CALCULATE VPD FROM LOCAL RASTERS
# =========================================================================
# Saturated Vapor Pressure (Es) = 0.6108 * exp((17.27 * Temp) / (Temp + 237.3))
# Vapor Pressure Deficit (VPD)  = Es - Ea (Actual Vapor Pressure)

vpd_months <- list()

for (m in 1:12) {
  # Get monthly temperature
  temp_m <- tavg_raster[[m]]
  
  # Check and apply temperature scaling factor correction if needed (Temp / 10)
  if (global(temp_m, "mean", na.rm = TRUE)[1, 1] > 50) {
    temp_m <- temp_m / 10
  }
  
  # Get actual vapor pressure (Ea)
  ea_m <- vapr_raster[[m]]
  
  # 1. Saturated vapor pressure (Es) in kPa using Tetens equation
  es_m <- 0.6108 * exp((17.27 * temp_m) / (temp_m + 237.3))
  
  # 2. Calculate VPD = Es - Ea
  vpd_m <- es_m - ea_m
  # Clamp negative values to absolute physical zero (using correct lowercase 'rcl')
  vpd_m <- classify(vpd_m, rcl = matrix(c(-Inf, 0, 0), ncol = 3, byrow = TRUE))
  
  vpd_months[[m]] <- vpd_m
}

# Combine all 12 monthly layers into a single SpatRaster stack
vpd_stack <- rast(vpd_months)


# =========================================================================
# CALCULATE SEASONAL METRICS & APPLY QGIS COMPATIBILITY FIXES
# =========================================================================

# --- 1. GROWING SEASON VPD (Months 4 to 9: April to September) ---
# Extract layers 4 to 9 and calculate the average
vpd_growing_season <- mean(vpd_stack[[4:9]])

# QGIS Fix: Define Coordinate Reference System (CRS) and clean band name
crs(vpd_growing_season) <- crs(tavg_raster)
names(vpd_growing_season) <- "vpd_growing_season_kPa"

# Save to disk using LZW compression and generate a .tfw world file
writeRaster(vpd_growing_season, 
            filename = file.path(climate_dir, "europe_vpd_growing_season_1km.tif"), 
            overwrite = TRUE,
            gdal = c("COMPRESS=LZW", "TFW=YES"))


# --- 2. SUMMER VPD (Months 6 to 8: June to August) ---
# Extract layers 6 to 8 and calculate the average
vpd_summer <- mean(vpd_stack[[6:8]])

# QGIS Fix: Define Coordinate Reference System (CRS) and clean band name
crs(vpd_summer) <- crs(tavg_raster)
names(vpd_summer) <- "vpd_summer_kPa"

# Save to disk using LZW compression and generate a .tfw world file
writeRaster(vpd_summer, 
            filename = file.path(climate_dir, "europe_vpd_summer_1km.tif"), 
            overwrite = TRUE,
            gdal = c("COMPRESS=LZW", "TFW=YES"))


# =========================================================================
# PART 4: EXTRACT AND CALCULATE ECO-CLIMATOLOGICAL INDEXES
# =========================================================================
# Extract values into your sites_clima table using X (longitude) and Y (latitude) columns.

message("--- PART 4: Extracting values and calculating final indexes ---")

# Prepare coordinate matrix (First column = X / Longitude, Second column = Y / Latitude)
coords_matrix <- as.matrix(sites_clima[, c("centroid_x", "centroid_y")])

# 1. Extract raw climate values using coordinates
# Extract all 12 monthly temperature bands at once to optimize performance
raw_tavg_all  <- terra::extract(tavg_raster, coords_matrix) 

# Extract bioclimatic variables (bio)
raw_temp_ann  <- terra::extract(bio_raster[[1]], coords_matrix)[, 1]    # BIO1 (Annual Mean Temp)
raw_temp_seas <- terra::extract(bio_raster[[4]], coords_matrix)[, 1]    # BIO4 (Temp Seasonality)
raw_temp_rang <- terra::extract(bio_raster[[7]], coords_matrix)[, 1]    # BIO7 (Temp Annual Range)
raw_prec_ann  <- terra::extract(bio_raster[[12]], coords_matrix)[, 1]   # BIO12 (Annual Precipitation)
raw_prec_seas <- terra::extract(bio_raster[[15]], coords_matrix)[, 1]   # BIO15 (Precipitation Seasonality)

# Extract newly created seasonal VPDs
raw_vpd_gs    <- terra::extract(vpd_growing_season, coords_matrix)[, 1]  # Growing season VPD

# 2. Extract January and July temperatures from our 12-month temperature matrix
raw_temp_jan  <- raw_tavg_all[, 1]   # Band 1: January
raw_temp_jul  <- raw_tavg_all[, 7]   # Band 7: July


# 3. Handle Scaling Factor for Temperatures (divide by 10 if values are > 50)
if (mean(raw_temp_jul, na.rm = TRUE) > 50) {
  # Apply scaling divisor to monthly temperatures matrix
  tavg_scaled <- raw_tavg_all / 10
  
  # Save scaled single variables to sites_clima
  sites_clima$temp_january_1km   <- raw_temp_jan / 10
  sites_clima$temp_july_1km      <- raw_temp_jul / 10
  sites_clima$temp_annual_1km    <- raw_temp_ann / 10
  sites_clima$temp_seasonality   <- raw_temp_seas / 100  # Convert SD*100 to actual SD in °C
  sites_clima$temp_annual_range  <- raw_temp_rang / 10
} else {
  tavg_scaled <- raw_tavg_all
  
  sites_clima$temp_january_1km   <- raw_temp_jan
  sites_clima$temp_july_1km      <- raw_temp_jul
  sites_clima$temp_annual_1km    <- raw_temp_ann
  sites_clima$temp_seasonality   <- raw_temp_seas / 100
  sites_clima$temp_annual_range  <- raw_temp_rang
}

# Assign precipitation and VPD directly (no scaling divisor needed)
sites_clima$prec_annual_1km      <- raw_prec_ann
sites_clima$prec_seasonality     <- raw_prec_seas
sites_clima$vpd_season_1km       <- raw_vpd_gs


# 4. Calculate Annual Temperature Sum above 5°C (GDD5)
# Define the number of days for each month (January to December)
days_in_months <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Create a temporary matrix to store GDD values for each month
gdd_matrix <- matrix(0, nrow = nrow(tavg_scaled), ncol = 12)

for (m in 1:12) {
  temp_m <- tavg_scaled[, m]
  # If monthly temperature is > 5°C, calculate: (Temp - 5) * days in month
  gdd_matrix[, m] <- ifelse(!is.na(temp_m) & temp_m > 5, (temp_m - 5) * days_in_months[m], 0)
}

# Sum all monthly GDD values to get the annual total
sites_clima$gdd_above_5 <- rowSums(gdd_matrix, na.rm = TRUE)


# 5. Calculate De Martonne Aridity Index: Precipitation / (Temperature + 10)
sites_clima$de_martonne_aridity  <- sites_clima$prec_annual_1km / (sites_clima$temp_annual_1km + 10)

message("--- ALL STEPS COMPLETED ---")
head(sites_clima)


# ==============================================================================
# VISUALISATION OF CLIMATIC VARIABLES
# ==============================================================================

library(ggplot2)
library(plotly)
library(scales)


# -------------------------------------------------------------------------
# 2D VISUALIZATION: Energy vs. Vapor Pressure Deficit
# -------------------------------------------------------------------------
# - It's linear relationship; doesn't make sense for explaining decomposition

p2d <- ggplot(sites_clima, aes(x = gdd_above_5, 
                               y = vpd_season, 
                               color = as.factor(his_fortype))) +
  
  geom_point(size = 3, alpha = 0.8) +
  
  scale_color_viridis_d(name = "Forest Type", option = "turbo") +
  
  labs(
    title = "Climatic Differentiation of Study Sites",
    subtitle = "Tradeoff between thermal energy (GDD5) and atmospheric drying stress (VPD)",
    x = "Growing Degree Days Above 5°C (GDD5)",
    y = "Growing Season VPD (April-October, kPa)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# Zobrazení grafu
print(p2d)

# Display the 2D plot
print(p2d)


# ------------------------------------------------------------------------------
# 2D CONTINENTALITY - PRECIPITATION GRADIENT
# ------------------------------------------------------------------------------


library(ggplot2)

p2d_clean <- ggplot(sites_clima, aes(x = temp_seasonality, 
                                     y = temp_annual, 
                                     color = as.factor(his_fortype), 
                                     )) +
  
  geom_point(alpha = 0.75) +
  
  
  scale_color_viridis_d(name = "Forest Type", option = "turbo") +
  
  
  labs(
    title = "Climate Space: Continentality vs. Mean Annual Temp.",
    subtitle = "Color shows forest type",
    x = "Temperature Seasonality (Standard Deviation * 100)",
    y = "Mean Annual Temp. (°C)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p2d_clean)

# -------------------------------------------------------------------------
# MATRIX OF X Y GRAPHS
# ------------------------------------------------------------------------------

library(GGally)

# Select only the specific variables you want to evaluate side-by-side
climate_sub <- sites_clima[, c("temp_seasonality", "prec_annual", 
                               "prec_seasonality", "de_martonne_aridity", "vpd_season")]

# Generate the pairs plot matrix
ggpairs(
  climate_sub,
  lower = list(continuous = wrap("points", alpha = 0.6, color = "#2c3e50")),
  upper = list(continuous = wrap("cor", size = 4.5, color = "darkred")),
  diag = list(continuous = wrap("densityDiag", fill = "#34495e", alpha = 0.5))
) +
  theme_bw()
