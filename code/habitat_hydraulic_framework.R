# =============================================================================
# Land Cover × Hydraulics × Depth-to-Groundwater Analysis
# Upper Mississippi River — Pools 11–19
#
# STRUCTURE:
#   PART 0: Settings and libraries
#   PART 1: Load and cache hydraulic rasters
#   PART 2: Load and cache land cover data (LTRM, NLCD, LANDFIRE)
#   PART 3: Load and cache depth-to-groundwater
#   PART 4: Build analysis datasets (sample rasters, extract to polygons)
#   PART 5: Inundation frequency analysis
#   PART 6: DTW analysis
#   PART 7: Combined inundation + DTW models
#   PART 8: Cross-classification of land cover systems
#   PART 9: Plots
#   PART 10: Export
#
# CACHING STRATEGY:
#   Expensive spatial operations write .rds or .tif files to cache/.
#   Each step checks for the cache file first and skips recomputation.
#   To force recomputation, delete the relevant cache file.
# =============================================================================

# =============================================================================
# PART 0: Settings and libraries
# =============================================================================

library(tidyverse)
library(sf)
library(terra)
library(ggdist)
library(patchwork)
library(dunn.test)
library(FedData)
library(rlandfire)

# ── User settings ─────────────────────────────────────────────────────────────

# Land cover class column in LTRM LULC shapefiles
lc_col <- "CLASS_15_N"

# Pools to load
pools <- 11:19

# Raster aggregation factor for extraction (2m → 20m = 100x fewer cells)
agg_factor <- 10

# Minimum fraction of polygon covered by valid hydraulic data
min_coverage <- 0.10

# Number of cells to sample from raster stacks
n_sample <- 2000000#200000

# Cache directory — all saved spatial intermediates go here
cache_dir <- "cache"
dir.create(cache_dir, showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

# Helper: load from cache if it exists, otherwise run expr and cache result
cache_rds <- function(path, expr) {
  if (file.exists(path)) {
    cat("  [cache] Loading", basename(path), "\n")
    return(readRDS(path))
  }
  result <- expr
  saveRDS(result, path)
  cat("  [saved]", basename(path), "\n")
  result
}

cache_tif <- function(path, expr, ...) {
  if (file.exists(path)) {
    cat("  [cache] Loading", basename(path), "\n")
    return(rast(path))
  }
  result <- expr
  writeRaster(result, path, overwrite = TRUE, ...)
  cat("  [saved]", basename(path), "\n")
  result
}

# Inundation frequency level order — used throughout
inund_levels <- c(
  "Permanent / baseflow",
  "5-yr event",
  "10-yr event",
  "25-yr event",
  "Not inundated (>25-yr)"
)

# Inundation color palette
inund_colors <- c(
  "Permanent / baseflow"   = "#1B6B3A",
  "5-yr event"             = "#2A7D6F",
  "10-yr event"            = "#D4850A",
  "25-yr event"            = "#C0392B",
  "Not inundated (>25-yr)" = "#B0BAC5"
)

scenario_colors <- c(
  "Baseflow (60k cfs)" = "#2A7D6F",
  "5-yr (167k cfs)"    = "#D4850A",
  "10-yr (190k cfs)"   = "#7B3F9E",
  "25-yr (241k cfs)"   = "#C0392B"
)

# =============================================================================
# PART 1: Load and cache hydraulic rasters
# =============================================================================
# Rasters are aggregated 2m → ~20m and kept in native CRS.
# Vectors are reprojected to raster CRS rather than reprojecting large rasters.
# Aggregated stacks are cached as .tif files.

bf_dir   <- "data/hydraulic/PhaseII_Baseflow_60kcfs"
fyr_dir  <- "data/hydraulic/PhaseII_5YR_167kCFS"
tyr_dir  <- "data/hydraulic/PhaseII_10Y_190kcfs"
tyfr_dir <- "data/hydraulic/PhaseII_25Y_241kcfs"

rast_paths <- list(
  bf = list(
    depth    = file.path(bf_dir,
                         "Depth_60k.Mississippi_UMR_Ph2_2m_Export_Terrain.Mississippi_UMR_Ph2_2m_Export_Terrain.tif"),
    velocity = file.path(bf_dir,
                         "Velocity.Mississippi_UMR_Ph2_2m_Export_Terrain.Mississippi_UMR_Ph2_2m_Export_Terrain.tif"),
    froude   = file.path(bf_dir,
                         "Froude.Mississippi_UMR_Ph2_2m_Export_Terrain.Mississippi_UMR_Ph2_2m_Export_Terrain.tif")
  ),
  fyr = list(
    depth    = file.path(fyr_dir,
                         "Depth_5Y.Mississippi_UMR_Ph2_2m_Export_Terrain.Mississippi_UMR_Ph2_2m_Export_Terrain.tif"),
    velocity = file.path(fyr_dir,
                         "Velocity_5Y.Mississippi_UMR_Ph2_2m_Export_Terrain.Mississippi_UMR_Ph2_2m_Export_Terrain.tif"),
    froude   = file.path(fyr_dir,
                         "Froude_5Y.Mississippi_UMR_Ph2_2m_Export_Terrain.Mississippi_UMR_Ph2_2m_Export_Terrain.tif")
  ),
  tyr = list(
    depth    = file.path(tyr_dir,
                         "Depth_10Y.Mississippi_UMR_Ph2_2m_Export_Terrain.Mississippi_UMR_Ph2_2m_Export_Terrain.tif"),
    velocity = file.path(tyr_dir,
                         "Velocity_10Y.Mississippi_UMR_Ph2_2m_Export_Terrain.Mississippi_UMR_Ph2_2m_Export_Terrain.tif"),
    froude   = file.path(tyr_dir,
                         "Froude_10Y.Mississippi_UMR_Ph2_2m_Export_Terrain.Mississippi_UMR_Ph2_2m_Export_Terrain.tif")
  ),
  tyfr = list(
    depth    = file.path(tyfr_dir,
                         "Depth_25Y.Mississippi_UMR_Ph2_2m_Export_Terrain.Mississippi_UMR_Ph2_2m_Export_Terrain.tif"),
    velocity = file.path(tyfr_dir,
                         "Velocity_25Y.Mississippi_UMR_Ph2_2m_Export_Terrain.Mississippi_UMR_Ph2_2m_Export_Terrain.tif"),
    froude   = file.path(tyfr_dir,
                         "Froude_25Y.Mississippi_UMR_Ph2_2m_Export_Terrain.Mississippi_UMR_Ph2_2m_Export_Terrain.tif")
  )
)

scenario_labels <- c(
  bf   = "Baseflow (60k cfs)",
  fyr  = "5-yr (167k cfs)",
  tyr  = "10-yr (190k cfs)",
  tyfr = "25-yr (241k cfs)"
)

load_and_agg <- function(path, label, fact = agg_factor) {
  cat("  Aggregating", label, "...\n")
  r     <- rast(path)
  r_agg <- aggregate(r, fact = fact, fun = "mean", na.rm = TRUE)
  rm(r); gc()
  r_agg
}

load_stack <- function(scenario_key, cache_path) {
  cache_tif(cache_path, {
    s <- c(
      load_and_agg(rast_paths[[scenario_key]]$depth,
                   paste(scenario_key, "depth")),
      load_and_agg(rast_paths[[scenario_key]]$velocity,
                   paste(scenario_key, "velocity")),
      load_and_agg(rast_paths[[scenario_key]]$froude,
                   paste(scenario_key, "froude"))
    )
    names(s) <- c("depth", "velocity", "froude")
    s
  }, gdal = c("COMPRESS=LZW", "BIGTIFF=IF_NEEDED"))
}

cat("Loading hydraulic stacks (cached after first run)...\n")
bf_stack   <- load_stack("bf",   file.path(cache_dir, "bf_stack_agg.tif"))
fyr_stack  <- load_stack("fyr",  file.path(cache_dir, "fyr_stack_agg.tif"))
tyr_stack  <- load_stack("tyr",  file.path(cache_dir, "tyr_stack_agg.tif"))
tyfr_stack <- load_stack("tyfr", file.path(cache_dir, "tyfr_stack_agg.tif"))

# Convenience references to individual depth layers
depth_bf   <- bf_stack[["depth"]]
depth_fyr  <- fyr_stack[["depth"]]
depth_tyr  <- tyr_stack[["depth"]]
depth_tyfr <- tyfr_stack[["depth"]]

# Raster CRS — all stacks share the same CRS
rast_crs <- crs(bf_stack)
cat("Raster CRS:", crs(bf_stack, describe = TRUE)$name, "\n")

# Hydraulic model footprint: 1 where ANY scenario has valid depth
footprint <- cache_tif(
  file.path(cache_dir, "hydraulic_footprint.tif"),
  ifel(!is.na(depth_bf) | !is.na(depth_fyr) |
         !is.na(depth_tyr) | !is.na(depth_tyfr), 1, NA)
)
cat("Footprint non-NA cells:", global(footprint, "notNA")[[1]], "\n")

# =============================================================================
# PART 2: Load and cache land cover data
# =============================================================================

# ── 2a. LTRM LULC polygons ────────────────────────────────────────────────────

lulc_clip <- cache_rds(file.path(cache_dir, "lulc_clip.rds"), {
  cat("Loading LTRM LULC shapefiles...\n")
  
  lulc_raw <- map(pools, function(p) {
    pool_str <- sprintf("%02d", p)
    path     <- sprintf("data/lulc/2020_LCU_Pool_%s/lc_2020_p%s.shp",
                        pool_str, pool_str)
    if (!file.exists(path)) { cat("  Not found:", path, "\n"); return(NULL) }
    cat("  Pool", p, "...\n")
    st_read(path, quiet = TRUE) %>%
      st_transform(crs = rast_crs) %>%
      mutate(pool = p)
  }) %>% compact() %>% bind_rows()
  
  # Clip to raster extent by centroid coordinates
  xy <- st_coordinates(st_centroid(lulc_raw$geometry))
  lulc_raw %>%
    mutate(cx = xy[, 1], cy = xy[, 2]) %>%
    filter(cx >= ext(bf_stack)$xmin, cx <= ext(bf_stack)$xmax,
           cy >= ext(bf_stack)$ymin, cy <= ext(bf_stack)$ymax) %>%
    select(-cx, -cy) %>%
    mutate(poly_id  = row_number(),
           lc_class = .data[[lc_col]])
})

cat("LTRM polygons:", nrow(lulc_clip), "\n")
cat("Classes:\n"); print(table(lulc_clip$lc_class))

# ── 2b. NLCD ──────────────────────────────────────────────────────────────────

nlcd_masked <- cache_tif(
  file.path(cache_dir, "nlcd_masked.tif"),
  {
    cat("Downloading NLCD (may take a few minutes)...\n")
    nlcd_raw <- get_nlcd(template = tyfr_stack, label = "hydraulic_nlcd",
                         year = 2021)
    
    fp_nlcd <- project(footprint, crs(nlcd_raw), method = "near")
    fp_nlcd <- resample(fp_nlcd, nlcd_raw, method = "near")
    mask(nlcd_raw, fp_nlcd)
  },
  datatype = "INT1U",
  gdal = c("COMPRESS=LZW", "BIGTIFF=IF_NEEDED")
)
cat("NLCD masked cells:", global(nlcd_masked, "notNA")[[1]], "\n")

# NLCD label lookup
nlcd_labels <- tribble(
  ~nlcd_code, ~nlcd_label,
  11, "Open Water",          21, "Developed, Open Space",
  22, "Developed, Low",      23, "Developed, Medium",
  24, "Developed, High",     31, "Barren Land",
  41, "Deciduous Forest",    42, "Evergreen Forest",
  43, "Mixed Forest",        52, "Shrub/Scrub",
  71, "Grassland",           81, "Pasture/Hay",
  82, "Cultivated Crops",    90, "Woody Wetlands",
  95, "Emergent Wetlands"
)

# ── 2c. LANDFIRE EVT ──────────────────────────────────────────────────────────

evt_masked_path <- file.path(cache_dir, "evt_masked.tif")

if (file.exists(evt_masked_path)) {
  cat("  [cache] Loading evt_masked.tif\n")
  evt_masked <- rast(evt_masked_path)
} else {
  cat("Downloading LANDFIRE EVT...\n")
  
  bbox_wgs84 <- tyfr_stack %>% project("EPSG:4326") %>% ext()
  aoi <- c(as.character(bbox_wgs84$xmin), as.character(bbox_wgs84$ymin),
           as.character(bbox_wgs84$xmax), as.character(bbox_wgs84$ymax))
  
  path <- tempfile(fileext = ".zip")
  
  lf <- landfireAPIv2(products = "LF2024_EVT", aoi = aoi,
                      projection = 5070, resolution = 30,
                      path = path,
                      email = "stowe.edward@gmail.com")
  
  lf_dir <- file.path(tempdir(), "lf")
  utils::unzip(path, exdir = lf_dir)
  
  evt_raw <- terra::rast(list.files(lf_dir, pattern = ".tif$", 
                               full.names = TRUE, 
                               recursive = TRUE))
  
  #evt_raw <- rast(list.files(lf_dir, pattern = "\\.tif$",
   #                          full.names = TRUE, recursive = TRUE)[1])
  
  fp_evt <- project(footprint, crs(evt_raw), method = "near")
  fp_evt <- resample(fp_evt, evt_raw, method = "near")
  evt_masked <- mask(evt_raw, fp_evt)
  
  writeRaster(evt_masked, evt_masked_path,
              datatype = "INT2U",
              gdal     = c("COMPRESS=LZW", "BIGTIFF=IF_NEEDED"),
              overwrite = TRUE)
  cat("  [saved] evt_masked.tif\n")
}

cat("EVT masked cells:", global(evt_masked, "notNA")[[1]], "\n")

# EVT attribute table
evt_attrs <- cats(evt_masked)[[1]] %>% as_tibble()
cat("EVT attribute columns:", names(evt_attrs), "\n")

# =============================================================================
# PART 3: Load and cache depth-to-groundwater
# =============================================================================

dtw_footprint <- cache_tif(
  file.path(cache_dir, "dtw_footprint.tif"),
  {
    cat("Processing DTW raster...\n")
    dtw_raw  <- rast("data/hydraulic/dtw_v3/base/LZW_DTW_base.tif")
    dtw_trim <- trim(dtw_raw)
    rm(dtw_raw); gc()
    
    # Aggregate to ~30m for alignment with land cover data
    fact     <- max(1, round(30 / res(dtw_trim)[1]))
    dtw_30m  <- aggregate(dtw_trim, fact = fact, fun = "mean", na.rm = TRUE)
    rm(dtw_trim); gc()
    
    # Reproject to raster CRS and mask to hydraulic footprint
    dtw_proj <- project(dtw_30m, rast_crs, method = "bilinear",
                        threads = TRUE)
    rm(dtw_30m); gc()
    
    fp_aligned <- resample(footprint, dtw_proj, method = "near")
    mask(dtw_proj, fp_aligned)
  },
  gdal = c("COMPRESS=LZW", "PREDICTOR=2", "TILED=YES", "BIGTIFF=IF_NEEDED")
)

cat("DTW footprint cells:", global(dtw_footprint, "notNA")[[1]], "\n")
cat("DTW value range:", as.numeric(global(dtw_footprint, "range",
                                          na.rm = TRUE)), "\n")

# Transform DTW:
#   Original: positive = surface water depth (ft), negative = depth to GW (ft)
#   Target:   negative = saturated/inundated, positive = depth to GW (m)
#   Formula:  multiply by -0.3048 (sign flip + ft to m)

dtw_footprint <- dtw_footprint * -0.3048

cat("DTW transformed summary (m, negative = saturated):\n")
print(global(dtw_footprint, c("min", "mean", "median", "max"), na.rm = TRUE))

# Check if DTW needs unit conversion (values stored ×-100 in some exports)
#dtw_median <- global(dtw_footprint, fun = "mean", na.rm = TRUE)[[1]]
#dtw_scale  <- if (abs(dtw_median) > 50) -1/100 else 1
#if (dtw_scale != 1) {
#  cat("  Applying DTW unit conversion (×", dtw_scale, ")\n")
#  dtw_footprint <- dtw_footprint * dtw_scale
#}

# DTW depth class breaks (in meters)
# Updated DTW breaks and labels (meters, sign-flipped)
dtw_breaks <- c(-Inf, -1, -0.3, 0, 0.5, 1, 2, Inf)
dtw_labels <- c("Deep surface water (< -1 m)",
                "Shallow surface water (-1–-0.3 m)",
                "Near surface (-0.3–0 m)",
                "0–0.5 m",
                "0.5–1 m",
                "1–2 m",
                "> 2 m")

dtw_colors <- c(
  "Deep surface water (< -1 m)"      = "#08306B",
  "Shallow surface water (-1–-0.3 m)"= "#2171B5",
  "Near surface (-0.3–0 m)"          = "#6BAED6",
  "0–0.5 m"                          = "#1B6B3A",
  "0.5–1 m"                          = "#7BC67E",
  "1–2 m"                            = "#D4850A",
  "> 2 m"                            = "#8B4513"
)
# =============================================================================
# PART 4: Build analysis datasets
# =============================================================================

# ── 4a. Raster sample: NLCD × EVT × inundation × DTW ─────────────────────────
# Single large sample used for all raster-based analyses

rast_sample <- cache_rds(file.path(cache_dir, "raster_sample.rds"), {
  cat("Building raster sample dataset...\n")
  
  # Align NLCD and EVT to hydraulic stack CRS and resolution
  nlcd_aligned <- project(nlcd_masked, rast_crs, method = "near",
                          threads = TRUE) %>%
    resample(aggregate(depth_tyfr, fact = agg_factor), method = "near")
  
  evt_aligned <- project(evt_masked, rast_crs, method = "near",
                         threads = TRUE) %>%
    resample(aggregate(depth_tyfr, fact = agg_factor), method = "near")
  
  dtw_aligned <- resample(dtw_footprint,
                          aggregate(depth_tyfr, fact = agg_factor),
                          method = "bilinear")
  
  soil_t_aligned <- project(texture_raster, rast_crs, method = "near",
                         threads = TRUE) %>%
    resample(aggregate(depth_tyfr, fact = agg_factor), method = "near")
  
  soil_d_aligned <- project(drainage_raster, rast_crs, method = "near",
                            threads = TRUE) %>%
    resample(aggregate(depth_tyfr, fact = agg_factor), method = "near")
  
  sample_stack <- c(
    nlcd_aligned,
    evt_aligned,
    dtw_aligned,
    soil_t_aligned,
    soil_d_aligned,
    aggregate(depth_bf,   fact = agg_factor, fun = "mean"),
    aggregate(depth_fyr,  fact = agg_factor, fun = "mean"),
    aggregate(depth_tyr,  fact = agg_factor, fun = "mean"),
    aggregate(depth_tyfr, fact = agg_factor, fun = "mean")
  )
  names(sample_stack) <- c("nlcd", "evt", "dtw", "soil_t", "soil_d",
                           "depth_bf", "depth_fyr",
                           "depth_tyr", "depth_tyfr")
  
  set.seed(42)
  rast_sample <- spatSample(sample_stack, size = n_sample, method = "random",
             na.rm = FALSE) %>%
    as_tibble() %>%
    filter(#!is.na(nlcd) | 
      !is.na(evt) | 
        !is.na(dtw)) %>%
    mutate(
      # Inundation flags
      inund_bf   = !is.na(depth_bf),
      inund_fyr  = !is.na(depth_fyr),
      inund_tyr  = !is.na(depth_tyr),
      inund_tyfr = !is.na(depth_tyfr),
      # Inundation frequency class
      inund_class = factor(case_when(
        inund_bf                   ~ "Permanent / baseflow",
        !inund_bf  & inund_fyr     ~ "5-yr event",
        !inund_fyr & inund_tyr     ~ "10-yr event",
        !inund_tyr & inund_tyfr    ~ "25-yr event",
        TRUE                       ~ "Not inundated (>25-yr)"
      ), levels = inund_levels),
      # DTW depth class
      dtw_class = factor(
        cut(dtw, breaks = dtw_breaks, labels = dtw_labels,
            include.lowest = TRUE),
        levels = dtw_labels
      ),
      # Labels
      nlcd_label = nlcd_labels$nlcd_label[
        match(nlcd, nlcd_labels$nlcd_code)
      ],
      evt_group  = evt_attrs$EVT_GP_N[match(evt, evt_attrs$Value)]
    )
})

cat("Raster sample rows:", nrow(rast_sample), "\n")

# ── 4b. Polygon extraction: LTRM polygons × hydraulics × DTW ──────────────────

extract_poly <- function(lulc_pool, stack, scenario_label, pool_id) {
  cat("  Pool", pool_id, scenario_label, "...\n")
  v  <- vect(lulc_pool)
  ex <- terra::extract(stack, v, fun = mean, na.rm = TRUE, bind = FALSE)
  ex_cov <- terra::extract(stack[["depth"]], v,
                           fun = function(x, ...) mean(!is.na(x)),
                           bind = FALSE)
  lulc_pool %>%
    st_drop_geometry() %>%
    mutate(depth_mean    = ex[, 2],
           velocity_mean = ex[, 3],
           froude_mean   = ex[, 4],
           pct_covered   = ex_cov[, 2],
           scenario      = scenario_label)
}

lulc_hydro_raw <- cache_rds(file.path(cache_dir, "lulc_hydro_raw.rds"), {
  cat("Extracting hydraulics to LTRM polygons (slow — cached after first run)...\n")
  # Aggregate stacks before extraction
  bf_agg   <- aggregate(bf_stack,   fact = 1, fun = "mean")
  fyr_agg  <- aggregate(fyr_stack,  fact = 1, fun = "mean")
  tyr_agg  <- aggregate(tyr_stack,  fact = 1, fun = "mean")
  tyfr_agg <- aggregate(tyfr_stack, fact = 1, fun = "mean")
  
  map_dfr(unique(lulc_clip$pool), function(p) {
    pl <- filter(lulc_clip, pool == p)
    bind_rows(
      extract_poly(pl, bf_agg,   scenario_labels["bf"],   p),
      extract_poly(pl, fyr_agg,  scenario_labels["fyr"],  p),
      extract_poly(pl, tyr_agg,  scenario_labels["tyr"],  p),
      extract_poly(pl, tyfr_agg, scenario_labels["tyfr"], p)
    )
  })
})

# DTW extraction at polygon centroids (faster than polygon mean for DTW)
lulc_dtw <- cache_rds(file.path(cache_dir, "lulc_dtw.rds"), {
  cat("Extracting DTW to LTRM polygon centroids...\n")
  centroids <- lulc_clip %>%
    st_centroid() %>%
    st_transform(crs = crs(dtw_footprint))
  ex <- terra::extract(dtw_footprint, vect(centroids), bind = FALSE)
  lulc_clip %>%
    st_drop_geometry() %>%
    select(poly_id, lc_class) %>%
    mutate(dtw = ex[, 2])
})

# Clean: mask low-coverage hydraulics, join DTW, add inundation class
lulc_hydro <- lulc_hydro_raw %>%
  mutate(
    depth_mean    = if_else(pct_covered < min_coverage, NA_real_, depth_mean),
    velocity_mean = if_else(pct_covered < min_coverage, NA_real_, velocity_mean),
    froude_mean   = if_else(pct_covered < min_coverage, NA_real_, froude_mean),
    scenario      = factor(scenario, levels = unname(scenario_labels))
  )

# Inundation class per polygon
inundation <- lulc_hydro %>%
  select(poly_id, lc_class, pool, scenario, pct_covered) %>%
  mutate(inundated = pct_covered >= min_coverage) %>%
  group_by(poly_id, lc_class, pool) %>%
  summarize(
    inund_bf   = any(inundated[scenario == scenario_labels["bf"]],   na.rm = TRUE),
    inund_fyr  = any(inundated[scenario == scenario_labels["fyr"]],  na.rm = TRUE),
    inund_tyr  = any(inundated[scenario == scenario_labels["tyr"]],  na.rm = TRUE),
    inund_tyfr = any(inundated[scenario == scenario_labels["tyfr"]], na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(
    inund_class = factor(case_when(
      inund_bf                  ~ "Permanent / baseflow",
      !inund_bf  & inund_fyr   ~ "5-yr event",
      !inund_fyr & inund_tyr   ~ "10-yr event",
      !inund_tyr & inund_tyfr  ~ "25-yr event",
      TRUE                      ~ "Not inundated (>25-yr)"
    ), levels = inund_levels)
  ) %>%
  left_join(lulc_dtw, by = c("poly_id", "lc_class")) %>%
  mutate(
    dtw_class = factor(
      cut(dtw, breaks = dtw_breaks, labels = dtw_labels,
          include.lowest = TRUE),
      levels = dtw_labels
    )
  )

# Order LTRM classes by inundation (most to least)
class_order <- inundation %>%
  count(lc_class, inund_class) %>%
  group_by(lc_class) %>%
  mutate(pct = n / sum(n)) %>%
  filter(inund_class == "Permanent / baseflow") %>%
  arrange(pct) %>%
  pull(lc_class)
class_order <- c(class_order,
                 setdiff(unique(inundation$lc_class), class_order))

cat("LTRM polygon extraction complete.\n")
cat("Inundation class counts:\n")
print(table(inundation$inund_class))

# =============================================================================
# PART 5: Inundation frequency analysis
# =============================================================================

# ── 5a. LTRM LULC × inundation ───────────────────────────────────────────────

ltrm_inund_props <- inundation %>%
  mutate(lc_class = factor(lc_class, levels = class_order)) %>%
  filter(!str_detect(lc_class, "Developed|Road")) %>%
  count(lc_class, inund_class) %>%
  group_by(lc_class) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

# ── 5b. NLCD × inundation ────────────────────────────────────────────────────

nlcd_inund_props <- rast_sample %>%
  filter(!is.na(nlcd_label),
         !str_detect(nlcd_label, "Developed")) %>%
  count(nlcd_label, inund_class) %>%
  group_by(nlcd_label) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

nlcd_order <- nlcd_inund_props %>%
  filter(inund_class == "Not inundated (>25-yr)") %>%
  arrange(pct) %>% pull(nlcd_label)
nlcd_order <- c(nlcd_order,
                setdiff(unique(nlcd_inund_props$nlcd_label), nlcd_order))

# ── 5c. LANDFIRE EVT × inundation ────────────────────────────────────────────

evt_grp_min_n <- 20   # minimum pixels to include a group

evt_inund_props <- rast_sample %>%
  filter(!is.na(evt_group),
         !str_detect(evt_group, "Developed|Crop|Urban")) %>%
  count(evt_group, inund_class) %>%
  group_by(evt_group) %>%
  mutate(pct = n / sum(n) * 100,
         total_n = sum(n)) %>%
  ungroup() %>%
  filter(total_n >= evt_grp_min_n)

evt_order <- evt_inund_props %>%
  filter(inund_class == "Not inundated (>25-yr)") %>%
  arrange(pct) %>% pull(evt_group)
evt_order <- c(evt_order,
               setdiff(unique(evt_inund_props$evt_group), evt_order))

# ── 5d. Statistical tests ─────────────────────────────────────────────────────

cat("\nChi-square: inundation ~ LTRM class\n")
print(chisq.test(table(inundation$lc_class, inundation$inund_class)))

# Kruskal-Wallis effect size across scenarios
kw_effect <- function(data, var, group) {
  kw   <- kruskal.test(as.formula(paste(var, "~", group)), data = data)
  k    <- n_distinct(data[[group]], na.rm = TRUE)
  n    <- sum(!is.na(data[[var]]))
  tibble(variable = var, H = round(kw$statistic, 1), k = k, n = n,
         eta_sq = round((kw$statistic - k + 1) / (n - k), 3))
}

kw_results <- map_dfr(levels(lulc_hydro$scenario), function(sc) {
  d <- filter(lulc_hydro, scenario == sc, !is.na(depth_mean),
              !is.na(lc_class))
  bind_rows(
    kw_effect(d, "depth_mean",    "lc_class"),
    kw_effect(d, "velocity_mean", "lc_class")
  ) %>% mutate(scenario = sc)
})
cat("\nKruskal-Wallis eta-squared:\n"); print(kw_results)

# =============================================================================
# PART 6: DTW analysis
# =============================================================================

# ── 6a. DTW × LTRM LULC ──────────────────────────────────────────────────────

ltrm_dtw_summary <- inundation %>%
  filter(!is.na(dtw), !str_detect(lc_class, "Developed|Road")) %>%
  mutate(lc_class = factor(lc_class, levels = class_order)) %>%
  group_by(lc_class) %>%
  summarize(
    n             = n(),
    dtw_median    = round(median(dtw, na.rm = TRUE), 2),
    dtw_q25       = round(quantile(dtw, 0.25, na.rm = TRUE), 2),
    dtw_q75       = round(quantile(dtw, 0.75, na.rm = TRUE), 2),
    pct_saturated = round(mean(dtw < 0, na.rm = TRUE) * 100, 1),
    .groups       = "drop"
  ) %>%
  arrange(dtw_median)

cat("\nDTW summary by LTRM class:\n"); print(ltrm_dtw_summary)

ltrm_dtw_props <- inundation %>%
  filter(!is.na(dtw_class), !str_detect(lc_class, "Developed|Road")) %>%
  mutate(lc_class = factor(lc_class, levels = class_order)) %>%
  count(lc_class, dtw_class) %>%
  group_by(lc_class) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

# ── 6b. DTW × NLCD and EVT (from raster sample) ──────────────────────────────

nlcd_dtw_props <- rast_sample %>%
  filter(!is.na(nlcd_label), !is.na(dtw_class),
         !str_detect(nlcd_label, "Developed")) %>%
  count(nlcd_label, dtw_class) %>%
  group_by(nlcd_label) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

evt_dtw_props <- rast_sample %>%
  filter(!is.na(evt_group), !is.na(dtw_class),
         !str_detect(evt_group, "Developed|Crop|Urban")) %>%
  count(evt_group, dtw_class) %>%
  group_by(evt_group) %>%
  mutate(pct = n / sum(n) * 100, total_n = sum(n)) %>%
  ungroup() %>%
  filter(total_n >= evt_grp_min_n)

# Kruskal-Wallis: DTW ~ land cover class
kw_dtw_ltrm <- kruskal.test(dtw ~ lc_class, data = inundation)
k  <- n_distinct(inundation$lc_class, na.rm = TRUE)
n  <- sum(!is.na(inundation$dtw))
eta2_dtw_ltrm <- (kw_dtw_ltrm$statistic - k + 1) / (n - k)

cat("\nKruskal-Wallis: DTW ~ LTRM class\n")
cat("H =", round(kw_dtw_ltrm$statistic, 1),
    " eta-sq =", round(eta2_dtw_ltrm, 3), "\n")

# =============================================================================
# LANDFIRE EVT classes in DTW × Inundation frequency space
# =============================================================================

library(ggrepel)   # for non-overlapping labels

# ── 1. Build summary: median DTW and inundation score per EVT group ───────────
# Assign a numeric score to inundation frequency (lower = wetter)
inund_score <- c(
  "Permanent / baseflow"   = 1,
  "5-yr event"             = 2,
  "10-yr event"            = 3,
  "25-yr event"            = 4,
  "Not inundated (>25-yr)" = 5
)

# Higher resolution
evt_space <- rast_sample %>%
  filter(!is.na(evt),
         !is.na(dtw),
         !is.na(inund_class),
         !str_detect(evt, "Developed|Crop|Urban")) %>%
  mutate(inund_score = inund_score[as.character(inund_class)]) %>%
  group_by(evt) %>%
  summarize(
    n             = n(),
    dtw_median    = median(dtw,         na.rm = TRUE),
    dtw_q25       = quantile(dtw, 0.25, na.rm = TRUE),
    dtw_q75       = quantile(dtw, 0.75, na.rm = TRUE),
    inund_median  = median(inund_score, na.rm = TRUE),
    inund_q25     = quantile(inund_score, 0.25, na.rm = TRUE),
    inund_q75     = quantile(inund_score, 0.75, na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  filter(n >= 20) %>%   # drop very rare classes
  mutate(
    # Scale rectangle size: area proportional to log(n) for readability
    rect_size = sqrt(n / max(n)),
    # Half-widths for geom_rect
    x_half    = (inund_q75 - inund_q25) / 2,
    y_half    = (dtw_q75   - dtw_q25)   / 2,
    # Short label for display
    label_short = str_wrap(evt, 20)
  )

# ── 2. Assign broad color groups ──────────────────────────────────────────────
# Color by vegetation life form — adjust regex to match your EVT_GP_N values
evt_space <- evt_space %>%
  mutate(
    life_form = case_when(
      str_detect(evt, regex("forest|woodland", ignore_case = TRUE))   ~ "Forest/Woodland",
      str_detect(evt, regex("shrub|scrub",     ignore_case = TRUE))   ~ "Shrubland",
      str_detect(evt, regex("herb|grass|sedge|meadow",
                                  ignore_case = TRUE))                     ~ "Herbaceous",
      str_detect(evt, regex("riparian|floodplain",
                                  ignore_case = TRUE))                     ~ "Riparian",
      str_detect(evt, regex("wetland|marsh|bog|fen",
                                  ignore_case = TRUE))                     ~ "Wetland",
      str_detect(evt, regex("water|aquatic",   ignore_case = TRUE))   ~ "Open Water",
      TRUE                                                                   ~ "Other"
    ),
    life_form = factor(life_form,
                       levels = c("Open Water", "Wetland", "Riparian",
                                  "Forest/Woodland", "Shrubland",
                                  "Herbaceous", "Other"))
  ) %>%
  mutate(evt2 = str_remove(evt, "North-Central Interior|Eastern Great Plains|Eastern Cool Temperate|Northern & Central"),
         dtw_q25 = case_when(evt == "Open Water" ~ -2.4,
                             TRUE ~ dtw_q25)) #%>%
  #left_join(evt_attrs, by = c("evt"= "EVT_NAME")) %>%
  #mutate(evt2 = paste(EVT_PHYS, life_form))


# Lower resolution
evt_space <- rast_sample %>%
  filter(!is.na(evt),
         !is.na(dtw),
         !is.na(inund_class),
         !str_detect(evt, "Developed|Crop|Urban|Agri|Pasture|Hayland")) %>%
  mutate(inund_score = inund_score[as.character(inund_class)]) %>%
  left_join(evt_attrs, by = c("evt"= "EVT_NAME")) %>%
  mutate(EVT_GP_N = paste(EVT_PHYS, EVT_CLASS)) %>%
  group_by(EVT_GP_N) %>%
  summarize(
    n             = n(),
    dtw_median    = median(dtw,         na.rm = TRUE),
    dtw_q25       = quantile(dtw, 0.25, na.rm = TRUE),
    dtw_q75       = quantile(dtw, 0.75, na.rm = TRUE),
    inund_median  = median(inund_score, na.rm = TRUE),
    inund_q25     = quantile(inund_score, 0.25, na.rm = TRUE),
    inund_q75     = quantile(inund_score, 0.75, na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  filter(n >= 20) %>%   # drop very rare classes
  mutate(
    # Scale rectangle size: area proportional to log(n) for readability
    rect_size = sqrt(n / max(n)),
    # Half-widths for geom_rect
    x_half    = (inund_q75 - inund_q25) / 2,
    y_half    = (dtw_q75   - dtw_q25)   / 2,
    # Short label for display
    label_short = str_wrap(EVT_GP_N, 20)
  )

# ── 2. Assign broad color groups ──────────────────────────────────────────────
# Color by vegetation life form — adjust regex to match your EVT_GP_N values
evt_space <- evt_space %>%
  mutate(
    life_form = case_when(
      str_detect(EVT_GP_N, regex("forest|woodland|hardwood", ignore_case = TRUE))   ~ "Forest/Woodland",
      str_detect(EVT_GP_N, regex("shrub|scrub",     ignore_case = TRUE))   ~ "Shrubland",
      str_detect(EVT_GP_N, regex("herb|grass|sedge|meadow",
                            ignore_case = TRUE))                     ~ "Herbaceous",
      str_detect(EVT_GP_N, regex("riparian|floodplain",
                            ignore_case = TRUE))                     ~ "Riparian",
      str_detect(EVT_GP_N, regex("wetland|marsh|bog|fen",
                            ignore_case = TRUE))                     ~ "Wetland",
      str_detect(EVT_GP_N, regex("water|aquatic",   ignore_case = TRUE))   ~ "Open Water",
      TRUE                                                                   ~ "Other"
    ),
    life_form = factor(life_form,
                       levels = c("Open Water", "Wetland", "Riparian",
                                  "Forest/Woodland", "Shrubland",
                                  "Herbaceous", "Other"))
  ) %>%
  mutate(evt2 = case_when(str_detect(EVT_GP_N, "Open Water") ~ "Open Water",
                          str_detect(EVT_GP_N, "Exotic Herb") ~ "Exotic Meadow",
                          str_detect(EVT_GP_N, "Riparian Open") ~ "Floodplain Forest",
                          str_detect(EVT_GP_N, "Riparian Closed") ~ "Swamp Forest",
                          str_detect(EVT_GP_N, "Hardwood Closed") ~ "Upland Hardwood",
                          str_detect(EVT_GP_N, "Conifer-Hardwood") ~ "Conifer-Hardwood",
                          str_detect(EVT_GP_N, "Hardwood Open") ~ "Savanna",
                          TRUE ~ EVT_GP_N),
         evt2 = str_replace(evt2, "Herbaceous - grassland", ~ "Herbaceous"),
         dtw_q25 = case_when(EVT_GP_N == "Open Water Non-vegetated" ~ -2.4,
                             TRUE ~ dtw_q25)
         )
  

life_form_colors <- c(
  "Open Water"      = "#08519C",
  "Wetland"         = "#2CA25F",
  "Riparian"        = "#74C476",
  "Forest/Woodland" = "#006D2C",
  "Herbaceous"       = "#FEC44F",
  #"Herbaceous"      = "#ADDD8E",
  "Other"           = "#BDBDBD"
)

# ── 3. X-axis: inundation frequency labels ────────────────────────────────────
inund_x_labels <- c(
  "1" = "Permanent\n(baseflow)",
  "2" = "5-yr\nevent",
  "3" = "10-yr\nevent",
  "4" = "25-yr\nevent",
  "5" = ">25-yr\n(rare/never)"
)

# ── 4. Plot ───────────────────────────────────────────────────────────────────

ggplot(evt_space) +
  # Background shading to guide interpretation
  annotate("rect", xmin = 0.5, xmax = 1.8, ymin = -Inf, ymax = Inf,
           fill = "#EBF5FB", alpha = 1) +
  #annotate("text", x = 1.15, y = max(evt_space$dtw_q75, na.rm = TRUE) * 0.95,
   #        label = "Permanently\nwet", size = 3, color = "#2471A3",
    #       fontface = "italic") +
  
  # Rectangles: IQR of inundation score (x) × IQR of DTW (y)
  # Width and height reflect spread; alpha reflects relative frequency
  geom_rect(aes(
    xmin  = inund_median - 2*x_half * rect_size,
    xmax  = inund_median + 2*x_half * rect_size,
    ymin  = -(dtw_q75),   # flip: deeper GW = lower on plot
    ymax  = -(dtw_q25),
    fill  = life_form,
    alpha = rect_size
  ), color = "black", linewidth = 0.5) +
  
  # Centre point
  geom_jitter(aes(x = inund_median, y = -dtw_median,
                 size = n, color = life_form), width = .05,
             shape = 21, fill = "white", stroke = 1.2) +
  
  # Non-overlapping labels
  geom_label_repel(aes(x = inund_median, y = -dtw_median,
                       label = evt2, #label_short,
                       ), #life_form),
                   color = "black",
                   size           = 3.5,
                   fill           = alpha("white", 0.8),
                   label.size     = 0.2,
                   box.padding    = 0.4,
                   point.padding  = 0.3,
                   max.overlaps   = 20,
                   segment.color  = "grey60",
                   segment.size   = 0.3,
                   show.legend    = FALSE) +
  
  # Axes
  scale_x_continuous(
    #limits = c(0.5, 5.5),
    breaks = 1:5,
    labels = inund_x_labels,
    position = "top",
    #name   = expression("← Wetter            "~ bold("Inundation frequency") ~"            Drier →")
  ) +
  scale_y_continuous(
    #name   = expression("↑ Deeper groundwater          "~  bold("Depth to groundwater (m)") ~ "            Surface water"), # ↓",
    labels = ~ abs(.x),   # show positive depth values on flipped axis,
    expand = expansion(mult = 0)
  ) +
  scale_fill_manual(values = life_form_colors,  name = "Vegetation type") +
  scale_color_manual(values = life_form_colors, name = "Vegetation type") +
  scale_alpha_continuous(range = c(0.3, 0.8), guide = "none") +
  #scale_size_continuous(
  #  name   = "Relative area\n(n pixels)",
  #  range  = c(1, 6),
  #  labels = scales::label_number(scale = 1e-3, suffix = "k")
  #) +
  
  # Reference line: ground surface
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.6) +
  #annotate("text", x = 1.3, y = 0.15, label = "Ground surface",
   #        size = 2.8, color = "grey40", hjust = 1) +
  
  labs(
    x = "Inundation frequency" , y = "Depth to groundwater (m)",
    #title    = "LANDFIRE EVT classes in hydrologic space",
    #subtitle = "Position = median inundation frequency × median depth to groundwater\nRectangle height = IQR of DTW  ·  Size scaled by relative area"
  ) +
  theme_classic(base_size = 16) +
  theme(
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    #axis.title.x       = element_text(size = 10, margin = margin(t = 8)),
    #axis.title.y       = element_text(size = 10, margin = margin(r = 8)),
    legend.position    = "None",
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 9, color = "grey40")
  )

ggsave("output/evt_hydrologic_space.png",
       width = 12, height = 8, dpi = 200)




# =============================================================================
# PART 7: Combined inundation + DTW models
# =============================================================================
# Do inundation frequency and DTW together explain more about land cover
# than either alone? Use multinomial logistic regression as a simple test.

library(nnet)

# Prepare combined dataset from raster sample
model_dat <- rast_sample %>%
  filter(!is.na(nlcd), !is.na(inund_class), !is.na(dtw)) %>%
  mutate(
    nlcd_label  = factor(nlcd),
    inund_class = factor(inund_class, levels = inund_levels),
    # Numeric predictors for regression
    inund_num  = as.numeric(inund_class),  # ordinal encoding
    dtw_scaled = scale(dtw)[, 1]
  )

set.seed(42)
train_idx  <- sample(nrow(model_dat), 0.7 * nrow(model_dat))
train_dat  <- model_dat[ train_idx, ]
test_dat   <- model_dat[-train_idx, ]

# Model 1: inundation only
m_inund <- multinom(nlcd ~ inund_num, data = train_dat,
                    trace = FALSE, MaxNWts = 10000)

# Model 2: DTW only
m_dtw   <- multinom(nlcd ~ dtw_scaled, data = train_dat,
                    trace = FALSE, MaxNWts = 10000)

# Model 3: inundation + DTW
m_both  <- multinom(nlcd ~ inund_num + dtw_scaled, data = train_dat,
                    trace = FALSE, MaxNWts = 10000)

preds <- predict(m_both, test_dat)

test_w_preds <- test_dat %>% mutate(preds = predict(m_both, test_dat)) %>%
  relocate(preds, .after = nlcd) %>%
  mutate(nlcd = as.character(nlcd),
         preds = as.character(preds),
         pred_success = nlcd == preds) %>%
  group_by(nlcd, pred_success) %>%
  tally()

# Test set accuracy
acc <- function(model, test) {
  preds <- predict(model, newdata = test)
  # Convert both to character before comparing to avoid factor level mismatch
  mean(as.character(preds) == as.character(test$nlcd), na.rm = TRUE)
}

model_comparison <- tibble(
  model    = c("Inundation only", "DTW only", "Inundation + DTW"),
  accuracy = c(acc(m_inund, test_dat),
               acc(m_dtw,   test_dat),
               acc(m_both,  test_dat)),
  AIC      = c(AIC(m_inund), AIC(m_dtw), AIC(m_both))
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat("\nModel comparison — predicting NLCD class:\n")
print(model_comparison)

### EVT

# Prepare combined dataset from raster sample
model_dat <- rast_sample %>%
  filter(!is.na(evt), !is.na(inund_class), !is.na(dtw),
         !str_detect(evt, "Developed|Crop|Urban|Agri|Pasture|Hayland")) %>%
  mutate(
    evt  = factor(evt),
    inund_class = factor(inund_class, levels = inund_levels),
    # Numeric predictors for regression
    inund_num  = as.numeric(inund_class),  # ordinal encoding
    dtw_scaled = scale(dtw)[, 1]
  )

set.seed(42)
train_idx  <- sample(nrow(model_dat), 0.7 * nrow(model_dat))
train_dat  <- model_dat[ train_idx, ]
test_dat   <- model_dat[-train_idx, ]

# Model 1: inundation only
m_inund <- multinom(evt ~ inund_num, data = train_dat,
                    trace = FALSE, MaxNWts = 10000)

# Model 2: DTW only
m_dtw   <- multinom(evt ~ dtw_scaled, data = train_dat,
                    trace = FALSE, MaxNWts = 10000)

# Model 3: inundation + DTW
m_both  <- multinom(evt ~ inund_num + dtw_scaled, data = train_dat,
                    trace = FALSE, MaxNWts = 10000)

# Test set accuracy
acc <- function(model, test) {
  preds <- predict(model, newdata = test)
  # Convert both to character before comparing to avoid factor level mismatch
  mean(as.character(preds) == as.character(test$evt), na.rm = TRUE)
}

model_comparison <- tibble(
  model    = c("Inundation only", "DTW only", "Inundation + DTW"),
  accuracy = c(acc(m_inund, test_dat),
               acc(m_dtw,   test_dat),
               acc(m_both,  test_dat)),
  AIC      = c(AIC(m_inund), AIC(m_dtw), AIC(m_both))
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat("\nModel comparison — predicting evt class:\n")
print(model_comparison)

# Same thing but with random forest
library(randomForest)

# Prepare combined dataset from raster sample
model_dat <- rast_sample %>%
  filter(
    evt !="Open Water",
    !is.na(soil_d),
    !is.na(soil_t),
    !is.na(evt),
    !is.na(inund_class),
    !is.na(dtw),
    !str_detect(evt, "Developed|Crop|Urban|Agri|Pasture|Hayland")
  ) %>%
  left_join(evt_attrs, by = c("evt"="EVT_NAME")) %>%
  mutate(
    evt = factor(EVT_GP_N),
    inund_class = factor(inund_class, levels = inund_levels),
    inund_num = as.numeric(inund_class),
    dtw_scaled = as.numeric(scale(dtw)),
  )

library(caret)

set.seed(42)

train_idx <- createDataPartition(
  model_dat$evt,
  p = 0.7,
  list = FALSE
)[, 1]

train_dat <- droplevels(model_dat[train_idx, ])
test_dat  <- model_dat[-train_idx, ]

# Fit random forest models
set.seed(42)

# Model 1: inundation only
m_inund <- randomForest(
  evt ~ inund_num,
  data = train_dat,
  ntree = 500,
  importance = TRUE
)

# Model 2: DTW only
m_dtw <- randomForest(
  evt ~ dtw_scaled,
  data = train_dat,
  ntree = 500,
  importance = TRUE
)

# Model 3: inundation + DTW
m_both <- randomForest(
  evt ~ inund_num + dtw_scaled + soil_d + soil_t,
  data = train_dat,
  ntree = 500,
  importance = TRUE
)

preds <- predict(m_both, test_dat)

test_w_preds <- test_dat %>% mutate(preds = predict(m_both, test_dat)) %>%
  relocate(preds, .after = evt) %>%
  mutate(pred_success = evt == preds) %>%
  group_by(evt, pred_success) %>%
  tally()

# Assess
library(yardstick)

# Calculate test-set accuracy
acc <- function(model, test) {
  preds <- predict(model, newdata = test)
  mean(as.character(preds) == as.character(test$evt))
}

# Calculate multiclass AUC using predicted class probabilities
auc_multiclass <- function(model, test) {
  
  probs <- predict(model, newdata = test, type = "prob")
  
  # Keep only classes represented in the test set and probability matrix
  test_classes <- intersect(
    unique(as.character(test$evt)),
    colnames(probs)
  )
  
  test_auc <- test %>%
    mutate(evt = factor(evt, levels = colnames(probs))) %>%
    filter(as.character(evt) %in% test_classes) %>%
    select(evt)
  
  probs <- probs[as.character(test$evt) %in% test_classes, test_classes,
                 drop = FALSE]
  
  # Rename probability columns for yardstick
  prob_dat <- as.data.frame(probs)
  names(prob_dat) <- paste0(".pred_", make.names(test_classes))
  
  test_auc <- test_auc %>%
    mutate(across(evt, ~ factor(.x, levels = test_classes))) %>%
    bind_cols(prob_dat)
  
  # Macro-averaged multiclass ROC AUC
  yardstick::roc_auc(
    test_auc,
    truth = evt,
    !!!rlang::syms(names(prob_dat)),
    estimator = "hand_till"
  )$.estimate
}

model_comparison <- tibble(
  model = c(
    "Inundation only",
    "DTW only",
    "Inundation + DTW"
  ),
  accuracy = c(
    acc(m_inund, test_dat),
    acc(m_dtw, test_dat),
    acc(m_both, test_dat)
  ),
  AUC = c(
    auc_multiclass(m_inund, test_dat),
    auc_multiclass(m_dtw, test_dat),
    auc_multiclass(m_both, test_dat)
  ),
  OOB_error = c(
    m_inund$err.rate[m_inund$ntree, "OOB"],
    m_dtw$err.rate[m_dtw$ntree, "OOB"],
    m_both$err.rate[m_both$ntree, "OOB"]
  )
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

model_comparison

importance(m_both)

varImpPlot(m_both)
# =============================================================================
# PART 8: Cross-classification of land cover systems
# =============================================================================

# Rasterize LTRM to NLCD grid for cross-tab
lulc_rast <- cache_tif(
  file.path(cache_dir, "lulc_rasterized.tif"),
  {
    lc_lev  <- tibble(lc_class = sort(unique(lulc_clip$lc_class))) %>%
      mutate(lc_code = row_number())
    lulc_nlcd <- st_transform(lulc_clip, crs = crs(nlcd_masked)) %>%
      left_join(lc_lev, by = "lc_class")
    r <- rasterize(vect(lulc_nlcd), nlcd_masked,
                   field = "lc_code", fun = "max")
    names(r) <- "ltrm_code"
    r
  }
)

# Store lc_levels for joining later
lc_levels <- tibble(lc_class = sort(unique(lulc_clip$lc_class))) %>%
  mutate(lc_code = row_number())

cross_sample <- cache_rds(file.path(cache_dir, "cross_sample.rds"), {
  cat("Building cross-classification sample...\n")
  evt_aligned_nlcd <- resample(evt_masked, nlcd_masked, method = "near")
  names(evt_aligned_nlcd) <- "evt_code"
  names(nlcd_masked)      <- "nlcd_code"
  
  lc_stack <- c(lulc_rast, nlcd_masked, evt_aligned_nlcd)
  set.seed(42)
  spatSample(lc_stack, size = n_sample, method = "random",
             na.rm = FALSE) %>%
    as_tibble() %>%
    filter(!is.na(ltrm_code) | !is.na(nlcd_code) | !is.na(evt_code)) %>%
    left_join(lc_levels,   by = c("ltrm_code" = "lc_code")) %>%
    left_join(nlcd_labels, by = "nlcd_code") %>%
    left_join(evt_attrs %>% select(evt_code = Value, evt_group = EVT_GP_N),
              by = "evt_code")
})

# Pairwise cross-tabs
cross_ltrm_nlcd <- cross_sample %>%
  filter(!is.na(lc_class), !is.na(nlcd_label)) %>%
  count(lc_class, nlcd_label) %>%
  group_by(lc_class) %>%
  mutate(pct = n / sum(n) * 100) %>% ungroup()

cross_ltrm_evt <- cross_sample %>%
  filter(!is.na(lc_class), !is.na(evt_group)) %>%
  count(lc_class, evt_group) %>%
  group_by(lc_class) %>%
  mutate(pct = n / sum(n) * 100) %>% ungroup()

cross_nlcd_evt <- cross_sample %>%
  filter(!is.na(nlcd_label), !is.na(evt_group)) %>%
  count(nlcd_label, evt_group) %>%
  group_by(nlcd_label) %>%
  mutate(pct = n / sum(n) * 100) %>% ungroup()

cat("\nTop NLCD class per LTRM class:\n")
cross_ltrm_nlcd %>% group_by(lc_class) %>%
  slice_max(pct, n = 1) %>% print(n = 50)

# =============================================================================
# PART 9: Plots
# =============================================================================

# Helper: stacked bar plot for inundation frequency
inund_bar <- function(data, x_var, title, filter_str = NULL) {
  if (!is.null(filter_str))
    data <- filter(data, !str_detect(.data[[x_var]], filter_str))
  ord <- data %>% filter(inund_class == "Not inundated (>25-yr)") %>%
    arrange(pct) %>% pull(.data[[x_var]])
  ord <- c(ord, setdiff(unique(data[[x_var]]), ord))
  data %>%
    mutate(across(all_of(x_var), ~ factor(.x, levels = rev(ord)))) %>%
    ggplot(aes(x = .data[[x_var]], y = pct, fill = inund_class)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = inund_colors, name = "Inundation frequency") +
    scale_y_continuous(labels = scales::label_number(suffix = "%")) +
    coord_flip() +
    labs(title = title, x = NULL, y = "Proportion of area (%)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
          legend.text = element_text(size = 8)) +
    guides(fill = guide_legend(nrow = 2))
}

# Helper: stacked bar for DTW
dtw_bar <- function(data, x_var, title, filter_str = NULL) {
  if (!is.null(filter_str))
    data <- filter(data, !str_detect(.data[[x_var]], filter_str))
  dtw_ord <- data %>%
    group_by(.data[[x_var]]) %>%
    summarize(med_dtw = weighted.mean(
      as.numeric(dtw_class), n, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(med_dtw)) %>% pull(.data[[x_var]])
  data %>%
    mutate(across(all_of(x_var), ~ factor(.x, levels = rev(dtw_ord)))) %>%
    ggplot(aes(x = .data[[x_var]], y = pct, fill = dtw_class)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = dtw_colors, name = "Depth to groundwater") +
    scale_y_continuous(labels = scales::label_number(suffix = "%")) +
    coord_flip() +
    labs(title = title, x = NULL, y = "Proportion of area (%)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
          legend.text = element_text(size = 8)) +
    guides(fill = guide_legend(nrow = 2))
}

# Helper: cross-classification heatmap
cross_heatmap <- function(data, x_var, y_var, pct_var, title,
                          palette = "YlGnBu") {
  data %>%
    mutate(across(c(all_of(x_var), all_of(y_var)), ~ str_wrap(.x, 20))) %>%
    ggplot(aes(x = .data[[x_var]], y = .data[[y_var]],
               fill = .data[[pct_var]])) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = if_else(.data[[pct_var]] >= 5,
                                  paste0(round(.data[[pct_var]]), "%"), "")),
              size = 2.8, color = "white", fontface = "bold") +
    scale_fill_distiller(palette = palette, direction = 1,
                         name = "% overlap", limits = c(0, 100)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y  = element_text(size = 8),
          panel.grid   = element_blank(),
          legend.position = "right")
}

# ── Generate all plots ────────────────────────────────────────────────────────

p_ltrm_inund <- inund_bar(ltrm_inund_props, "lc_class",
                          "Inundation frequency — LTRM land cover")
p_nlcd_inund <- inund_bar(nlcd_inund_props, "nlcd_label",
                          "Inundation frequency — NLCD",
                          "Developed")
p_evt_inund  <- inund_bar(evt_inund_props,  "evt_group",
                          "Inundation frequency — LANDFIRE EVT")

p_ltrm_dtw   <- dtw_bar(ltrm_dtw_props, "lc_class",
                        "DTW class — LTRM land cover")
p_nlcd_dtw   <- dtw_bar(nlcd_dtw_props, "nlcd_label",
                        "DTW class — NLCD", "Developed")
p_evt_dtw    <- dtw_bar(evt_dtw_props,  "evt_group",
                        "DTW class — LANDFIRE EVT")

p_cross_ltrm_nlcd <- cross_heatmap(
  filter(cross_ltrm_nlcd, !str_detect(lc_class, "Developed|Road|Agri"),
         !str_detect(nlcd_label, "Developed|Crop")),
  "nlcd_label", "lc_class", "pct",
  "LTRM × NLCD cross-classification", "YlGnBu")

p_cross_ltrm_evt <- cross_heatmap(
  filter(cross_ltrm_evt, !str_detect(lc_class, "Developed|Road|Agri"),
         !str_detect(evt_group, "Developed|Crop")),
  "evt_group", "lc_class", "pct",
  "LTRM × LANDFIRE EVT cross-classification", "YlOrRd")

p_cross_nlcd_evt <- cross_heatmap(
  filter(cross_nlcd_evt, !str_detect(nlcd_label, "Developed|Crop"),
         !str_detect(evt_group, "Developed|Crop")),
  "evt_group", "nlcd_label", "pct",
  "NLCD × LANDFIRE EVT cross-classification", "PuBuGn")

# Kruskal-Wallis effect size progression
p_kw_effect <- kw_results %>%
  mutate(
    scenario = factor(scenario, levels = unname(scenario_labels)),
    variable = recode(variable,
                      depth_mean    = "Depth",
                      velocity_mean = "Velocity")
  ) %>%
  ggplot(aes(x = scenario, y = eta_sq, color = variable, group = variable)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c(Depth = "#1B3A5C", Velocity = "#D4850A"),
                     name = NULL) +
  scale_y_continuous(limits = c(0, NA),
                     labels = scales::label_number(accuracy = 0.01)) +
  labs(title    = "Hydraulic discrimination of land cover across flow scenarios",
       subtitle = "Kruskal-Wallis eta-squared",
       x = "Flow scenario", y = "Eta-squared") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 20, hjust = 1))

# ── Save all plots ────────────────────────────────────────────────────────────

ggsave("output/inund_ltrm.png",     p_ltrm_inund,   width = 9, height = 7,  dpi = 150)
ggsave("output/inund_nlcd.png",     p_nlcd_inund,   width = 9, height = 6,  dpi = 150)
ggsave("output/inund_evt.png",      p_evt_inund,    width = 9, height = 8,  dpi = 150)
ggsave("output/dtw_ltrm.png",       p_ltrm_dtw,     width = 9, height = 7,  dpi = 150)
ggsave("output/dtw_nlcd.png",       p_nlcd_dtw,     width = 9, height = 6,  dpi = 150)
ggsave("output/dtw_evt.png",        p_evt_dtw,      width = 9, height = 8,  dpi = 150)
ggsave("output/cross_ltrm_nlcd.png",p_cross_ltrm_nlcd, width = 11, height = 8, dpi = 150)
ggsave("output/cross_ltrm_evt.png", p_cross_ltrm_evt,  width = 11, height = 8, dpi = 150)
ggsave("output/cross_nlcd_evt.png", p_cross_nlcd_evt,  width = 11, height = 8, dpi = 150)
ggsave("output/kw_effect_size.png", p_kw_effect,    width = 7,  height = 5,  dpi = 150)

# Combined inundation panels
(p_ltrm_inund | p_nlcd_inund) / p_evt_inund
ggsave("output/inund_all_systems.png", width = 14, height = 12, dpi = 150)

# Combined DTW panels
(p_ltrm_dtw | p_nlcd_dtw) / p_evt_dtw
ggsave("output/dtw_all_systems.png", width = 14, height = 12, dpi = 150)

# =============================================================================
# PART 10: Export tables
# =============================================================================

write_csv(inundation,       "output/ltrm_inundation.csv")
write_csv(ltrm_dtw_summary, "output/ltrm_dtw_summary.csv")
write_csv(model_comparison, "output/model_comparison_inund_dtw.csv")
write_csv(kw_results,       "output/kw_effect_sizes.csv")
write_csv(cross_ltrm_nlcd,  "output/cross_ltrm_nlcd.csv")
write_csv(cross_ltrm_evt,   "output/cross_ltrm_evt.csv")
write_csv(cross_nlcd_evt,   "output/cross_nlcd_evt.csv")

cat("\n=== Done ===\n")
cat("Cache files in:", cache_dir, "\n")
cat("Outputs in: output/\n")
cat("\nTo re-run any section without recomputing expensive steps,\n")
cat("simply source this script — cached .rds and .tif files will be loaded.\n")
cat("To force recomputation of a step, delete its file in", cache_dir, "\n")