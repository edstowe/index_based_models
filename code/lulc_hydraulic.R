# =============================================================================
# LULC vs. Hydraulic Variables — Pools 11–19
# Four flow scenarios: baseflow (60k), 5-yr (167k), 10-yr (190k), 25-yr (241k)
# Hydraulic variables: depth, velocity, froude
# Includes: inundation frequency analysis by land cover class
# =============================================================================

library(tidyverse)
library(sf)
library(terra)
library(ggdist)
library(patchwork)
library(dunn.test)

# =============================================================================
# USER SETTINGS — change here only
# =============================================================================

# Land cover class column — update if using a different classification
lc_col <- "CLASS_15_N"

# Pools to load
pools <- 11:19

# Raster aggregation factor (10 = 2m → 20m; reduces extraction time ~100x)
agg_factor <- 10

# Minimum fraction of polygon with valid hydraulic data to retain
min_coverage <- 0.10

# =============================================================================
# SECTION 1: Raster paths
# =============================================================================

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

# =============================================================================
# SECTION 2: Load and aggregate rasters (no reprojection)
# =============================================================================
# Vectors will be reprojected to raster CRS instead — much faster.
# Aggregate from 2m to ~20m before extraction to reduce memory/time.

load_and_agg <- function(path, label, fact = agg_factor) {
  cat("  Loading and aggregating", label, "...\n")
  r     <- rast(path)
  r_agg <- aggregate(r, fact = fact, fun = "mean", na.rm = TRUE)
  rm(r); gc()
  r_agg
}

cat("Loading baseflow rasters...\n")
bf_stack <- c(
  load_and_agg(rast_paths$bf$depth,    "BF depth"),
  load_and_agg(rast_paths$bf$velocity, "BF velocity"),
  load_and_agg(rast_paths$bf$froude,   "BF froude")
)
names(bf_stack) <- c("depth", "velocity", "froude")

cat("Loading 5-yr rasters...\n")
fyr_stack <- c(
  load_and_agg(rast_paths$fyr$depth,    "5YR depth"),
  load_and_agg(rast_paths$fyr$velocity, "5YR velocity"),
  load_and_agg(rast_paths$fyr$froude,   "5YR froude")
)
names(fyr_stack) <- c("depth", "velocity", "froude")

cat("Loading 10-yr rasters...\n")
tyr_stack <- c(
  load_and_agg(rast_paths$tyr$depth,    "10YR depth"),
  load_and_agg(rast_paths$tyr$velocity, "10YR velocity"),
  load_and_agg(rast_paths$tyr$froude,   "10YR froude")
)
names(tyr_stack) <- c("depth", "velocity", "froude")

cat("Loading 25-yr rasters...\n")
tyfr_stack <- c(
  load_and_agg(rast_paths$tyfr$depth,    "25YR depth"),
  load_and_agg(rast_paths$tyfr$velocity, "25YR velocity"),
  load_and_agg(rast_paths$tyfr$froude,   "25YR froude")
)
names(tyfr_stack) <- c("depth", "velocity", "froude")

# Grab raster CRS from one raster — all share the same CRS
rast_crs <- crs(bf_stack)
cat("\nRaster CRS:\n"); print(crs(bf_stack, describe = TRUE))

# =============================================================================
# SECTION 3: Load LULC shapefiles
# =============================================================================
# Reproject each pool to raster CRS BEFORE binding to avoid CRS mismatch error

cat("\nLoading LULC shapefiles...\n")

lulc_raw <- map(pools, function(p) {
  pool_str <- sprintf("%02d", p)
  path     <- sprintf("data/lulc/2020_LCU_Pool_%s/lc_2020_p%s.shp",
                      pool_str, pool_str)
  if (!file.exists(path)) {
    cat("  File not found:", path, "\n"); return(NULL)
  }
  cat("  Loading pool", p, "...\n")
  st_read(path, quiet = TRUE) %>%
    st_transform(crs = rast_crs) %>%   # reproject before binding
    mutate(pool = p)
}) %>%
  compact() %>%
  bind_rows()

cat("\nLULC columns:", names(lulc_raw), "\n")
cat("Land cover classes (", lc_col, "):\n")
print(table(lulc_raw[[lc_col]]))

# =============================================================================
# SECTION 4: Clip LULC to raster extent using coordinate filter
# =============================================================================
# Avoids CRS object mismatch from st_bbox; adds poly_id for later pivoting

rast_xmin <- ext(bf_stack)$xmin
rast_xmax <- ext(bf_stack)$xmax
rast_ymin <- ext(bf_stack)$ymin
rast_ymax <- ext(bf_stack)$ymax

lulc_centroids_xy <- st_coordinates(st_centroid(lulc_raw$geometry))

lulc_clip <- lulc_raw %>%
  mutate(
    cx = lulc_centroids_xy[, 1],
    cy = lulc_centroids_xy[, 2]
  ) %>%
  filter(cx >= rast_xmin, cx <= rast_xmax,
         cy >= rast_ymin, cy <= rast_ymax) %>%
  select(-cx, -cy) %>%
  mutate(
    poly_id  = row_number(),            # unique ID for pivot_wider later
    lc_class = .data[[lc_col]]
  )

cat("\nPolygons overlapping raster extent:", nrow(lulc_clip), "\n")
cat("By pool:\n"); print(table(lulc_clip$pool))
cat("By land cover class:\n"); print(table(lulc_clip$lc_class))

# =============================================================================
# SECTION 5: Extraction function
# =============================================================================

extract_by_pool <- function(lulc_pool, stack, scenario_label, pool_id) {
  cat("  Extracting pool", pool_id, scenario_label, "...\n")
  v  <- vect(lulc_pool)
  ex <- terra::extract(stack, v, fun = mean, na.rm = TRUE, bind = FALSE)
  
  lulc_pool %>%
    st_drop_geometry() %>%
    mutate(
      depth_mean    = ex[, 2],
      velocity_mean = ex[, 3],
      froude_mean   = ex[, 4],
      # pct_covered: fraction of polygon cells with valid depth
      pct_covered   = terra::extract(
        bf_stack[["depth"]], v,
        fun = function(x, ...) mean(!is.na(x)),
        bind = FALSE)[, 2],
      scenario      = scenario_label
    )
}

# =============================================================================
# SECTION 6: Run extraction — all scenarios, all pools
# =============================================================================

lulc_hydro <- map_dfr(unique(lulc_clip$pool), function(p) {
  pool_lulc <- filter(lulc_clip, pool == p)
  bind_rows(
    extract_by_pool(pool_lulc, bf_stack,   scenario_labels["bf"],   p),
    extract_by_pool(pool_lulc, fyr_stack,  scenario_labels["fyr"],  p),
    extract_by_pool(pool_lulc, tyr_stack,  scenario_labels["tyr"],  p),
    extract_by_pool(pool_lulc, tyfr_stack, scenario_labels["tyfr"], p)
  )
})

gc()

# =============================================================================
# SECTION 7: Coverage summary and filtering
# =============================================================================

cat("\nCoverage by land cover class (baseflow):\n")
lulc_hydro %>%
  filter(scenario == scenario_labels["bf"]) %>%
  group_by(lc_class) %>%
  summarize(
    n_polygons   = n(),
    mean_pct_cov = round(mean(pct_covered, na.rm = TRUE) * 100, 1),
    n_any_cov    = sum(pct_covered > 0, na.rm = TRUE),
    .groups      = "drop"
  ) %>%
  arrange(desc(mean_pct_cov)) %>%
  print()

# Filter to polygons with meaningful coverage at baseflow
# (if covered at baseflow, higher flows will also cover them)
covered_polys <- lulc_hydro %>%
  filter(scenario == scenario_labels["bf"],
         pct_covered >= min_coverage) %>%
  pull(poly_id)

lulc_hydro_clean <- lulc_hydro %>%
  filter(poly_id %in% covered_polys) %>%
  mutate(scenario = factor(scenario, levels = unname(scenario_labels)))

cat("\nPolygons retained (>=", min_coverage * 100, "% coverage):",
    n_distinct(lulc_hydro_clean$poly_id), "\n")
cat("By class:\n")
print(table(filter(lulc_hydro_clean,
                   scenario == scenario_labels["bf"])$lc_class))

# Order land cover classes by median baseflow depth
class_order <- lulc_hydro_clean %>%
  filter(scenario == scenario_labels["bf"]) %>%
  group_by(lc_class) %>%
  summarize(med = median(depth_mean, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(lc_class)

lulc_hydro_clean <- mutate(lulc_hydro_clean,
                           lc_class = factor(lc_class, levels = class_order))

# =============================================================================
# SECTION 8: Summary statistics
# =============================================================================

lulc_summary <- lulc_hydro_clean %>%
  group_by(lc_class, scenario) %>%
  summarize(
    n               = n(),
    depth_median    = round(median(depth_mean,    na.rm = TRUE), 2),
    depth_iqr       = round(IQR(depth_mean,       na.rm = TRUE), 2),
    velocity_median = round(median(velocity_mean, na.rm = TRUE), 3),
    velocity_iqr    = round(IQR(velocity_mean,    na.rm = TRUE), 3),
    froude_median   = round(median(froude_mean,   na.rm = TRUE), 3),
    .groups         = "drop"
  )

cat("\nSummary by class and scenario:\n")
print(lulc_summary)

# =============================================================================
# SECTION 9: Inundation frequency analysis
# =============================================================================
# For each polygon, determine the lowest flow scenario at which it is inundated
# (pct_covered >= min_coverage). This gives a coarse inundation frequency class.
# Scenarios are ordered by increasing discharge:
#   baseflow < 5-yr < 10-yr < 25-yr
# A polygon inundated only at 25-yr is rarely flooded; one inundated at
# baseflow is permanently or frequently inundated.

inund_levels <- c(
  "Permanent / baseflow",    # inundated at baseflow
  "5-yr event",              # first inundated at 5-yr
  "10-yr event",             # first inundated at 10-yr
  "25-yr event",             # first inundated at 25-yr
  "Not inundated (>25-yr)"   # never inundated across all scenarios
)

inundation <- lulc_hydro %>%
  filter(poly_id %in% lulc_clip$poly_id) %>%   # all polygons, not just covered
  select(poly_id, lc_class, pool, scenario, pct_covered) %>%
  mutate(
    inundated = pct_covered >= min_coverage,
    scenario  = factor(scenario, levels = unname(scenario_labels))
  ) %>%
  group_by(poly_id, lc_class, pool) %>%
  summarize(
    inund_bf   = any(inundated[scenario == scenario_labels["bf"]],   na.rm = TRUE),
    inund_fyr  = any(inundated[scenario == scenario_labels["fyr"]],  na.rm = TRUE),
    inund_tyr  = any(inundated[scenario == scenario_labels["tyr"]],  na.rm = TRUE),
    inund_tyfr = any(inundated[scenario == scenario_labels["tyfr"]], na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(
    inund_class = case_when(
      inund_bf              ~ "Permanent / baseflow",
      !inund_bf & inund_fyr ~ "5-yr event",
      !inund_fyr & inund_tyr  ~ "10-yr event",
      !inund_tyr & inund_tyfr ~ "25-yr event",
      TRUE                  ~ "Not inundated (>25-yr)"
    ),
    inund_class = factor(inund_class, levels = inund_levels),
    lc_class    = factor(lc_class,    levels = class_order)
  )

cat("\nInundation frequency by land cover class:\n")
print(table(inundation$lc_class, inundation$inund_class))

# Proportion of each land cover class in each inundation frequency category
inund_props <- inundation %>%
  group_by(lc_class, inund_class) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(lc_class) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

# =============================================================================
# SECTION 10: Plots
# =============================================================================

scenario_colors <- c(
  "Baseflow (60k cfs)"  = "#2A7D6F",
  "5-yr (167k cfs)"     = "#D4850A",
  "10-yr (190k cfs)"    = "#7B3F9E",
  "25-yr (241k cfs)"    = "#C0392B"
)

inund_colors <- c(
  "Permanent / baseflow"  = "#1B6B3A",
  "5-yr event"            = "#2A7D6F",
  "10-yr event"           = "#D4850A",
  "25-yr event"           = "#C0392B",
  "Not inundated (>25-yr)"= "#B0BAC5"
)

# ── Plot 1: Depth by land cover class × scenario ─────────────────────────────

p_depth <- ggplot(lulc_hydro_clean,
                  aes(x = lc_class, y = depth_mean, fill = scenario)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.65,
               position = position_dodge(0.75)) +
  scale_fill_manual(values = scenario_colors, name = "Flow scenario") +
  scale_y_log10(labels = scales::label_number(accuracy = 0.1)) +
  coord_flip() +
  labs(title    = "Mean polygon depth by land cover class",
       subtitle = "Four flow scenarios  |  Log scale",
       x = NULL, y = "Mean depth (m, log scale)") +
  theme_minimal(base_size = 12) +
  theme(legend.position    = "bottom",
        panel.grid.major.y = element_blank())

# ── Plot 2: Velocity by land cover class × scenario ──────────────────────────

p_vel <- ggplot(lulc_hydro_clean,
                aes(x = lc_class, y = velocity_mean, fill = scenario)) +
  geom_boxplot(outliers = F, alpha = 0.8, width = 0.65,
               position = position_dodge(0.75)) +
  scale_fill_manual(values = scenario_colors, name = "Flow scenario") +
  scale_y_log10(labels = scales::label_number(accuracy = 0.001)) +
  coord_flip() +
  labs(title    = "Mean polygon velocity by land cover class",
       subtitle = "Four flow scenarios  |  Log scale",
       x = NULL, y = "Mean velocity (m/s, log scale)") +
  theme_minimal(base_size = 12) +
  theme(legend.position    = "bottom",
        panel.grid.major.y = element_blank())

# ── Plot 3: Depth increase from baseflow across scenarios ─────────────────────

p_delta <- lulc_hydro_clean %>%
  select(poly_id, lc_class, scenario, depth_mean) %>%
  group_by(lc_class, scenario) %>%
  mutate(within_group_id = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from  = scenario,
              values_from = depth_mean,
              id_cols     = c(poly_id, lc_class, within_group_id)) %>%
  pivot_longer(
    cols      = all_of(unname(scenario_labels[-1])),  # all except baseflow
    names_to  = "scenario",
    values_to = "depth_event"
  ) %>%
  rename(depth_bf = all_of(unname(scenario_labels["bf"]))) %>%
  filter(!is.na(depth_bf), !is.na(depth_event)) %>%
  mutate(
    depth_increase = depth_event - depth_bf,
    scenario       = factor(scenario, levels = unname(scenario_labels[-1]))
  ) %>%
  ggplot(aes(x = lc_class, y = depth_increase, fill = scenario)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.65,
               position = position_dodge(0.75)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(
    values = scenario_colors[unname(scenario_labels[-1])],
    name   = "Event scenario"
  ) +
  coord_flip() +
  labs(title    = "Depth increase above baseflow by event scenario",
       subtitle = "Positive = deeper at higher flow  |  Classes ordered by baseflow depth",
       x = NULL, y = "Δ depth above baseflow (m)") +
  theme_minimal(base_size = 12) +
  theme(legend.position    = "bottom",
        panel.grid.major.y = element_blank())

# ── Plot 4: Inundation frequency by land cover class (stacked bar) ────────────

p_inund_stack <- ggplot(inund_props,
                        aes(x = lc_class, y = pct, fill = inund_class)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = inund_colors, name = "Inundation frequency") +
  scale_y_continuous(labels = scales::label_number(suffix = "%")) +
  coord_flip() +
  labs(title    = "Inundation frequency class by land cover type",
       subtitle = paste0("Based on ", min_coverage * 100,
                         "% polygon coverage threshold across four flow scenarios"),
       x = NULL, y = "Proportion of polygons (%)") +
  theme_minimal(base_size = 12) +
  theme(legend.position    = "bottom",
        panel.grid.major.y = element_blank(),
        legend.text        = element_text(size = 8))

# ── Plot 5: Inundation frequency × depth (hydraulic profile at each freq) ─────

p_inund_depth <- lulc_hydro_clean %>%
  filter(scenario == scenario_labels["bf"]) %>%
  left_join(select(inundation, poly_id, inund_class), by = "poly_id") %>%
  filter(!is.na(inund_class)) %>%
  ggplot(aes(x = inund_class, y = depth_mean, fill = inund_class)) +
  stat_halfeye(adjust = 0.6, width = 0.5, .width = 0,
               justification = -0.2, point_colour = NA, alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6,
               color = "grey30") +
  scale_fill_manual(values = inund_colors, guide = "none") +
  scale_y_log10(labels = scales::label_number(accuracy = 0.1)) +
  coord_flip() +
  labs(title    = "Baseflow depth by inundation frequency class",
       subtitle = "Log scale  |  Shows hydraulic profile of each inundation class",
       x = NULL, y = "Mean baseflow depth (m, log scale)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

# ── Assemble and save ─────────────────────────────────────────────────────────

# Main hydraulic panels
(p_depth / p_vel) | p_delta
ggsave("output/lulc_hydraulics_by_scenario.png",
       width = 15, height = 10, dpi = 150)

# Inundation panels
p_inund_stack / p_inund_depth
ggsave("output/lulc_inundation_frequency.png",
       width = 10, height = 10, dpi = 150)

# =============================================================================
# SECTION 11: Statistical tests
# =============================================================================

for (sc in levels(lulc_hydro_clean$scenario)) {
  dat_sc <- filter(lulc_hydro_clean, scenario == sc, !is.na(depth_mean))
  cat("\n=== Kruskal-Wallis:", sc, "===\n")
  kw_d <- kruskal.test(depth_mean    ~ lc_class, data = dat_sc)
  kw_v <- kruskal.test(velocity_mean ~ lc_class, data = dat_sc)
  cat("  Depth    p =", signif(kw_d$p.value, 3), "\n")
  cat("  Velocity p =", signif(kw_v$p.value, 3), "\n")
}

# Inundation frequency × land cover: chi-square test
cat("\nChi-square: inundation class ~ land cover class\n")
print(chisq.test(table(inundation$lc_class, inundation$inund_class)))

# =============================================================================
# SECTION 12: Export
# =============================================================================

write_csv(lulc_hydro_clean, "output/lulc_hydraulic_extracted.csv")
write_csv(lulc_summary,     "output/lulc_hydraulic_summary.csv")
write_csv(inundation,       "output/lulc_inundation_frequency.csv")
write_csv(inund_props,      "output/lulc_inundation_proportions.csv")

cat("\nDone.\n")