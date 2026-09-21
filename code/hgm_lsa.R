library(tidyverse)
library(sf)
library(leaflet)

# UTM Zone 15, NAD 27
# Transform to Lat Lon
lsa_raw <- sf::st_read(dsn = "umr_hgm/umrmvp_lsa/geomorph27.shp", crs = 26715) %>%
  st_transform(crs = 4326)

length(unique(lsa_raw$LSA_FULL))
length(unique(lsa_raw$LSA_NAME))

# Define new categories (fewer than original 44)
classes_lsa <- st_drop_geometry(lsa_raw) %>% 
  select(LSA_FULL, LSA_NAME) %>% 
  distinct()%>%
  mutate(
    lsa_cat = case_when(
      
      # ------------------------------------------------------------
      # 1. MAIN CHANNEL
      # ------------------------------------------------------------
      LSA_NAME %in% c(
        "MCMR", "MCXR", "MCCR", "MCZR", "MCRR", "MCIR", "MCWR",
        "MLB", "MVU", "MVM", "MVL", "MVC"
      ) ~ "Main Channel",
      
      # ------------------------------------------------------------
      # 2. MAIN-CHANNEL ISLAND
      # ------------------------------------------------------------
      LSA_FULL == "MAIN CHANNEL ISLAND" |
        LSA_FULL == "MAIN CHANNEL ISLAND WISCONSIN RIVER" ~
        "Main-Channel Island",
      
      # ------------------------------------------------------------
      # 3. INACTIVE/MINOR CHANNEL
      # ------------------------------------------------------------
      LSA_NAME %in% c("IVU", "ILB", "IVL", "IC") ~
        "Inactive/Minor Channel",
      
      # ------------------------------------------------------------
      # 4. MINOR-CHANNEL ISLAND
      # ------------------------------------------------------------
      LSA_FULL == "MINOR CHANNEL ISLAND" ~
        "Minor-Channel Island",
      
      # ------------------------------------------------------------
      # 5. MAIN-CHANNEL LAKE/BACKWATER
      # ------------------------------------------------------------
      LSA_FULL == "MAIN CHANNEL VERTICAL ACCRETION LAKE" ~
        "Main-Channel Lake/Backwater",
      
      # ------------------------------------------------------------
      # 6. MINOR-CHANNEL LAKE/BACKWATER
      # ------------------------------------------------------------
      LSA_FULL == "INACTIVE/MINOR CHANNEL VERTICAL ACCRETION LAKE" ~
        "Minor-Channel Lake/Backwater",
      
      # ------------------------------------------------------------
      # 7. OPEN WATER
      # ------------------------------------------------------------
      LSA_NAME == "W" ~
        "Open Water",
      
      # ------------------------------------------------------------
      # 8. TRIBUTARY CHANNEL
      # ------------------------------------------------------------
      LSA_FULL == "TRIBUTARY CHANNEL" |
        LSA_FULL == "VERMILLION RIVER AND ASSOCIATED LAKES" |
        LSA_NAME == "TC" ~
        "Tributary Channel",
      
      # ------------------------------------------------------------
      # 9. TRIBUTARY FLOODPLAIN
      # ------------------------------------------------------------
      LSA_FULL == "TRIBUTARY FLOODPLAIN" ~
        "Tributary Floodplain",
      
      # ------------------------------------------------------------
      # 10. TRIBUTARY MEANDER BELT
      # ------------------------------------------------------------
      LSA_FULL == "TRIBUTARY MEANDER BELT" ~
        "Tributary Meander Belt",
      
      # ------------------------------------------------------------
      # 11. TRIBUTARY FAN/DELTA
      # ------------------------------------------------------------
      LSA_FULL %in% c(
        "TRIBUTARY ALLUVIAL FAN",
        "TRIBUTARY FAN/DELTA"
      ) ~
        "Tributary Fan/Delta",
      
      # ------------------------------------------------------------
      # 12. TRIBUTARY MARSH
      # ------------------------------------------------------------
      LSA_FULL == "TRIBUTARY MARSH" ~
        "Tributary Marsh",
      
      # ------------------------------------------------------------
      # 13. STREAM SCARP
      # ------------------------------------------------------------
      LSA_FULL %in% c(
        "TRIBUTARY STREAM SCARP",
        "GLACIAL STREAM SCARP"
      ) ~
        "Stream Scarp",
      
      # ------------------------------------------------------------
      # 14. GLACIAL TERRACE
      # ------------------------------------------------------------
      LSA_NAME %in% c("GTM", "GTH", "GTL") ~
        "Glacial Terrace",
      
      # ------------------------------------------------------------
      # 15. VALLEY-SIDE / UPLAND
      # ------------------------------------------------------------
      LSA_FULL %in% c(
        "VALLEY SIDE COLLUVIAL SLOPES",
        "UPLAND HILLTOPS"
      ) ~
        "Valley-Side/Upland",
      
      # ------------------------------------------------------------
      # 16. LACUSTRINE SHORELINE
      # ------------------------------------------------------------
      LSA_NAME %in% c("LSB", "LSC") ~
        "Lacustrine Shoreline",
      
      # ------------------------------------------------------------
      # 17. EOLIAN
      # ------------------------------------------------------------
      LSA_FULL == "EOLIAN DUNES OVER TERRACE" ~
        "Eolian Dunes",
      
      # ------------------------------------------------------------
      # 18. MODIFIED LAND
      # ------------------------------------------------------------
      LSA_FULL == "MADE/MODIFIED LAND" ~
        "Modified Land",
      
      # ------------------------------------------------------------
      # ABANDONED GLACIAL CHANNEL
      # ------------------------------------------------------------
      LSA_FULL == "GLACIAL STREAM CHANNEL - ABANDONED" ~
        "Glacial Terrace",
      
      # Missing
      is.na(LSA_NAME) ~ "Unknown",
      
      # Anything not otherwise classified
      TRUE ~ "Other"
    )
  ) %>%
  mutate(
    habitat_cat = case_when(
      
      # Modern Channel
      LSA_NAME %in% c("MCMR") ~
        "Modern Channel",
      
      # Modern Backwater
      LSA_NAME %in% c("MVC", "W") ~
        "Modern Backwater",
      
      # Active Floodplain - Wet
      LSA_NAME %in% c("MVM", "MVU", "TF", "TFD") ~
        "Active Floodplain - Wet",
      
      # Active Floodplain - Dry
      LSA_NAME %in% c("MLB", "MCI", "TY", "TAF") ~
        "Active Floodplain - Dry",
      
      # Paleo-Floodplain - Wet
      LSA_NAME %in% c("IVM", "IVS", "IVU", "IC") ~
        "Paleo-Floodplain - Wet",
      
      # Paleo-Floodplain - Dry
      LSA_NAME %in% c("ILB", "ICI") ~
        "Paleo-Floodplain - Dry",
      
      # Natural Levees
      LSA_NAME %in% c("UH", "MVL") ~
        "Natural Levees",
      
      # Colluvial Slope
      LSA_NAME %in% c("TS", "VCS") ~
        "Colluvial Slope",
      
      # Glacial Terrace
      LSA_NAME %in% c(
        "EDT", "GSC", "GSS", "GTH", "GTL", "GTM", "TVT"
      ) ~
        "Glacial Terrace",
      
      # Missing
      is.na(LSA_NAME) ~ "Unknown",
      
      # Anything not represented in the table
      TRUE ~ "Other"
    )
  )

# Full dataset
lsa_cats <- left_join(lsa_raw, classes_lsa) 

# Areas 
lsa_sums <- st_drop_geometry(lsa_cats) %>%
  group_by(lsa_cat) %>%
  summarize(area = sum(AREA)) %>%
  arrange(desc(area))

# Count
lsa_count <- st_drop_geometry(lsa_cats) %>%
  group_by(lsa_cat) %>%
  summarize(count = n()) %>%
  arrange(desc(count))

# Plot
#Using st_crop with a named numeric vector (xmin, ymin, xmax, ymax)
#bbox <- c(xmin = -92.78, ymin = 44.64, xmax = -92.65, ymax = 44.68)
#lsa_sub <- st_crop(lsa_cats, bbox)

# Create a color palette function
pal <- colorFactor(
  palette = "viridis",          
  domain = lsa_cats$habitat_cat
)

# 5. Build and render the leaflet map
leaflet(data = lsa_cats) %>%
  addTiles() %>%                # Adds the default OpenStreetMap background
  addPolygons(
    fillColor = ~pal(habitat_cat), # Dynamically colors the polygons
    weight = 1,                 # Border thickness
    opacity = 1,                # Border opacity
    color = "white",            # Border color
    fillOpacity = 0.7,          # Transparency of polygon fill
    highlightOptions = highlightOptions(
      weight = 3,
      color = "#666",
      fillOpacity = 0.9,
      bringToFront = TRUE
    ),
    label = ~as.character(habitat_cat) # Shows value on hover
  ) %>%
  leaflet::addCircles(data = forest_sf) %>%
  addLegend(
    pal = pal, 
    values = ~habitat_cat, 
    opacity = 0.7, 
    title = "Legend Title", 
    position = "bottomright"
  )



# Transform to Lat Lon
audubon_raw <- sf::st_read(dsn = "audubon/AudubonUMR_BFS_Points_2014-2026/AudubonUMR_BFS_Points_20142026.shp") #%>%
  st_transform(crs = 4326)

USACE_forest_data <- read_csv("winona_fixed_2025.csv")

forest_sf <- st_as_sf(USACE_forest_data, coords = c("EASTING", "NORTHING"))

forest_sf <- st_read("Winona_fixed_2025/Winona_fixed_2025.shp") %>% 
  st_transform(crs = 4326)



# Plot subset
# Plot
#Using st_crop with a named numeric vector (xmin, ymin, xmax, ymax)
bbox <- c(xmin = -92.78, ymin = 44.64, xmax = -92.65, ymax = 44.68)
lsa_sub <- st_crop(lsa_cats, bbox) 

lsa_sub$area <- st_area(lsa_sub)

library(units)

lsa_sub$area_acres <- set_units(lsa_sub$area, "acre")

sub_areas <- lsa_sub %>%
  group_by(habitat_cat)%>% 
  summarize(area = sum(area_acres)) %>%
  mutate(habitat_cat = fct_reorder(habitat_cat, area)) 

lsa_sub1 <- lsa_sub %>%
  mutate(habitat_cat = factor(habitat_cat, levels = levels(sub_areas$habitat_cat)))

bar_plot <- ggplot(sub_areas)+
  geom_col(aes(x = area, y = habitat_cat, fill = habitat_cat)) +
  scale_fill_viridis_d()

map_plot <- ggplot(lsa_sub1)+
  geom_sf(aes(fill = habitat_cat)) +
  scale_fill_viridis_d()

#library(patchwork)
map_plot + bar_plot + patchwork::plot_layout(guides = "collect")


library(tidyverse)
library(sf)
library(units)
library(patchwork)

# --- Processing remains identical ---
bbox <- c(xmin = -92.78, ymin = 44.64, xmax = -92.65, ymax = 44.68)
lsa_sub <- st_crop(lsa_cats, bbox) 
lsa_sub$area <- st_area(lsa_sub)
lsa_sub$area_acres <- set_units(lsa_sub$area, "acre")

sub_areas <- lsa_sub %>%
  group_by(habitat_cat)%>% 
  summarize(area = sum(area_acres)) %>%
  mutate(habitat_cat = fct_reorder(habitat_cat, area)) 

ordered_levels <- levels(sub_areas$habitat_cat)

sub_areas <- sub_areas %>%
  mutate(habitat_cat = factor(habitat_cat, levels = ordered_levels))

lsa_sub1 <- lsa_sub %>%
  mutate(habitat_cat = factor(habitat_cat, levels = ordered_levels))

# --- Plotting with the Glyphs Aligned ---
bar_plot <- ggplot(sub_areas)+
  geom_col(aes(x = area, y = habitat_cat, fill = habitat_cat), key_glyph = "rect") +
  scale_fill_viridis_d() 

map_plot <- ggplot(lsa_sub1)+
  # FIX: Change the glyph here to match geom_col's rectangular key
  geom_sf(aes(fill = habitat_cat), key_glyph = "rect") + 
  scale_fill_viridis_d()

# This will now perfectly blend into a single, unified legend
(map_plot + bar_plot) + 
  patchwork::plot_layout(guides = "collect") & 
  # This moves the final collected legend to the bottom
  theme(legend.position = "bottom")
