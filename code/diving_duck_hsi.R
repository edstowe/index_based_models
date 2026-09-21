library(tidyverse)
library(sf)

#Calculate diving duck HSI score at ~12,000 locations
#Assess overall patterns of HSI scores
#Are there informative scores?
  #Examine:
  #Correlation of scores
#Sensitivity of models outcomes to different variables
#Presence of redundant variables
#Alternative HSI with reduced variable set (e.g., without invariant data)
#Explore spatial proxies

raw_sav<- read_csv("data/ltrm_veg_srs_data_0611093645/ltrm_vegsrs_data.csv")

aqa_dat <- st_read("data/aqa/aqa_2010_lvl3_011918.shp") %>%
  mutate(across(where(is.numeric), ~ na_if(.x, -9999))) %>%
  select(uniq_id, Acres)

# Filter to recent years
sav_sub1 <- raw_sav %>%
  mutate(date = as.Date(DATE, "%m/%d/%Y"),
         YEAR = format(date, "%Y")) %>%
  relocate(date, .after = DATE) %>%
filter(YEAR >= 2010 & YEAR <= 2019)

# Make spatial dataframe from survey data
sav_sf <- sav_sub1 %>%
  select(BARCODE, EAST_U, NORTH_U) %>%
  distinct() %>%
  st_as_sf(coords = c("EAST_U", "NORTH_U"), crs = st_crs(26915))

# Intersect survey and aquatic area data
sav_aqa <- st_intersection(sav_sf, aqa_dat) %>%
  st_drop_geometry() %>%
  rename(aqa_uniq_id = uniq_id)

# Filter to just ones with size
sav_sub <- left_join(sav_sub1, sav_aqa) %>%
  filter(!is.na(Acres))

# 1. Size
sav_size <- sav_sub %>% select(BARCODE, Acres)

# 2. Depth
sav_depth <- sav_sub %>%
  group_by(BARCODE) %>%
  summarize(avg_depth = mean(c(DEPTH1, DEPTH2, DEPTH3, DEPTH4, DEPTH5, DEPTH6)))

# 3. Percent SAV cover
sav_rake_cov <- sav_sub %>%
  group_by(BARCODE) %>%
  summarize(avg_density = mean(c(DENSITY1, DENSITY2, DENSITY3, DENSITY4, DENSITY5, DENSITY6)))

# 4. Key SAV species: VAAM3 = wild celery,  POPE6 = sago pondweed, 
#NLPW = narrow-leaved pondweeds, POCR3 = curly pondweed, PONO2 = longleaf pondweed, 
#POZO = flatstem pondweed, ZAPA = horned pondweed 

barcode_df <- select(sav_sub, BARCODE) %>% distinct()

sav_food_cov <- sav_sub %>% 
  select(BARCODE, SPPCD, RAKE1:RAKE6) %>%
  filter(SPPCD %in% c("VAAM3",  "POPE6", 
                      "NLPW", "POCR3", "PONO2", 
                      "POZO", "ZAPA")) %>%
  group_by(BARCODE) %>%
  summarize(RAKE1 = sum(RAKE1),
            RAKE2 = sum(RAKE2),
            RAKE3 = sum(RAKE3),
            RAKE4 = sum(RAKE4),
            RAKE5 = sum(RAKE5),
            RAKE6 = sum(RAKE6),
            n_sp = length(unique(SPPCD))) %>%
  rowwise() %>%
  mutate(sav_food_cov = median(c_across(RAKE1:RAKE6), na.rm = TRUE)) %>%
  ungroup()

sav_food_cov_all <- left_join(barcode_df, sav_food_cov)

# 5. Percent emergent veg. cover
emerg_cover <- sav_sub %>%
  select(BARCODE, COV_E) %>%
  distinct()
  
# 6. Key emergent species: SARI = stiff arrowhead; SCVA = softstem bulrush; ZIAQ = wild rice
emerg_food_cov <- sav_sub %>% 
  select(BARCODE, SPPCD, COVSPP) %>%
  filter(SPPCD %in% c("SARI", "SCVA", "ZIAQ")) %>%
  group_by(BARCODE) %>%
  summarize(COVSPP_key = sum(COVSPP),
            n_sp = length(unique(SPPCD)))

# 7. Invert. pops

# 8.  Disturbance (i.e., hunting)


#Data explore
test <- select(sav_sub, BARCODE, SPPCD, DENSITY1:DENSITY6, VISUAL1:VISUAL6, RAKE1:RAKE6) %>%
  + filter(SPPCD == "SARI")
test2 <- select(sav_sub, BARCODE, SPPCD, DENSITY1:DENSITY6, VISUAL1:VISUAL6, RAKE1:RAKE6) %>% filter(BARCODE %in% test$BARCODE)
test3 <- test2 %>% group_by(BARCODE) %>% summarize(DENSITY1 = length(unique(DENSITY1)),
                                                   vis = length(unique(VISUAL1)),
                                                   rake = length(unique(RAKE1)))

test4 <- select(sav_sub, BARCODE, SPPCD, COVSPP, DENSITY1:DENSITY6, VISUAL1:VISUAL6, RAKE1:RAKE6)

# DENSITY == total SAV on rake track (invariant across rows)
# VISUAL == presence/abs for species in vis. sample  (1,0)
# RAKE == density each species on rake



table(sav_sub$SPPCD) %>% sort(decreasing = T)

sav_sp <- sav_sub %>%
  select(SPPCD, VEG_S) %>%
  distinct() %>% 
  group_by(SPPCD) %>%
  filter(n()==1,
         VEG_S == "S")

e_sp <- sav_sub %>%
  select(SPPCD, VEG_E) %>%
  distinct() %>% 
  group_by(SPPCD) %>%
  filter(n()==1,
         VEG_E == "E")







# SAV species per code from delaney and larson
test <- sav_sub %>%
  select(SPPCD, COV_E) %>%
  distinct() %>%
  filter(COV_E == 0,
         SPPCD %in% c("CEDE4", "CHAR", "ELCA7", "MYSI", "MYSP2", "NAFL", "NAGU", "NAMI", "POCR3", "NLPW",
         "PONO2", "POPE6", "PORI2", "POZO", "RALO2", "UTMA", "VAAM3", "ZAPA", "ZODU"),
         depth = mean(DEPTH1, DEPTH2, DEPTH3, DEPTH4, DEPTH5))


# Size of water body

# Water depth

# Percent sub. cover

# Species sub cover



# Percent emerg cover



# 



# Hunting polygons Steve WInter


# Each survey:

# Emergent or not
# How determine emergent?
E = emergent present
BUT IT APPEARS A SITE IS EITHER ALL E (OR ALL S) OR NONE
HOW SEPARATE BY SPECIES? ARE ANY BOTH E AND S

# Which species
# Other pondweeds? 4 OTHERS:NARROW LEAVED, CURLY, LONGLEAF, FLATSTEM, HORNED

# SAV or not
# How determine SAV?
S = Submersed present
WHAT ABOUT COVER S?

# How define aquatic beds?
# Percent cover, presence of key species

# How define backwater
# Strata or aquatic areas?

# OTHER THINGS TO FILTER BY
# PROJECT CODE?
# Years? Season?
# USE PLANT DENSITY (RAKE GRABS?)
  
  # OTHER DIVING DUCK VARIABLES
  # SIZE OF WATERBODY
  # WATER DEPTH (PERCENT AREA 18 INCHES TO FIVE FEET)
  # KEY GROUPS OF INVERTS (SPHAERIID, GASTROPOD, HEXAGENIA, AMPHIDOA, CHIRONMID)
  # DISTURBANCE (HUNTING ACCESS)





