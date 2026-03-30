
library(terra)
library(tidyterra)
library(geodata)
library(dplyr)
library(ggplot2)
library(landscapemetrics)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)




#SECTION 0: define crs

WGSCRS <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ATACRS <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"



?landscapemetrics
#----------------------------------------------------------------------------------------------------------------------------------

#SECTION 1: LOAD RASTERS

#present lichen rasters

p_ag <- rast("C:/Users/Student/OneDrive/Documents/EDINBURGH/Dissertation/Maxent_Outputs/AG/P_AG_JK/Acarospora_gwynnii_MaxEnt_layers.asc")
p_ud <- rast("C:/Users/Student/OneDrive/Documents/EDINBURGH/Dissertation/Maxent_Outputs/UD/P_UD_JK/Umbilicaria_decussata_MaxEnt_layers.asc")
p_ua <- rast("C:/Users/Student/OneDrive/Documents/EDINBURGH/Dissertation/Maxent_Outputs/UA/P_UA_JK/Usnea_antarctica_MaxEnt_layers.asc")

#future lichen rasters
f_ag <- rast("C:/Users/Student/OneDrive/Documents/EDINBURGH/Dissertation/Maxent_Outputs/AG/F_AG_JK/Acarospora_gwynnii_F_MaxEnt_layers.asc")
f_ud <- rast("C:/Users/Student/OneDrive/Documents/EDINBURGH/Dissertation/Maxent_Outputs/UD/F_UD_JK/Umbilicaria_decussata_F_MaxEnt_layers.asc")
f_ua <- rast("C:/Users/Student/OneDrive/Documents/EDINBURGH/Dissertation/Maxent_Outputs/UA/F_UA_JK/Usnea_antarctica_F_MaxEnt_layers.asc")

# 1. Ensure your original raster has the correct CRS
# Replace 'ATACRS' with the actual CRS string if not already defined

w <- p_ag

crs(w) <- ATACRS
  latlon <- project(w, WGSCRS, method = "bilinear")
  
  
  
#2: Create nice plots of SDMS

#Get the coastline data (scale 10 is high resolution)
coast <- ne_coastline(scale = "medium", returnclass = "sf")

coast <- st_make_valid(coast)

f <- ggplot() +
  # The Raster (MaxEnt)
  geom_spatraster(data = latlon) +
  
  # The Coastline
  geom_sf(data = coast, fill = NA, color = "black", linewidth = 0.4) +
  
  # The "Bending" projection
  coord_sf(crs = "+proj=laea +lat_0=-90 +lon_0=-65 +datum=WGS84", 
           xlim = c(-80, -50), ylim = c(-78, -60),
           expand = FALSE,
           default_crs = st_crs(4326)) +
   scale_fill_viridis_c(option = "viridis", na.value = "transparent") +
  #different type of colour system to try
  #scale_fill_gradientn(colors = c("#0000ff", "#0036ff", "#006cff", "#00a2ff", "#00d7ff", "#00fff1", "#00ffbb", "#00ff86", "#00ff50", "#00ff1a", "#1bff00", "#51ff00", "#87ff00", "#bcff00", "#f2ff00", "#ffd600", "#ffa100", "#ff6b00", "#ff3500", "#ff0000"), na.value = "transparent") +
  theme_minimal() +
  labs(title = "MaxEnt Present Suitability: AG",
       fill = "Probability")
plot(f)

ggsave(
  filename = "newP_A_gwynnii.pdf",
  plot = f,
  width = 10,
  height = 8,
  dpi = 600,
  bg = "white",
)

#--------------------------------------------------------------------------------------------------
# SECTION 2: converting to binary vector, so anaylsis can begin

#set threshold value. this now sets rasters to vector. from continous data ie likelihood - to binary data.
#thresholds from maxentresults excel ss (in maxent_output folder - sensitivity to cumulative threshold column)



binary_p_ud <- p_ud >= 0.498
binary_p_ua <- p_ua >= 0.5551
binary_p_ag <- p_ag >= 0.2105


binary_f_ud <- f_ud >= 0.498
binary_f_ua <- f_ua >= 0.5551
binary_f_ag <- f_ag >= 0.2105



check_landscape(binary_p_ud)

crs(binary_f_ud) <- ATACRS



# set species of interest binary_p_""

current_binary <- binary_p_ag#ua#ag
future_binary <- binary_f_ag#ua#ag

#sum of all ice free pixels for current and future
total_current_pixels <- global(current_binary, fun = "notNA")$notNA
total_future_pixels <- global(future_binary, fun = "notNA")$notNA
#percentage change
ice_change <- ((total_future_pixels - total_current_pixels)/ total_current_pixels) * 100
#7% increase in ice-free area ??? definitely wrong, and think im correct about resolution


#sum of binary TRUE pixels for present and future
current_sum <- global(current_binary, "sum", na.rm = TRUE)$sum
future_sum <- global(future_binary, "sum", na.rm = TRUE)$sum

#calculate percentage change
percentage_change <- ((future_sum - current_sum) / current_sum) * 100


#need to use these values to figure out percentage gain and percentage loss. however these figures dont line up due to their being a differnce in pixels.
#says to tell r to ignore na values as currently with crosstab its treating na to 0s or 1s as a shift. 
# this might not make sense though as the total from the changetable is lower (11018 than either the total_future pixels 13370 OR even current 12485 )
#so its missing / ignoring pixels for some reason. maybe need to figure out if theres loss between the total pixels. 
#aka the 7% is net not gross - this would be strange

#-------------------------------------------------------------------------------------------------------------------------------
#SECTION 3: GAIN, LOSS AND UNCHANGED TABLES

# understand habitat shift table
combined <- c(current_binary, future_binary)
names(combined) <- c("Current", "Future")

plot(combined)

#table with current, future and n columns. 
#0 = no presence. 1 = presence. n = number of pixels that this is true for. 
#0 - 0 = not found in either. #1 - 0 = habitat loss. 
#0 - 1 = habitat gain. 1 - 1 = persisten in habitat.
change_table <- crosstab(combined, long = TRUE)

print(change_table)



#ICE MASK, TO UNDERSTAND GAINS AND LOSS OF LICHENS RELATED TO ICE GAIN AND LOSS

#mask for ice cells in current sdm which are gained in the future sdm
gain_mask <- is.na(current_binary) & !is.na(future_binary)

#test to check same as gain in change table
total_new_land <- global(gain_mask, "sum", na.rm = TRUE)$sum
print(total_new_land)


f_ice_gain <- future_binary %>% 
  mask(gain_mask, maskvalues = 0, updatevalue = NA)

lichen_gained <- global(f_ice_gain == 1, "sum", na.rm = TRUE)$sum
lichen_never <- global(f_ice_gain == 0, "sum", na.rm = TRUE)$sum



#mask for pixels present in current and na in future
loss_mask <- !is.na(current_binary) & is.na(future_binary)
plot(loss_mask)

total_land_lost <- global(loss_mask, "sum", na.rm = TRUE)$sum


f_ice_inc <- current_binary %>% 
  mask(loss_mask, maskvalues = 0, updatevalue = NA)

lichen_lost <- global(f_ice_inc == 1, "sum", na.rm = TRUE)$sum
lichen_never2 <- global(f_ice_inc == 0, "sum", na.rm = TRUE)$sum


#Testing for NA cells in current vs future etc

current_binary_test <- current_binary
values(current_binary_test) <- ifelse(is.na(values(current_binary_test)), 0, 1)

future_binary_test <- future_binary
values(future_binary_test) <- ifelse(is.na(values(future_binary_test)), 0, 1)

combined_test <- c(current_binary_test, future_binary_test)
names(combined_test) <- c("Current", "Future")

change_table <- crosstab(combined_test, long = TRUE)
print(change_table)

# ud 
never = 2000 + 937 + 374
# ud 
gain = 669 + 1424
# ud 
loss = 1919 + 1093
# ud 
persist = 6430 

# ua 
never = 3863 + 1371 + 1181
# ua 
gain = 4575 + 990
# ua 
loss = 304 + 286
# ua 
persist = 2276

# ag 
never = 2772 + 1039 + 301
# ag 
gain = 83 + 1322
# ag 
loss = 4403 + 1166
# ag persist = 3760 

#I want to look at only the cells that are na in current and 1 in future.
#I want to look at only the cells that are 1 in current and na in future.

5192*899.46



#--------------------------------------------------------------------------------------------

# SECTION 3: EXPLORING PATCH DYNAMICS 

#pick targets
target_current <- binary_p_ag
target_future <- binary_f_ag

#  Check Landscape Integrity
# This verifies if patches are valid for analysis

crs(target_future) <- ATACRS

current_clean <- terra::classify(target_current, cbind(0, NA))
future_clean <- terra::classify(target_future, cbind(0, NA))


check_ta <- calculate_lsm(current_clean, what = "lsm_l_ta")
print(check_ta)


check_f <- calculate_lsm(future_clean, what = "lsm_l_ta")
print(check_f)

#ud current = 8492701
#ud future = 7666097

#ua c = 2577852
#ua f = 7052665

#ag c = 8391062
#ag f = 4645711



# Calculate Key Patch Metrics
# NP = Number of Patches
# MPS = Mean Patch Area
# ENN_MN = Mean Distance to Nearest Neighbor


metrics_to_calculate <- c("lsm_l_np", "lsm_l_area_mn", "lsm_l_enn_mn")

current_results <- calculate_lsm(current_clean, what = metrics_to_calculate)
future_results <- calculate_lsm(future_clean, what = metrics_to_calculate)

# 5. Compare the results
comparison_ud <- bind_rows(
  current_results %>% mutate(period = "Current"),
  future_results %>% mutate(period = "Future")
)

print(comparison_ud)

#--------------------------------------------------------------------------------------------

#SECTION 4: PERSITENCE MAPPING AND SHARED OVERLAP


persistence_ud <- (binary_p_ud == 1) & (binary_f_ud == 1)
persistence_ua <- (binary_p_ua == 1) & (binary_f_ua == 1)
persistence_ag <- (binary_p_ag == 1) & (binary_f_ag == 1)

shared_persistence <- persistence_ud + persistence_ua + persistence_ag
shared_future <- binary_f_ag + binary_f_ua + binary_f_ud
shared_current <- binary_p_ag + binary_p_ua + binary_p_ud

persistence_counts <- freq(shared_persistence)
print(persitence_counts)

future_counts <- freq(shared_future)
print(future_counts)

current_counts <- freq(shared_current)
print(current_counts)

# basic shared maps
png("Shared_persistence.png", width = 2000, height = 2000, res = 300)

plot(shared_persistence, 
     col = c("lightblue", "yellow", "green", "magenta"), 
     main = "Number of Persisting Lichen species")

dev.off()

#
png("Shared_future.png", width = 2000, height = 2000, res = 300)

plot(shared_future,
     col = c("lightblue", "yellow", "green", "magenta"), 
     main = "Number of future Lichen Species")

dev.off()

png("Shared_current.png", width = 2000, height = 2000, res = 300)

plot(shared_current,
     col = c("lightblue", "yellow", "green", "darkgreen"), 
     main = "Number of current Lichen Species")
dev.off()

# NICE SHARED MAPS

run <- shared_persistence

crs(run) <- ATACRS

shared_latlon <- project(run, WGSCRS)

shared_latlon <- as.factor(shared_latlon)

plot(shared_latlon)

#levels(ag_latlon_factor) <- data.frame(UD=c(1,2,3),
                                      # category=c("Unsuitable","Marginal",))

# 1. Get the coastline data (scale 10 is high resolution)
coast <- ne_coastline(scale = "medium", returnclass = "sf")

coast <- st_make_valid(coast)

shared <- ggplot() +
  geom_spatraster(data = shared_latlon) +
  geom_sf(data = coast, fill = NA, color = "black", linewidth = 0.4) +
  # The "Bending" projection
  coord_sf(crs = "+proj=laea +lat_0=-90 +lon_0=-65 +datum=WGS84", 
           xlim = c(-80, -50), ylim = c(-78, -60),
           expand = FALSE,
           default_crs = st_crs(4326)) +
  scale_fill_manual(values = c("0" = "#e1e1e1", "1" = "#FFFF00","2" = "#58bbed","3" = "#FF00FF"), 
                    na.translate = FALSE, 
                    name = "Number of lichens",
  ) +
  theme_minimal() +
  labs(title = "MaxEnt Persistence Suitability: Shared",
       fill = "Probability")
plot(shared)



ggsave(
  filename = "persistence_shared.pdf",
  plot = shared,
  width = 10,
  height = 8,
  dpi = 600,
  bg = "white",
)



### ------------------------------------------------------------------------------------------------------------
#SECTION 5: GAIN AND LOSS MAPS

## first need to na all pixels that stay zeros on the gain maps

## add the current and future together - make combined. 
### only select the zeroes - mask this onto current and future. this hides the zeroes 

## then do future - current. -1 = loss, 0 = unchanged, 1 = gain.

combined_ud <- binary_f_ud + binary_p_ud
combined_ud[combined_ud > 0] <- NA

plot(combined_ud)

mask_p_ud <- binary_p_ud %>% 
  mask(combined_ud, maskvalues = 0, updatevalue = NA)
mask_f_ud <- binary_f_ud %>% 
  mask(combined_ud, maskvalues = 0, updatevalue = NA)

gains_ud <- mask_f_ud - mask_p_ud
plot(gains_ud)

#repeat for ag
combined_ag <- binary_f_ag + binary_p_ag
combined_ag[combined_ag > 0] <- NA

plot(combined_ag)

mask_p_ag <- binary_p_ag %>% 
  mask(combined_ag, maskvalues = 0, updatevalue = NA)
mask_f_ag <- binary_f_ag %>% 
  mask(combined_ag, maskvalues = 0, updatevalue = NA)

plot(mask_p_ag)
plot(mask_f_ag)

gains_ag <- mask_f_ag - mask_p_ag
plot(gains_ag)

#repeat for ua
combined_ua <- binary_f_ua + binary_p_ua
combined_ua[combined_ua > 0] <- NA

plot(combined_ua)

mask_p_ua <- binary_p_ua %>% 
  mask(combined_ua, updatevalue = NA)
mask_f_ua <- binary_f_ua %>% 
  mask(combined_ua, maskvalues = 0, updatevalue = NA)

gains_ua <- mask_f_ua - mask_p_ua
plot(gains_ua)

#new ag
combined_ag <- binary_f_ag + binary_p_ag
combined_ag[combined_ag < 1] <- NA


mask_p_ag <- binary_p_ag %>% 
  mask(combined_ag)
mask_f_ag <- binary_f_ag %>% 
  mask(combined_ag)

gains_ag <- mask_f_ag - mask_p_ag
plot(gains_ag)

#new ua
combined_ua <- binary_f_ua + binary_p_ua
combined_ua[combined_ua < 1] <- NA


mask_p_ua <- binary_p_ua %>% 
  mask(combined_ua)
mask_f_ua <- binary_f_ua %>% 
  mask(combined_ua)

gains_ua <- mask_f_ua - mask_p_ua
plot(gains_ua)

#new ud
combined_ud <- binary_f_ud + binary_p_ud
combined_ud[combined_ud < 1] <- NA


mask_p_ud <- binary_p_ud %>% 
  mask(combined_ud)
mask_f_ud <- binary_f_ud %>% 
  mask(combined_ud)

gains_ud <- mask_f_ud - mask_p_ud
plot(gains_ud)

# nice maps for gains


#set to species of interest
gains <- gains_ud#ud#ua

#ensure in correct crs
crs(gains) <- ATACRS

#convert to lon,lat crs for plotting map
gains <- project(gains, WGSCRS) %>% 
  as.factor()

#check plot
plot(gains)


#nice map with coastline underlayed
coast <- ne_coastline(scale = "medium", returnclass = "sf")

coast <- st_make_valid(coast)

changes_heh <- ggplot() +
  geom_spatraster(data = gains) +
  geom_sf(data = coast, fill = NA, color = "black", linewidth = 0.4) +
  # The "Bending" projection
  coord_sf(crs = "+proj=laea +lat_0=-90 +lon_0=-65 +datum=WGS84", 
           xlim = c(-80, -50), ylim = c(-78, -60),
           expand = FALSE,
           default_crs = st_crs(4326)) +
  scale_fill_manual(values = c("-1" = "#ff1300", "0" = "#FFFF00","1" = "#43f531"), 
                    na.translate = FALSE, 
                    name = "Distribution shifts",
                    labels = c("Loss", "Unchanged", "Gain")
  ) +
  theme_minimal() +
  labs(title = "MaxEnt distribution shifts: ud",
       fill = "Probability")
plot(changes_heh)

#save file to species name
ggsave(
  filename = "gains_ud.pdf",
  plot = changes_heh,
  width = 10,
  height = 8,
  dpi = 600,
  bg = "white",
)

