

library(terra)
library(geodata)
library(dplyr)

# 0: Define projections

WGSCRS <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ATACRS <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"

#Steps 1-3 Loading and manipulating the three data sets

#--------------------------------------------------------------------------------------

## 1: convert BIOCLIM data to Antartic (ATA) Coordinate reference system (CRS)
# 1:1 LOAD ALL 19 bioclim VARIABLES

Bioclim_Global <- worldclim_global(var = "bio", res = 2.5, path = tempdir(), version = "2.1")
plot(Bioclim_Global[[1]], main = "Bio 1: Annual Mean Temp")  

# Check the names to see you have wc2.1_10m_bio_1 to bio_19
# broken print(names(Bioclim_Global ....

# 1:2 Define extent towards Antarctic peninsula (AP)
#(might have to swap to antarctic crs first and then figure out new extent coords in the reference system - this is because the ice-cover data is already in ATA reference system and so cannot apply this extent to it)
#currently having to swap ice data to WGCRS, apply this extent and then convert back to ATA
AP_Ext <- ext(x = c(-80, -55, -78, -60))

# 1:3: CROP ALL 19 LAYERS AT ONCE
AP_Bioclim <- crop(Bioclim_Global, AP_Ext)

# 1:4 Visualize just the first one (Annual Mean Temperature) to check
plot(AP_Bioclim[[1]], main = "Bio 1: Annual Mean Temp")

# 1:5 convert to CORRECT CRS (PROJECT)
ATA_Bioclim_Proj <- project(AP_Bioclim, ATACRS)

# 1:6 Plot the result (test with bioclim 1, temp)
plot(ATA_Bioclim_Proj[[1]], main = "Projected Bio 1")

# 1. Define the output file path and name
# Note: .tif is the standard format for geospatial data
output_path <- "ATA_Bioclim_Proj.tif"

# 2. Write the raster to disk
# we use overwrite=TRUE in case you run the code more than once
writeRaster(ATA_Bioclim_Proj, filename = output_path, overwrite = TRUE)

message("Success! Your clean variables are saved to: ", getwd())

ATA_Bioclim_Proj <- rast("ATA_Bioclim_Proj.tif")

#--------------------------------------------------------------------


#  2: Apply Ice mask
# 2:1 load ice data
Ice_Free_Area <- vect("C:/Users/Student/OneDrive/Documents/EDINBURGH/Dissertation/Modelling/Data/Penisula_IceCover/add_rock_outcrop_medium_res_polygon_v7_11.shp")
                      #C:/Users/Student/OneDrive - University of Edinburgh/Documents/EDINBURGH/Dissertation/Modelling/Data/Penisula_IceCover/add_rock_outcrop_medium_res_polygon_v7_11.shp) 
plot(Ice_Free_Area)

  
#This is a work around. will convert the extent to ATACRS later instead. but need to know coordinates for antartic peninsula in that reference system first.
#2.2 convert to wgscrs so its same as crs as ATA_ext
Ice_Free_WG <- project(Ice_Free_Area, WGSCRS)
plot(Ice_Free_WG)

#2.3 crop using AP_ext in world crs
AP_Ice_WG <- crop(Ice_Free_WG, AP_Ext)
plot(Ap_Ice_WG)

#2.4 use project function to convert back to ATACRS
AP_Ice_ATA <- project(AP_Ice_WG, ATACRS)
plot(AP_Ice_ATA)

### LOAD ATA_BIOCLIM_PROJ
ATA_Bioclim_Proj <- rast("ATA_Bioclim_Proj.tif")

#2.5 use mask function to hide ice covered area

                  
Bioclim_masked <- ATA_Bioclim_Proj %>% 
  mask(AP_Ice_ATA)

#2.6 plots to compare the difference
plot(ATA_Bioclim_Proj[[1]], main = "Projected Bio 1")
plot(Bioclim_masked[[1]], main = "Projected Bio 1")

plot(Bioclim_masked)


#2.7 save projection layer folders for maxent input

my_path <- "C:/Users/Student/OneDrive - University of Edinburgh/Documents/EDINBURGH/Dissertation/Modelling/Code/Lichen_Diss/MaxEnt_layers/"

# create present layers folder
for (i in 1:nlyr(Bioclim_masked)) {
  
  # Get the layer name from the raster stack
  layer_name <- names(Bioclim_masked)[i]
  
  # Combine path + name + extension
  # This creates: "C:/.../Projection_Layers/bio1.asc"
  full_destination <- paste0(my_path, layer_name, ".asc")
  
  # Save the file
  writeRaster(Bioclim_masked[[i]], 
              filename = full_destination, 
              filetype = "AAIGrid", 
              NAflag = -9999, 
              datatype = "FLT4S", 
              overwrite = TRUE)
  
  print(paste("Saved to:", full_destination))
}




#----------------------------------------------------------------------------------------------------------------------------------------

#convert Lichen data in vector form to get it to same format as env values. (same extent and projection)

### 3: Lichens
#3:1: load antarctic dataset and isolate lichens
Lichen <- AAS_4296_Biodiversity_IceFree_Antarctica_DB %>% 
  filter(phylum %in% c("Ascomycota")) %>%
  rename(species = scientificName)



#3:2: turn into vector (convert spatially)
Vect_Lichen <- vect(as.data.frame(Lichen), geom=c("decimalLongitude", "decimalLatitude"), WGSCRS)

#3:3: crop to the extent of Antarctic peninsula
ATA_Lichen <- crop(Vect_Lichen, AP_Ext)

plot(ATA_Lichen)

#3:4: project to correct ATA CRS
ATA_Lichen_Proj <- project(ATA_Lichen, ATACRS)

#3:5: visualise
plot(ATA_Bioclim_Proj[[1]], main = "Projected Bio 1")
points(ATA_Lichen_Proj)


#--------------------------------------------------------------------------------------------------

### 4: Combine the 2 data sets (environment parameters and lichens)

#4.1 step convert both back to data frames
Lichen_DF <- as.data.frame(ATA_Lichen_Proj, geom = "XY") %>% 
  select(species, x, y)

# this gives environmental values for where the lichens where found and only where the lichens were found. Also using the ice_masked bioclim just in case any of the coords where in the ice. This will create an Na value in the row, this will also occur if the coordinate is in the sea.
Env_Values <- extract(Bioclim_masked, ATA_Lichen_Proj, ID = FALSE)


#4.2 combine these two data sets
Full_Dataset <- cbind(Lichen_DF, Env_Values)
head(Full_Dataset)
summary(Full_Dataset$wc2.1_2.5m_bio_1)

#------------------------------------------------------------------------------------

### 5: create two data sets one for species of interest and one for the background points.
#THIS IS THE PART THAT CHANGES FOR EACH SPECIES OF INTEREST

### separate out species of interest. end with two data sets in exact same format. large dataset, hide species names this acts as targeted background.

#5.1 separate out target species

target_species <- "Acarospora gwynnii"
#Acarospora gwynnii
#Usnea antarctica
#Umbilicaria decussata

#5.2 target data set

Target_A_g <- Full_Dataset %>% 
  filter(species == target_species) %>% 
  na.omit()

head(Target_A_g)
summary(Target_A_g$wc2.1_2.5m_bio_1)

#test to visaulise where points would be

Target_proj <- Lichen %>% 
  filter(species %in% c("Acarospora gwynnii"))

Target_proj <- vect(as.data.frame(Target_proj), geom=c("decimalLongitude", "decimalLatitude"), WGSCRS) %>% 
  crop(AP_Ext) 

Target_proj <- project(Target_proj, ATACRS)


plot(ATA_Bioclim_Proj[[1]], main = "Projected Bio 1")
points(Target_proj)

#removing na as it will not be read by maxent and will break it. It could also be converted to -9999. By using bioclim_asc, which has already had na converted to -9999 but also need to make change any Na's that have occured from being in the water to -9999 as well.

#5.3 background dataset
Background_A_g <- Full_Dataset %>% 
  filter(species != target_species) %>% 
  mutate(species = "background") %>% 
  na.omit()


#checking that we don't have any Na's or -9999
summary(Background_A_g)

#5.4 save the two new dataset's as CSV files ready for inputting into maxent. Species of interest data on left, background data on right. 

write.csv(Target_A_g, "Presence_AG.csv", row.names = FALSE)
write.csv(Background_A_g, "Absence_AG.csv", row.names = FALSE)

#------------------------------------------------------------------------------------------------------------------------------------


## GO TO MAXENT, Target_Data saved CSV file input on left, BAckground data CSV file inout on right. Ice-free masked environmental layer projetion input into project layer directory. SORTED!!!!
## Repeat for other species.

## Have 4,500 background points. Research needed to be done on which species to pick and whether there is restraints on presence data required.

### Secondly how to test model. can model outputs be tested in r after doing them manually or is this a reason to make it automatic and use the dismo package.
# jackknifing potentially - need to ask hannah
# AIC values etc to understand model fit.

#check asc

#checking maxent results.
#should probably implement dismo package so that it is done automatically rather than having to manually open maxent and input.

library(terra)
result <- rast("C:/Users/Student/Documents/edinburgh/R_studio/Maxent_outputs/Acarospora_gwynnii_MaxEnt_Layers.asc")

setMinMax(result)

# 3. Check the values now
print(minmax(result))

# 4. Check if the map is just all one value (e.g., all 0 or all 1)
# This calculates the global range properly
global(result, fun="range", na.rm=TRUE)


#--------------------------------------------------------------------------------------------------------------------------------------
#### FUTURE

WGSCRS <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ATACRS <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"


# 6. Load environmental layers for the future. create ASC file. mask with future ice cover data

####6.1 import future environmental layers

  
Future_Bioclim <- cmip6_world(
  model = "EC-Earth3-Veg", 
  ssp = "585", 
  time = "2081-2100", 
  var = "bioc", 
  res = 2.5, 
  path = tempdir()
)


print(Future_Bioclim)

#define extent again
AP_Ext <- ext(x = c(-80, -55, -78, -60))

#crop
  F_AP_Bioclim <- crop(Future_Bioclim, AP_Ext)

#project
F_ATA_Bioclim_Proj <- project(F_AP_Bioclim, ATACRS)
plot(F_ATA_Bioclim_Proj)

###6.2 import future ice free data. 8.5 scenario - best fit

F_Ice_free_area <- vect("C:/Users/Student/OneDrive - University of Edinburgh/Documents/EDINBURGH/Dissertation/Data/Penisula_IceCover/AAS_4297_Future_Ice-free_Layers/AAS_4297_Ice_Free_Shapefiles_8.5_Best") 
plot(F_Ice_free_area)

F_Ice_Free_WG <- project(F_Ice_free_area, WGSCRS)
plot(F_Ice_Free_WG)

#2.3 crop using AP_ext in wgscrs
F_AP_Ice_WG <- crop(F_Ice_Free_WG, AP_Ext)
plot(F_AP_Ice_WG)

#2.4 use project function to convert back to ATACRS
F_AP_Ice_ATA <- project(F_AP_Ice_WG, ATACRS)
plot(F_AP_Ice_ATA)


#6.3 mask
F_asc <- F_ATA_Bioclim_Proj %>% 
  mask(F_AP_Ice_ATA)



names(F_asc) <- c("wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_2", "wc2.1_2.5m_bio_3", "wc2.1_2.5m_bio_4", "wc2.1_2.5m_bio_5", "wc2.1_2.5m_bio_6", "wc2.1_2.5m_bio_7", "wc2.1_2.5m_bio_8", "wc2.1_2.5m_bio_9", "wc2.1_2.5m_bio_10", "wc2.1_2.5m_bio_11", "wc2.1_2.5m_bio_12", "wc2.1_2.5m_bio_13", "wc2.1_2.5m_bio_14", "wc2.1_2.5m_bio_15", "wc2.1_2.5m_bio_16", "wc2.1_2.5m_bio_17", "wc2.1_2.5m_bio_18", "wc2.1_2.5m_bio_19")


#6.4 save as folder

F_path <- "C:/Users/Student/OneDrive - University of Edinburgh/Documents/EDINBURGH/Dissertation/Code/Lichen_Diss/F_MaxEnt_layers/"


for (i in 1:nlyr(F_asc)) {
  
  # Get the layer name from the raster stack
  layer_name <- names(F_asc)[i]
  
  # Combine path + name + extension
  # This creates: "C:/.../Projection_Layers/bio1.asc"
  full_destination <- paste0(F_path, layer_name, ".asc")
  
  # Save the file
  writeRaster(F_asc[[i]], 
              filename = full_destination, 
              filetype = "AAIGrid", 
              NAflag = -9999, 
              datatype = "FLT4S", 
              overwrite = TRUE)
  
  print(paste("Saved to:", full_destination))
}

#### very little difference in ice cover between future and present. the colour changes slighlty in certain areas. I think might be to do with resolution need to check what resolution is in J. LEE 2017 paper. 
# Might also be good idea to rename everything to bioclim1 so its standardised to something simple rather than the wc_2.5 (which represents resolution as resolution might be changing)
## also change worldclim data and lichen data to ATACRS and crop using extent in ATACRS rather than WGSCRS. This is rather than swapping ice cover data from ATACRS to WGSCRS to then crop to then convert back to ATACRS
