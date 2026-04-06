### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Note: this script really only needs to be run when updating the master
# shapefile that has the invaded counties/county-equivalents in the US
# and Canada
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 



# Load in county shapefile -----------------------------------------------------
#
#
#---
# Obtained from USDA FS on June 23, 2025 (2024 data)
# Obtained from USDA FS on January 20, 2026 (2025 data)
BLD_cty_2024 <- read_excel("spreadsheets/county_level_BLD_data/US BLD Counties by Year 2024.xlsx")
nrow(BLD_cty_2024)
BLD_cty <- read_excel("spreadsheets/county_level_BLD_data/US BLD Counties by Year 2025.xlsx")
nrow(BLD_cty)
BLD_cty$State_province <- toupper(BLD_cty$State_province)
BLD_cty <- as_tibble(BLD_cty)
head(BLD_cty)

# correcting misspelled or "mislabeled" counties
# misspelling
BLD_cty[BLD_cty$County=="Hurong","County"] <- "Huron"
# Quinte West is located in Hastings County
BLD_cty[BLD_cty$County=="Quinte West","County"] <- "Hastings"
# Mississauga is located in Peel County
BLD_cty[BLD_cty$County=="Mississauga","County"] <- "Peel"
# removing underscore
BLD_cty[BLD_cty$County=="Anne_Arundel","County"] <- "Anne Arundel"
# misspelling
BLD_cty[BLD_cty$County=="Muskingham","County"] <- "Muskingum"
# DE assumed to Delaware County in PA
BLD_cty[BLD_cty$County=="DE","County"] <- "Delaware"
# misspelling
BLD_cty[BLD_cty$County=="Montogomery","County"] <- "Montgomery"
# missing an apostrophe
BLD_cty[BLD_cty$County=="Prince Georges","County"] <- "Prince George's"
# Washington DC
BLD_cty[BLD_cty$County=="Washington DC","County"] <- "District Of Columbia"
BLD_cty[BLD_cty$State_province=="WASHINGTON DC","State_province"] <- "DC"

# this county was missing year, so we input it (confirmed year by looking at USDA FS maps online)
BLD_cty[which(BLD_cty$County=="Kent" & BLD_cty$State_province=="DE"),"BLD_Year"] <- 2025

# this county is missing - determined from later visual inspection of maps
df_BedfordPA <- data.frame(State_province="PA", County="Bedford",  BLD_Year=2021)
BLD_cty <- rbind.data.frame(BLD_cty,df_BedfordPA)

# remove duplicates
# Counties should not be repeated within a state, so create an ID variable
# that identifies unique pairings of states X counties
BLD_cty$ID <- paste(BLD_cty$State_province,BLD_cty$County, sep="_")
BLD_cty$ID_YR <- paste(BLD_cty$State_province,BLD_cty$County,BLD_cty$BLD_Year, sep="_")
nrow(BLD_cty)

# the are repeated observations of "ID" - a variable indicating pairings
# of State_province and County. This code sorts the observations by invasion year
# and then takes the first occurrence of a unique county Id
BLD_cty <- BLD_cty[order(BLD_cty$BLD_Year),]
BLD_cty_first <- BLD_cty[!duplicated(BLD_cty$ID), ]


# these are harder to find duplicates/problems, let's 'em remove now:
# De Haven is an unincorporated community in northern Frederick County,
# and we have Frederick County in the data. So, delete "De Haven county"
# We have Monmouth New Jersey - invaded in 2021 - whereas monmouth New Jersey
# was invaded in 2022 - delete the 2022, lower case observation
# there is a "Test" county in OH - remove
BLD_cty_first <- BLD_cty_first[BLD_cty_first$County %!in% c("De Haven county", "monmouth", "Test"), ]
nrow(BLD_cty_first)
length(unique(BLD_cty_first$ID)) == nrow(BLD_cty_first) # confirmation duplicates removed
length(unique(BLD_cty_first$ID_YR))== nrow(BLD_cty_first) # confirmation duplicates removed

# check number of US counties for QAQC
BLD_US_only <- BLD_cty_first[which(BLD_cty_first$State_province  != "ONTARIO"),]
nrow(BLD_US_only)


# load in FIPS identifiers
data("fips_codes", package = "tidycensus")
head(fips_codes)
fips_clean <- fips_codes %>%
  mutate(
    State_province = state,
    County = county %>%
      str_remove_all(regex(" county| parish| borough| census area| municipality| city and borough", ignore_case = TRUE)) %>%
      str_trim() %>%
      str_to_title()
  )


# create a copy of database
BLD_all <- BLD_cty_first

# combine BLD data with FIPS database 
BLD_all_joined <- BLD_all %>%
  left_join(fips_clean, by = c("State_province", "County"))

# combine columns for state and county codes (5 digit FIPS)
BLD_all_joined$full_fips <- ifelse(is.na(BLD_all_joined$state_code) | is.na(BLD_all_joined$county_code),
                                   NA, paste0(BLD_all_joined$state_code, BLD_all_joined$county_code))

#
BLD_all_joined[is.na(BLD_all_joined$full_fips),]
unique(BLD_all_joined[is.na(BLD_all_joined$full_fips),"State_province"])

# FIPS fully assigned to all US counties?
sum(table(BLD_all_joined$full_fips)) == nrow(BLD_US_only)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 




# Load in US shapefile ---------------------------------------------------------
#
#
#---
# https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
US_countiesall <- st_read(dsn="gis_data/US_shapefiles/county", layer = "cb_2018_us_county_500k")
# plot(st_geometry(US_countiesall))
US_counties <- US_countiesall %>% filter(!STATEFP  %in% c("15", "02", "60", "66", "69", "72", "78"))
# plot(st_geometry(US_counties))

# get FIPS set up
US_counties$full_fips <- paste0(US_counties$STATEFP, US_counties$COUNTYFP) 

# make sure all FIPS have 5 digits
BLD_all_joined$full_fips <- str_pad(BLD_all_joined$full_fips, width = 5, side = "left", pad = "0")
table(nchar(BLD_all_joined$full_fips))

# combine invasion data with shapefile
combined_sf_matched_US <- US_counties %>%
  left_join(BLD_all_joined, by = "full_fips")
# plot(st_geometry(combined_sf_matched_US))
table(combined_sf_matched_US$BLD_Year)
sum(table(combined_sf_matched_US$BLD_Year)) # US counties invaded
sum(table(BLD_all_joined[which(BLD_all_joined$full_fips %in% combined_sf_matched_US$full_fips),"State_province"]))

# Load in CAN shapefile --------------------------------------------------------
#
#
#---
# https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/index2021-eng.cfm?year=21
CAN_countiesall1 <- st_read("gis_data/CAN_shapefiles/census_divisions/lcd_000b21a_e.shp")
# plot(st_geometry(CAN_countiesall1))

# change CDNAME to County and remove leading/trailing whitespace
CAN_counties_all <- CAN_countiesall1 %>%
  rename(County = CDNAME) %>%
  mutate(
    County = str_trim(str_to_title(County))
  )

# Load in Canada "county" codes
# https://www23.statcan.gc.ca/imdb/p3VD.pl?Function=getVD&TVD=134850&CVD=134851&CPV=35&CST=01012001&CLV=1&MLV=3
# Manually transribed "Census division" column into "County" column in Excel
CAN_cty <- read_excel("spreadsheets/CSD_codes/CSD_codes.xlsx")
CAN_cty <- CAN_cty[,c("CDUID", "County")]

# get invaded counties in Canada
table(BLD_all_joined$State_province) # just Ontario has been invaded
BLD_CAN_only <- BLD_all_joined[which(BLD_all_joined$State_province  == "ONTARIO"),]
nrow(BLD_CAN_only) + sum(table(combined_sf_matched_US$BLD_Year)) == nrow(BLD_all_joined)

# combine invasion data with county codes
BLD_CAN_only <- BLD_CAN_only %>%
  left_join(CAN_cty, by = "County")

# combine Canada shapefile with codes/invasion data
combined_sf_matched_CAN <- CAN_counties_all %>%
  left_join(BLD_CAN_only, by = "CDUID")

# QAQC: check all invaded Canadian counties transferred 
sum(table(combined_sf_matched_CAN$BLD_Year)) == sum(table(BLD_CAN_only$CDUID))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


       

# Combine and save shapefiles --------------------------------------------------
#
#
#---
# convert Canada shapefile to same projection
combined_sf_matched_CAN <- st_transform(combined_sf_matched_CAN, st_crs(combined_sf_matched_US))

# select pertinent columns in each shapefile, harmonize names
US_fin_sf <- combined_sf_matched_US  %>% dplyr::select(full_fips, BLD_Year, State_province, County, geometry)
CAN_fin_sf <- combined_sf_matched_CAN %>% dplyr::select(full_fips=CDUID, BLD_Year, State_province, County=County.y, geometry)

# combine shapefules, project to Albers Equal Area
USCAN <- bind_rows(US_fin_sf, CAN_fin_sf)
USCAN_albers <- st_transform(USCAN, crs=5070)

# QAQC: did all invaded counties make it?
sum(table(USCAN_albers$BLD_Year)) == nrow(BLD_US_only)+nrow(BLD_CAN_only)

# Shorten column names and export shapefile
USCAN_albers_export <- USCAN_albers %>% dplyr::select(FIPS=full_fips, 
                                                 BLD_YR=BLD_Year, 
                                                 ST_PR=State_province, 
                                                 CTY=County, 
                                                 geometry)


#plot(st_geometry(USCAN_albers))
st_write(USCAN_albers_export, "gis_data/invasion_data_USCAN/USCAN_albers.shp", delete_layer = TRUE)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 






# Aggregate raster of beech distribution --------------------------------------------------
#
#
#---
#
# this file must be downloaded or requested from Sam - too large for storage on GitHub
# https://data.fs.usda.gov/geodata/rastergateway/bigmap/index.php
beech_rast <- terra::rast('raw_beech_raster/Hosted_AGB_0531_2018_AMERICAN_BEECH_06052023045637.tif')

mean_na_rm <- function(data_vector) {
  # Calculate the mean while removing NA values
  result <- mean(data_vector, na.rm = TRUE)
  # Return the result
  return(result)
}

beech_rast_10km <- aggregate(beech_rast, fact = 33, fun=mean_na_rm)
plot(beech_rast)
plot(beech_rast_10km)
names(beech_rast_10km) <- "agb_tons_per_acre"
# convert from tons per acre to tons per hectare
beech_rast_10km$agb_tons_per_ha <- beech_rast_10km$agb_tons_per_acre*2.47105
beech_rast_10km_ha <- beech_rast_10km[["agb_tons_per_ha"]]

writeRaster(beech_rast_10km_ha, filename = 'gis_data/beech_distribution/beech_rast_10km.tif', overwrite=T)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 