### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# This script is used for calculating spread rates of beech leaf disease using
# three different approaches: effective range radius, distance regression, and
# boundary displacement
#
# IMPORTANT!!! Downloading the GIS data from Sam (or Zenodo) (directory gis_data) and storing
# it in the same directory as the other folders is necessary for this code to run
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 





# Load in shapefile of USCAN counties with invasion information ----------------
#
#
#---
source("r_scripts/1_packages.R")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 



# Load in shapefile of USCAN counties with invasion information ----------------
#
#
#---
USCAN_albers <- st_read("gis_data/invasion_data_USCAN/USCAN_ALBERS.shp")

# first detection location for calculating distance metrics
location_of_LakeOH <- data.frame(long = -81.2400, lat = 41.8200, NAME = "First detection")
# project into LatLon
pts_LakeOH <- st_as_sf(location_of_LakeOH, coords = c("long", "lat"), crs = 4326)
# project into Albers
pts_LakeOH_ALBERS <- st_transform(pts_LakeOH, 5070)

max_year_of_invasion <- max(USCAN_albers$BLD_YR, na.rm=T)
#plot(st_geometry(USCAN_albers), col = "lightgrey", border = "black", main = "US + Canada Map")
#plot(st_geometry(pts_LakeOH_ALBERS), col = "red", pch = 19, cex = 2, add = TRUE)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 






# Create a buffered area for mapping -------------------------------------------
#
#
#---
# ONLY EDIT/RUN IF INTENDING TO EDIT THE RESULTING FILES
# PLEASE DISCUSS BEFORE DOING SO
# computationally demanding 
if(FALSE){
  
  table(is.na(USCAN_albers$BLD_YR))
  
  BLD_invaded_areas <- USCAN_albers[which(!is.na(USCAN_albers$BLD_YR)),]
  # create a buffered polygon (invaded area + 500 km)
  counties_buff_500km <- st_buffer(BLD_invaded_areas, 500*1000)%>%  # 500 km
    st_union() %>% # unite to a geometry object
    st_sf() # make the geometry a data frame object
  # select counties within buffered area
  counties_within_buffer_pts <- st_intersection(st_centroid(USCAN_albers), counties_buff_500km)
  counties_within_buffer <- USCAN_albers[which(USCAN_albers$FIPS %in% counties_within_buffer_pts$FIPS),]
  nrow(counties_within_buffer)
  
  # crop shapefile by invaded counties (later used for plotting)
  BLD_ALBERS_cropped <- st_crop(USCAN_albers, BLD_invaded_areas)
  
  # save shapefiles
  #plot(st_geometry(counties_within_buffer))
  st_write(counties_within_buffer, "gis_data/buffered_invasion/buffered_invasion.shp", delete_layer = TRUE) 
  st_write(BLD_ALBERS_cropped, "gis_data/buffered_invasion/BLD_ALBERS_cropped.shp", delete_layer = TRUE) 
  beep(1)
}
# load in the buffered and cropped shapefiles (mostly used for plotting)
counties_within_buffer <- st_read("gis_data/buffered_invasion/buffered_invasion.shp")
BLD_ALBERS_cropped <- st_read("gis_data/buffered_invasion/BLD_ALBERS_cropped.shp")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


# Distances between counties, adjoining vs. isolated spread (setup) ------------
#
#
#---
# ONLY EDIT/RUN IF INTENDING TO EDIT THE RESULTING FILES
# PLEASE DISCUSS BEFORE DOING SO
if(FALSE){
  invaded_counties_all <- USCAN_albers %>%
    dplyr::filter(!is.na(BLD_YR))
  nrow(invaded_counties_all)
  
  # create spatial neighborhood for getting neighbors of each county
  # second order (queens case neighborhood)
  sec_order <- poly2nb(invaded_counties_all, queen = T, row.names=invaded_counties_all$FIPS)
  sec_order[[1]]
  
  ## For loop --------------------------------------------------
  # for loop determining distances between counties and whether
  # a county was isolated from or adjacent to invaded area 
  # should start in second year of invasion records, since 
  # counties invaded in the first year would
  # be definition be isolated
  # i <- 2013
  for(i in 2013:max_year_of_invasion){
    
    # current invaded counties in time i
    invaded_in_i <- invaded_counties_all[which(invaded_counties_all$BLD_YR %in% i),]
    
    # all counties invaded across previous years
    prev_invaded_counties <- invaded_counties_all[which(invaded_counties_all$BLD_YR %in% (2012:(i-1))),]
    
    
    if(nrow(invaded_in_i) > 0){ # if there counties invaded in year i, proceed into loop
      for(j in 1:nrow(invaded_in_i)){ # for each invaded county, loop through
        #j <- 1
        
        # get current county and all previously invaded points
        new_invaded_point <- invaded_in_i[j,"FIPS"]
        if(i == 2013){ # if second year of invasion, assume discovery location is nearest previously invaded location
          prev_invaded_points <- pts_LakeOH_ALBERS
        }else{ # otherwise, get all previous invaded counties
          prev_invaded_points <-  prev_invaded_counties
        }
        
        # get distance between current county and discovery location
        pts_discovery_location <- suppressWarnings(pointDistance(st_centroid(new_invaded_point),st_centroid(pts_LakeOH_ALBERS),lonlat=F))
        # get minimum distance to discovery location (DL = discovery location)
        invaded_in_i$DtoDL_km[j] <- min(pts_discovery_location)[1]/1000 # convert to km
        
        # find minimum distance to previous year's invasion boundary
        dist_vec_centroids <- suppressWarnings(pointDistance(st_centroid(new_invaded_point),st_centroid(prev_invaded_points),lonlat=F))
        invaded_in_i$D_BDY_km[j] <- min(dist_vec_centroids)[1]/1000
        
        # get closest invaded county to county j
        val <- which(dist_vec_centroids == min(dist_vec_centroids)[1], arr.ind = TRUE) # if tied, take first observation in the tie
        closest_county_BNDRY <-  prev_invaded_counties[val,]
        invaded_in_i$NR_CNTY[j] <- as.character(closest_county_BNDRY$FIPS)
        
        # determine whether a neighbor was invaded
        loc_in_vec <- which(invaded_counties_all$FIPS %in% invaded_in_i[j, "FIPS"]) # where is the current county located in data frame
        neighbs <- invaded_counties_all[sec_order[[loc_in_vec]],] # get current county's neighbors
        # determine if current county had neighbors in the previous year that were infested
        neighbs_prev <- prev_invaded_counties[which(prev_invaded_counties$FIPS %in% neighbs$FIPS),]
        # make neighbor assignment
        if(nrow(neighbs_prev) == 0){invaded_in_i$STAT[j] <- "iso"} else{ # if county j had previously invaded neighbors, assign iso
          invaded_in_i$STAT[j] <- "adj" # otherwise assign adj
        }}
      
      # specifying variables to extract and export from invaded_in_i dataframe
      vec_variables <- c("FIPS", "D_BDY_km", "DtoDL_km", "STAT", "BLD_YR", "NR_CNTY")
      
      if(i == 2013){ # if in first iteration of the loop, extract the columns the dataframe
        invaded_counties <- as.data.frame(invaded_in_i[, paste(vec_variables)])
      } else { # otherwise, append the dataframe
        invaded_counties <- rbind.data.frame(invaded_counties,as.data.frame(invaded_in_i[, paste(vec_variables)]))
      }
    }
    
    cat("Year", i,"out of", max_year_of_invasion,"\n") 
  }
  
  
  # Add data into shapefile --------------------------------------------------
  # create new shapefile
  BLD_sprd_ALBERS <- USCAN_albers
  # add empty columns, populate with NA
  BLD_sprd_ALBERS$STAT <- NA
  BLD_sprd_ALBERS$D_BDY_km <- NA
  BLD_sprd_ALBERS$DtoDL_km <- NA
  BLD_sprd_ALBERS$NR_CNTY <- NA
  
  # for loop extract data from invaded counties data frame
  for(i in 1:nrow(BLD_sprd_ALBERS)){
    curr_fips <- BLD_sprd_ALBERS$FIPS[i]
    if( length(which(invaded_counties$FIPS %in% curr_fips)) > 0 ){
      BLD_sprd_ALBERS$STAT[i] <- invaded_counties[which(invaded_counties$FIPS %in% curr_fips), "STAT"]
      BLD_sprd_ALBERS$D_BDY_km[i] <- invaded_counties[which(invaded_counties$FIPS %in% curr_fips), "D_BDY_km"]
      BLD_sprd_ALBERS$DtoDL_km[i] <- invaded_counties[which(invaded_counties$FIPS %in% curr_fips), "DtoDL_km"]
      BLD_sprd_ALBERS$NR_CNTY[i] <- invaded_counties[which(invaded_counties$FIPS %in% curr_fips), "NR_CNTY"]
    }
  }
  st_write(BLD_sprd_ALBERS, "gis_data/adjacent_isolated/BLD_sprd_ALBERS.shp", delete_layer=T)
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 









# Load in USCAN shapefile: states and provinces --------------------------------
#
#
#---
# mainly used for making pretty maps
# this chunk loads in, combines, and save the files
# ONLY EDIT/RUN IF INTENDING TO EDIT THE RESULTING FILES
# PLEASE DISCUSS BEFORE DOING SO
if(FALSE){
  # https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
  states = st_read("gis_data/US_shapefiles/states/cb_2018_us_state_20m.shp")
  states_albers <- st_transform(states, st_crs(BLD_ALBERS_cropped))
  states_albers <- states_albers %>% dplyr::filter(NAME %!in% c("Puerto Rico", "Hawaii")) %>% 
    dplyr::select(FIPPR = STATEFP,
                  GEOID = AFFGEOID,
                  NAME = NAME)
  #
  # https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/index2021-eng.cfm?year=21
  provinces = st_read("gis_data/CAN_shapefiles/provinces/lpr_000b21a_e.shp")
  provinces_albers <- st_transform(provinces, st_crs(BLD_ALBERS_cropped))
  provinces_albers <- provinces_albers %>%  dplyr::select(FIPPR = PRUID,
                                                          GEOID = DGUID,
                                                          NAME = PRENAME)
  
  
  states_provinces_albers <- bind_rows(states_albers, provinces_albers)
  #plot(st_geometry(states_provinces_albers))
  # shapefile cropped to invaded areas
  states_provinces_albers_crop <- st_crop(states_provinces_albers, BLD_invaded_areas)
  #plot(st_geometry(states_provinces_albers_crop))
  st_write(states_provinces_albers, "gis_data/USCAN_combined/states_provinces_albers.shp", delete_layer=T)
  st_write(states_provinces_albers_crop, "gis_data/USCAN_combined/states_provinces_albers_crop.shp", delete_layer=T)
}
states_provinces_albers = st_read("gis_data/USCAN_combined/states_provinces_albers.shp")
states_provinces_albers_crop = st_read("gis_data/USCAN_combined/states_provinces_albers_crop.shp")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 





# Figure: invaded counties -----------------------------------------------------
#
#
#---
my_BLD_colors <- scico(14, palette = 'managua')
plot(1:14, pch=15, cex=6, col=my_BLD_colors)
# https://www.r-bloggers.com/2019/12/inset-maps-with-ggplot2/
bld_bb = st_as_sfc(st_bbox(states_provinces_albers_crop))

### ### ### ### ### ### ### ###
# Inset map
### ### ### ### ### ### ### ###
inset <- ggplot() +
  geom_sf(data = states_provinces_albers, fill="white") + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data = bld_bb, fill = alpha("red",0.5), color = "red", size = 3) +
  coord_sf(datum=st_crs(BLD_ALBERS_cropped))

### ### ### ### ### ### ### ###
# Invaded area time series map
### ### ### ### ### ### ### ###
BLD_ALBERS_cropped$f_BLD_YR <- factor(BLD_ALBERS_cropped$BLD_YR)
length(table(BLD_ALBERS_cropped$f_BLD_YR))
invaded_area_continuous <- ggplot() +
  geom_sf(data = BLD_ALBERS_cropped, aes(fill = f_BLD_YR), color="light gray") + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data=states_provinces_albers_crop, fill="transparent", color="black", lwd=0.1)+
  scale_fill_manual(values = my_BLD_colors, 
                    name="", na.value = "white", na.translate = F)+ 
  #breaks = seq(2012, 2024, by = 3), 
  #  guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  theme(legend.position =c(0.9,0.4),legend.key.size = unit(0.2, 'cm'), # c(0.7,0.3)
        legend.text = element_text(size=8), legend.justification = "left")+
  geom_segment(
    data = pts_LakeOH_ALBERS,
    aes(x = st_coordinates(pts_LakeOH_ALBERS$geometry)[1]-90000, y = st_coordinates(pts_LakeOH_ALBERS$geometry)[2]+10000, 
        xend = st_coordinates(pts_LakeOH_ALBERS$geometry)[1], yend = st_coordinates(pts_LakeOH_ALBERS$geometry)[2],
    ),
    arrow = arrow(length = unit(0.1, "cm"), type = "closed"), # Adds the arrowhead
    linewidth = 1,
    color="red") +
  coord_sf(datum=st_crs(BLD_ALBERS_cropped))
#
# invaded_area_continuous
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 







# Figure: host distribution ----------------------------------------------------
#
#
#---
if(F){ # change to T if you want to plot beech distribution
  # Load biomass data for American beech (Fagus grandifolia)
  # original data from  (see 2_curation) https://data.fs.usda.gov/geodata/rastergateway/bigmap/index.php
  beech_agb <- rast('gis_data/beech_distribution/beech_rast_10km.tif') 
  beech_agb[beech_agb < 1] <- NA
  
  
  plot(beech_agb)
  states_provinces_albers_r <- st_transform(states_provinces_albers, st_crs(beech_agb))
  plot(st_geometry(states_provinces_albers_r), add=T, col="transparent", border="black")
  
  # Little's Range Maps Source
  # https://github.com/wpetry/USTreeAtlas/blob/main/shp/fagugran/fagugran.shx
  little_beech = st_read("gis_data/beech_distribution/fagugran.shp")
  little_beech <- st_set_crs(little_beech, "epsg:4326") 
  little_beech_albers <- st_transform(little_beech, st_crs(beech_agb))
  little_beech_albers <- little_beech_albers[which(little_beech_albers$CODE != 0),] # remove water 
  # remove points in Mexico
  little_beech_albers_crop <- st_crop(little_beech_albers, xmin= -380620 + 3555334*0.07, ymin= -381632.5, xmax= 2748212, ymax= 3173701)
  plot(st_geometry(little_beech_albers_crop), col="light green")
  plot(st_geometry(states_provinces_albers_r), add=T, col="transparent", border="black")
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 






# Adjacent versus isolated spread: analysis and figure -------------------------
#
#
#---
contigUSCAN_ALBERS_BLD <- st_read("gis_data/adjacent_isolated/BLD_sprd_ALBERS.shp")

summary(contigUSCAN_ALBERS_BLD)

BLD.adj.iso <- contigUSCAN_ALBERS_BLD[which(contigUSCAN_ALBERS_BLD$STAT %in% c("adj","iso")),]

contigUSCAN_ALBERS_BLD$JumpInv <- ifelse(contigUSCAN_ALBERS_BLD$STAT == "iso", "Non-contig.", contigUSCAN_ALBERS_BLD$STAT)
contigUSCAN_ALBERS_BLD$JumpInv <- ifelse(contigUSCAN_ALBERS_BLD$JumpInv == "adj", "Contig.", contigUSCAN_ALBERS_BLD$JumpInv)

contigUSCAN_ALBERS_BLD$color <- NA
contigUSCAN_ALBERS_BLD$color[contigUSCAN_ALBERS_BLD$JumpInv == "Contig."] <- "gray"
contigUSCAN_ALBERS_BLD$color[contigUSCAN_ALBERS_BLD$JumpInv == "Non-contig."] <- "tomato2"
contigUSCAN_ALBERS_BLD$color[is.na(contigUSCAN_ALBERS_BLD$JumpInv)] <- "white"


isolated_counties <- BLD.adj.iso[which(BLD.adj.iso$STAT == "iso"),]
summary(isolated_counties)
isolated_counties[which(isolated_counties$D_BDY_km %in% min(isolated_counties$D_BDY_km)),]
nrow(BLD.adj.iso[which(BLD.adj.iso$STAT == "iso"),])
nrow(BLD.adj.iso[which(BLD.adj.iso$STAT == "adj"),])
nrow(BLD.adj.iso[which(BLD.adj.iso$STAT %in% c("iso","adj")),])
isolated_counties[order(-isolated_counties$D_BDY_km),]
isolated_counties[order(isolated_counties$BLD_YR),]
isolated_counties[which(isolated_counties$ST_PR =="VA"),]
isolated_counties[which(isolated_counties$ST_PR =="NC"),]


### ### ### ### ### ### ### ###
# Adjacent vs isolated map
### ### ### ### ### ### ### ###
isoadj_cty <- BLD.adj.iso %>%  as.data.frame() %>% group_by(BLD_YR,STAT) %>% summarize(counties_invaded = n())
add_year0 <- data.frame(BLD_YR=2013, STAT="iso",  counties_invaded=0)
isoadj_cty <- rbind.data.frame(add_year0,isoadj_cty)

sum(table(BLD.adj.iso$STAT))
table(BLD.adj.iso$STAT)
year_county_summary <- as.data.frame(BLD.adj.iso) %>% group_by(BLD_YR) %>% summarize(counties_invaded = n())
summary(year_county_summary$counties_invaded)
stderr(year_county_summary$counties_invaded)

adj_col <- "#73C0E7" #my_BLD_colors[2]
iso_col <- "#E6A85A"
transp_gray <- alpha("light gray",0.3) 
BLD.adj.iso_fin.crp <- st_crop(contigUSCAN_ALBERS_BLD, BLD_ALBERS_cropped)
inv.isoadj_map <- ggplot() +
  geom_sf(data = BLD.adj.iso_fin.crp, aes(fill = JumpInv), color=transp_gray,lwd=0.0001) + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data=states_provinces_albers_crop, fill="transparent", color="black", lwd=0.1)+
  scale_fill_manual(values= c(adj_col,iso_col),
                    name="", na.translate=FALSE)+
  # theme(legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=8), legend.justification = "left",
  #       legend.position=c(0.3,0.25))
  theme(legend.position = c(0.8,0.40),legend.key.size = unit(0.2, 'cm'), # 
        legend.text = element_text(size=8), legend.justification = "left")

gg_inset_map1_isoadj = ggdraw() +
  draw_plot(inv.isoadj_map) +
  draw_plot(inset, x = 0.62, y = 0.05, width = 0.25, height = 0.25)


### ### ### ### ### ### ### ###
# Adjacent vs isolated area graph (OLD)
### ### ### ### ### ### ### ###

# 
# inv.isoadj_graph_AREA <- 
#   ggplot(isoadj_cty,  aes(x=BLD_YR, y=counties_invaded, fill=STAT)) +
#   geom_area( aes(fill=STAT), position="stack")+
#   labs(title="",x="Year", y = "Counties invaded")+theme_clean+
#   scale_fill_manual(values=c(adj_col, iso_col),labels = c("Contiguous", "Non-contiguous"))+
#   scale_x_continuous(expand = c(0, 0))+ # x axis is getting cut-off
#   scale_y_continuous(expand = c(0, 0))+ # x axis is getting cut-off
#   theme(legend.position="none",legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=8),
#         legend.justification = "left", legend.title = element_blank(),  plot.margin = margin(10, 10, 10, 10))+
#   annotate('text', x = 2014, y = 30, label =  paste("italic(n)", "~invaded~counties==", nrow(BLD.adj.iso)), parse = T, size=3, check_overlap = TRUE,  hjust = 0)+
#   geom_text(x = 2014, y = 26, label = "xx contiguous", parse = F, size=3,check_overlap = TRUE,  hjust = 0)+
#   geom_text(x = 2014, y = 22, label = "xx non-contiguous", parse = F, size=3,check_overlap = TRUE,  hjust = 0)
# 



### ### ### ### ### ### ### ###
# Adjacent vs isolated stacked bar
### ### ### ### ### ### ### ###
df_total <- isoadj_cty %>%
  group_by(BLD_YR) %>%
  summarise(total = sum(counties_invaded))

inv.isoadj_graph <- ggplot(isoadj_cty, aes(x=BLD_YR, y=counties_invaded, fill=STAT)) +
  # Create stacked bars
  geom_bar(stat="identity", position="stack") +
  # Add totals on top
  geom_text(data=df_total, 
            aes(x=BLD_YR , y=total, label=total, fill=NULL), 
            vjust=-0.5, 
            size=2, 
            fontface="bold") +
  scale_fill_manual(values=c(adj_col, iso_col),labels = c("Contiguous", "Non-contiguous"))+
  scale_x_continuous( breaks = seq(2012, 2025, by = 2))+ # x axis is getting cut-off
  scale_y_continuous(limits= c(0, 75), expand = c(0, 0))+ # x axis is getting cut-off
  # Formatting: Colors, titles, and theme
  theme_clean+
  theme(legend.position="none",legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=8),
        legend.justification = "left", legend.title = element_blank(),  plot.margin = margin(10, 10, 10, 10))+
  annotate('text', x = 2013, y = 60, label =  paste("italic(n)", "~invaded~counties==", nrow(BLD.adj.iso)), parse = T, size=3, check_overlap = TRUE,  hjust = 0)+
  geom_text(x = 2013, y = 55, label = "219 contiguous", parse = F, size=3,check_overlap = TRUE,  hjust = 0)+
  geom_text(x = 2013, y = 50, label = "82 non-contiguous", parse = F, size=3,check_overlap = TRUE,  hjust = 0)+
  labs(
    x="Year",
    y="Invaded counties",
    fill="")



### ### ### ### ### ### ### ###
# Isolated jump distances
### ### ### ### ### ### ### ###
iso_only <- BLD.adj.iso[which(BLD.adj.iso$STAT == "iso"),]
summary(iso_only$D_BDY_km)
round(mean(iso_only$D_BDY_km),0)
round(stderr(iso_only$D_BDY_km),0)
round(median(iso_only$D_BDY_km),0)
iso_only[order(iso_only$D_BDY_km, decreasing = T),]
quantile(iso_only$D_BDY_km, 0.99)

inv.isoadj_hist <-
  ggplot(iso_only, aes(x=D_BDY_km)) +
  geom_histogram(position="identity", col="black", fill=iso_col, boundary=0, bins=60)+
  #geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
  #           linetype="dashed")+
  labs(title="",x="Distance jumped (km)", y = "Count")+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0), limits= c(0,600))+ # x axis is getting cut-off
  scale_y_continuous(expand = c(0, 0), limits= c(0,15))+ # x axis is getting cut-off
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),  plot.margin = margin(10, 10, 10, 10))+
  geom_text(x = 200, y = 10, label = "Mean ± SE = 114 ± 10 km", parse = F, size=3,check_overlap = TRUE,  hjust = 0)+
  geom_text(x = 200, y = 8.8, label = "Median = 78 km", parse = F, size=3, check_overlap = TRUE,  hjust = 0)

resize.win(3.30709*2,7)

# make the patchwork graph
figure1_export <- (invaded_area_continuous | gg_inset_map1_isoadj)/
  (inv.isoadj_graph | inv.isoadj_hist) +
  plot_annotation(tag_levels = 'A')  # Places plots side-by-side
#ggsave("figures/figure1.tiff", plot = figure1_export, device = "pdf", width = 3.30709*2, height = 7, units = "in", dpi = 300)


figure1_export <- ggarrange(invaded_area_continuous, gg_inset_map1_isoadj, inv.isoadj_graph, inv.isoadj_hist, 
                            ncol = 2, nrow = 2,
                            labels = c("A", "B", "C", "D"))

# tiff("figures/figure1.tiff", width = 3.30709*2, height = 7, units = "in", res = 300, compression = "lzw")
# print(figure1_export)
# dev.off()
# beep(1)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 









# Effective range radius -------------------------------------------------------
#
#
#---
# calculate the area of each county
contigUSCAN_ALBERS_BLD$county_area <- st_area(contigUSCAN_ALBERS_BLD)

# create a data frame to be populated with year and area invaded 
annual_sprd <- data.frame(year=2012:max_year_of_invasion, area_invaded_km2=NA)
years.to.plot <- sort(unique(contigUSCAN_ALBERS_BLD$BLD_YR))
i <- years.to.plot[1]
# for loop that populates the annual_sprd data frame
for(i in years.to.plot){
  # all invaded counties as of year i
  curr.invaded <- contigUSCAN_ALBERS_BLD[which(contigUSCAN_ALBERS_BLD$BLD_YR %in% 2012:i),]
  # input current year
  annual_sprd[which(annual_sprd$year %in% i), "year"] <- i
  # input area invaded as of year i
  annual_sprd[which(annual_sprd$year==i), "area_invaded_km2"] <- sum(st_area(curr.invaded)/1e6) # convert from square m to square km
}

# equation for estimating radius of a circle (i.e., rearrange A = pi * r^2)
annual_sprd$radius <- sqrt(annual_sprd$area_invaded_km2/pi)
plot(radius~year, data=annual_sprd)

# treating area as circle
fit_sprd1 <- lm(radius ~ year, data=annual_sprd)
summary(fit_sprd1)
summary(fit_sprd1)
round(summary(fit_sprd1)$coef,2)
round(summary(fit_sprd1)$coef,2)[2,3]^2


fit_sprd1_NL <- lm(radius ~ year + I(year^2), data=annual_sprd)
summary(fit_sprd1_NL)
summary(fit_sprd1_NL)
round(summary(fit_sprd1_NL)$coef,2)
round(summary(fit_sprd1_NL)$coef,2)[2,3]^2



set.seed(12)

# breakpoint regression (broken stick model fit to data)
# determine if there is an inflection point
segmented.mod <- segmented(fit_sprd1, control = seg.control(display = FALSE))
summary(segmented.mod)

#
annual_sprd_BREAK_BEFORE <- annual_sprd[which(annual_sprd$year <= 2018),]
fit_rangeexpansion_BEFORE <- lm(radius ~ year, data=annual_sprd_BREAK_BEFORE)
summary(fit_rangeexpansion_BEFORE)
round(summary(fit_rangeexpansion_BEFORE)$coef,2)

#
annual_sprd_BREAK_AFTER <- annual_sprd[which(annual_sprd$year > 2018),]
fit_rangeexpansion_AFTER <- lm(radius ~ year, data=annual_sprd_BREAK_AFTER)
summary(fit_rangeexpansion_AFTER)
round(summary(fit_rangeexpansion_AFTER)$coef,2)

#
# graphing and analysis
# need this for adding equation information to graph
lm_eqn <- function(df){
  m <- lm(radius ~ year, df); # NOTE: response and predictor have to be specified here
  eq <- substitute(italic(y) == a + b * italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$adj.r.squared, digits = 2)))
  as.character(as.expression(eq));
}
## graphing --------------------------------------------------------------------
y_initial_ERR <- 430 # y location where text will be displayed on graph 
color_dots <- "#5A7DBA" # my_BLD_colors[11]# color of the points
min_x <- 2013 # minimum x value 

adjustment_btw_text <- 40
seg_g_txt_sz <- 2.5
dr_g_txt_sz <- 2.5
x_location_txt <- 2013.5
# graph for analysis of full time series
fig_ERR <- ggplot(data=annual_sprd, aes(x=year, y=radius)) +
  ylab("Radius (km)")+
  xlab("Year")+ theme_bw()+
  
  scale_x_continuous(breaks = seq(2012,2026,2), limits=c(2012,2027), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,500,100), limits=c(0,500), expand = c(0,0)) +
  
  geom_point(aes(x = (year), y = radius), col=color_dots, size=3)+
  stat_smooth(data = annual_sprd, method = "lm", col = "black", se=F)+
  geom_text(x = min_x+1, y = y_initial_ERR, label = lm_eqn(annual_sprd), parse = TRUE, size=dr_g_txt_sz, check_overlap = TRUE, hjust = 0)+
  #geom_segment(aes(x = 2004, y = y_initial_ERR, xend =  2004+0.8, yend = y_initial_ERR), colour = "black",size=0.8)+
  
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# graph for segmented regression

summary(annual_sprd)
fig_ERR_SEGMENT <- ggplot(data=annual_sprd, aes(x=year, y=radius)) +
  ylab("Radius (km)")+
  xlab("Year")+ theme_bw()+
  
  scale_x_continuous(breaks = seq(2012,2025,2), limits=c(2012,2026), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,500,100), limits=c(0,500), expand = c(0,0)) +
  
  geom_point(aes(x = (year), y = radius), col=color_dots, size=3)+
  
  stat_smooth(data = annual_sprd, method = "lm", col = "black", se=F)+
  geom_text(x = min_x+1.4, y = y_initial_ERR+adjustment_btw_text, label = lm_eqn(annual_sprd), parse = TRUE, size=seg_g_txt_sz, check_overlap = TRUE, hjust = 0)+
  annotate("segment", x = min_x, y = y_initial_ERR+adjustment_btw_text, xend =  min_x+1.2, yend = y_initial_ERR+adjustment_btw_text, colour = "black", size=0.8)+
  
  
  stat_smooth(data = annual_sprd_BREAK_BEFORE, method = "lm", col = "black", se=F, linetype="dotted")+
  geom_text(x = min_x+1.4, y = y_initial_ERR, label = lm_eqn(annual_sprd_BREAK_BEFORE), parse = TRUE, size=seg_g_txt_sz, check_overlap = TRUE, hjust = 0)+
  annotate("segment", x = min_x, y = y_initial_ERR, xend =  min_x+1, yend = y_initial_ERR, colour = "black", size=1.2, linetype="dotted")+
  
  stat_smooth(data = annual_sprd_BREAK_AFTER, method = "lm", col = "black", se=F, linetype="dashed")+
  geom_text(x = min_x+1.4, y = y_initial_ERR-adjustment_btw_text, label = lm_eqn(annual_sprd_BREAK_AFTER), parse = TRUE, size=seg_g_txt_sz,check_overlap = TRUE,  hjust = 0)+
  annotate("segment", x = min_x, y = y_initial_ERR-adjustment_btw_text, xend =  min_x+1.2, yend = y_initial_ERR-adjustment_btw_text, colour = "black",size=0.8, linetype="dashed")+
  
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"));fig_ERR_SEGMENT





lm_eqn2 <- function(df){
  m <- lm(radius ~ year + I(year^2), df); # NOTE: response and predictor have to be specified here
  eq <- substitute(italic(y) == a - b * italic(x)* + c * italic(x)^2*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(abs(unname(coef(m)[2])), digits = 4),
                        c = format(unname(coef(m)[3]), digits = 3),
                        r2 = format(summary(m)$adj.r.squared, digits = 2)))
  as.character(as.expression(eq));
}  


fig_ERR2 <- ggplot(data=annual_sprd, aes(x=year, y=radius)) +
  ylab("Radius (km)")+
  xlab("Year")+ theme_bw()+
  
  scale_x_continuous(breaks = seq(2012,2026,2), limits=c(2012,2027), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,500,100), limits=c(0,500), expand = c(0,0)) +
  
  geom_point(aes(x = (year), y = radius), col=color_dots, size=3)+
  stat_smooth(data = annual_sprd, method = "lm", formula = y ~ x + I(x^2), col = "black", se=F)+
  geom_text(x = x_location_txt-1, y = y_initial_ERR, label = lm_eqn2(annual_sprd), parse = TRUE, size=dr_g_txt_sz, check_overlap = TRUE, hjust = 0)+
  #annotate("segment", x = min_x, y = y_initial_ERR, xend =  min_x+1, yend = y_initial_ERR, colour = "black", size=0.8)+
  
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 






# Distance regression ----------------------------------------------------------
#
#
#---
# convert shapefile to a dataframe
sprd_df <-as.data.frame(contigUSCAN_ALBERS_BLD)
# remove non-invaded counties
sprd_df <- sprd_df[!is.na(sprd_df$BLD_YR),]
# remove the first year of distances from analysis
sprd_df <- sprd_df[which(sprd_df$BLD_YR > 2012),] # remove first year of spread
plot(DtoDL_km~BLD_YR, data=sprd_df)

# fit model, regress distance to DL on year of invasion
fit_sprd2 <- lm(DtoDL_km  ~ BLD_YR, data=sprd_df)
summary(fit_sprd2)
round(summary(fit_sprd2)$coef,2)
confint(fit_sprd2)

fit_sprd2_NL <- lm(DtoDL_km ~ BLD_YR + I(BLD_YR^2), data=sprd_df)
summary(fit_sprd2_NL)

fit_sprdln_NL <- lm(log(DtoDL_km) ~ BLD_YR, data=sprd_df)
summary(fit_sprd2_NL)

fit_sprdln2_NL <- lm(log(DtoDL_km) ~ BLD_YR + I(BLD_YR^2), data=sprd_df)
summary(fit_sprdln2_NL)
round(summary(fit_sprdln2_NL)$coef,2)
round(summary(fit_sprdln2_NL)$coef,4)

fit_sprd_logx_NL <- lm(DtoDL_km ~ log(BLD_YR), data=sprd_df)
summary(fit_sprd_logx_NL)



#
## graphing ---------------------------------------------------------------------
# need this for adding equation information to graph
lm_eqn_DR <- function(df){# NOTE: response and predictor have to be specified here
  m <- lm(DtoDL_km ~ BLD_YR, sprd_df);
  eq <- substitute(italic(y) == a + b * italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$adj.r.squared, digits = 2)))
  as.character(as.expression(eq));
}

# colors for graph
DR_col <- color_dots

summary(sprd_df$DtoDL_km)
summary(sprd_df$BLD_YR)


fig_DR <- ggplot(data=sprd_df, aes(x=BLD_YR, y=DtoDL_km)) +
  ylab("Distance (km)")+
  xlab("Year")+ theme_bw()+
  geom_jitter(aes(x = BLD_YR, y = DtoDL_km), position = position_jitter(width = 0.1, height = 0), alpha = 0.7, color=DR_col) +
  
  stat_smooth(data = sprd_df, fullrange=T, method = "lm", col = "black", se=F,  xseq = seq(2013,2025.5, length=100))+
  
  scale_x_continuous(breaks = seq(2013,2025,2), limits=c(2012.8,2025.5), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,1200,300), limits=c(0,1200), expand = c(0,0)) +
  
  #geom_text(x = 2005, y = 1200, label = expression("Spread = 43 ± 4 km/yr (95% CI: 35-51)"), size=6,check_overlap = TRUE, hjust = 0)+
  geom_text(x = x_location_txt, y = 1150, label = lm_eqn_DR(sprd_df), parse = TRUE, size=dr_g_txt_sz, check_overlap = TRUE, hjust = 0)+
  #geom_text(x = 2005, y = 875, label = expression(italic(F)["1,253"]*"= 144.9, "*italic(p)*" < 0.0001"),  size=6,check_overlap = TRUE, hjust = 0)+
  
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "mm"))





lm_eqn_DR2 <- function(df){
  m <- lm(log(DtoDL_km) ~ BLD_YR + I(BLD_YR^2), df); # NOTE: response and predictor have to be specified here
  eq <- substitute(atop(italic(log(y)) == a + b * italic(x)* - c * italic(x)^2, italic(r)^2~"="~r2 * phantom(~"                               ")), 
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(abs(unname(coef(m)[2])), digits = 4),
                        c = format(abs(unname(coef(m)[3])), digits = 3),
                        r2 = format(summary(m)$adj.r.squared, digits = 2)))
  as.character(as.expression(eq));
}  



grid <- data.frame(BLD_YR = seq(min(sprd_df$BLD_YR), max(sprd_df$BLD_YR), length.out = 1000))
grid$y_pred <- exp(predict(fit_sprdln2_NL, newdata = grid))

fig_DR2 <- ggplot(data=sprd_df, aes(x=BLD_YR, y=DtoDL_km)) +
  ylab("Distance (km)")+
  xlab("Year")+ theme_bw()+
  geom_jitter(aes(x = BLD_YR, y = DtoDL_km), position = position_jitter(width = 0.1, height = 0), alpha = 0.7, color=DR_col) +
  
  #  stat_smooth(data = sprd_df, fullrange=T, method = "lm", formula = log(y) ~ x + I(x^2), col = "black", se=F,  xseq = seq(2013,2025.5, length=100))+
  geom_line(data = grid, aes(y = y_pred), color = "black", linewidth = 1) +
  scale_x_continuous(breaks = seq(2013,2025,2), limits=c(2012.8,2025.5), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,1200,300), limits=c(0,1200), expand = c(0,0)) +
  
  #geom_text(x = 2005, y = 1200, label = expression("Spread = 43 ± 4 km/yr (95% CI: 35-51)"), size=6,check_overlap = TRUE, hjust = 0)+
  geom_text(x = x_location_txt-0.5, y = 1050, label = lm_eqn_DR2(sprd_df), parse = TRUE, size=dr_g_txt_sz,check_overlap = TRUE, hjust = 0)+
  #geom_text(x = 2005, y = 875, label = expression(italic(F)["1,253"]*"= 144.9, "*italic(p)*" < 0.0001"),  size=6,check_overlap = TRUE, hjust = 0)+
  
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "mm"))




resize.win(6.85,4)
fig_DR # sometime you get the below warning message if jittering pushes a point outside the graph
# "Removed 1 row containing missing values or values outside the scale range (`geom_point()`). "
# change jittering values or x axis limits if so
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 






# Boundary displacement --------------------------------------------------------
#
#
#---

# create the shapefile/sf for the wedges 
# custum functions were taken from stack overflow
# https://stackoverflow.com/questions/59328707/how-do-i-partition-a-circle-into-equal-multipolygon-slices-with-sf-and-r
st_wedge <- function(x,y,r,start,width,n=20){
  theta = seq(start, start+width, length=n)
  xarc = x + r*sin(theta)
  yarc = y + r*cos(theta)
  xc = c(x, xarc, x)
  yc = c(y, yarc, y)
  st_polygon(list(cbind(xc,yc)))   
}

st_wedges <- function(x, y, r, nsegs){
  width = (2*pi)/nsegs
  starts = (1:nsegs)*width
  polys = lapply(starts, function(s){st_wedge(x,y,r,s,width)})
  mpoly = st_cast(do.call(st_sfc, polys), "MULTIPOLYGON")
  mpoly
}

# wedge creation
# set up the data frame, first with years of the study
sprd_rad <- data.frame(year=2012:max_year_of_invasion)
# set the radii size
radii_n <- seq(0,360-22.5,22.5) # emanating every 22.5°
# create a data frame to populated (spread per year per radii)
sprd_rad[,paste("radii",radii_n, sep=".")] <- NA

summary(contigUSCAN_ALBERS_BLD$DtoDL_km) # farthest spread in <1,200km
wedge_length <- 1200 # km in which radii will extend in each direction
discx <- st_coordinates(pts_LakeOH_ALBERS)[1]
discy <- st_coordinates(pts_LakeOH_ALBERS)[2]
wedges = st_wedges(discx,discy,wedge_length*1000,length(radii_n))
st_crs(wedges) <- crs(contigUSCAN_ALBERS_BLD)
plot(wedges) # check em out
plot(wedges[1], add=T, col="red") # which edge is wedge 1
# make wedge 1 the first north and east wedge
wedges <- wedges[c(length(radii_n),1:length(radii_n)-1)]
plot(wedges[1], add=T, col="blue") # confirm it worked!

## for loop for calculating spread per wedge -----------------------------------
i <- 2012 # start in the first invasion year
for(i in years.to.plot){
  
  # get all invaded counties as of year i
  curr.invaded.all <- contigUSCAN_ALBERS_BLD[which(contigUSCAN_ALBERS_BLD$BLD_YR %in% min(na.omit(contigUSCAN_ALBERS_BLD$BLD_YR)):i),]
  
  # for each wedge, get the counties in which their CENTROIDS occur inside a wedge (to avoid county membership in multiple wedges)
  # this returns the row ID for each county in curr.invaded.all
  int = suppressWarnings(st_intersects(wedges, st_centroid(curr.invaded.all)))
  
  # for each wedge 
  for(r in 1:nrow(int)){
    # r <- 8
    curr_wedge_polys <- int[r][[1]]
    curr.invaded.wedge <- curr.invaded.all[curr_wedge_polys,]
    counties_to_evaluate <- curr.invaded.wedge[curr.invaded.wedge$CTY != "Lake", ]
    
    if(nrow(counties_to_evaluate)>0){ # if there are any invaded counties, calculate the farthest centroid reached
      max_val  <- max(counties_to_evaluate$DtoDL_km)} else{
        max_val <- 0}
    
    # create a data frame with the distances spread per radii per year
    if(r == 1){
      fin_dists <- max_val}else{
        fin_dists <- append(fin_dists,max_val)}
  }
  
  sprd_rad[which(sprd_rad$year == i), paste("radii",radii_n, sep=".") ] <-  fin_dists
  
}


# convert from farthest distance spread per wedge to distance spread between consecutive years
# that is, we need to know how far it moved each year, not how far away it was from the
# discovery location each year
i=2
for(i in 2:ncol(sprd_rad)){
  curr.col <- sprd_rad[,i]
  sprd_rad[,i] <- c(NA,((diff(curr.col))))
}
# remove 2012 (first year of detection)
sprd_rad <- na.omit(sprd_rad)
summary(sprd_rad)





## for loop to remove years after which the farthest point in  interval was reached ----
wedge_id <- data.frame(wedge_n = radii_n, year=NA)
# all counties per wedge (by wedge)
int_all = st_intersects(wedges, st_centroid(contigUSCAN_ALBERS_BLD))

# distance between each county and discovery location
contigUSCAN_ALBERS_BLD$dist_from_disc <- suppressWarnings(pointDistance(st_centroid(contigUSCAN_ALBERS_BLD),st_centroid(pts_LakeOH_ALBERS),lonlat=F))/1000 # km

for(r in 1:nrow(int_all)){
  # r <- 13 # radii_n[13]
  curr_wedge_polys <- int_all[r][[1]]
  curr.wedge <- contigUSCAN_ALBERS_BLD[curr_wedge_polys,] # polygons in current wedge
  
  # calculate the distances from each county to the discovery location
  if(nrow(curr.wedge)>0){
    max_val_cty <- curr.wedge[which(curr.wedge$dist_from_disc == max(curr.wedge$dist_from_disc)),]
    FIPS_FAR <- max_val_cty$FIPS
  }else{FIPS_FAR <- "00001"}
  
  # calculate the distances from each INVADED coutny to discovery location
  invaded_cty_wedge <- curr.wedge[!is.na(curr.wedge$D_BDY_km),]
  if(nrow(invaded_cty_wedge)>0){
    max_val_invaded_cty <- invaded_cty_wedge[which(invaded_cty_wedge$DtoDL_km == max(invaded_cty_wedge$DtoDL_km)),]
    FIPS_INV_FAR <- max_val_invaded_cty$FIPS
  }else{FIPS_INV_FAR <- "00002"}
  if(FIPS_INV_FAR == FIPS_FAR){ # if the farthest county in a wedge is invaded, find the year in which that happened
    wedge_id[r,"year"] <- max_val_invaded_cty$BLD_YR 
  } else{
    wedge_id[r,"year"] <- max_year_of_invasion
  }
}

#plot(st_geometry(contigUSCAN_ALBERS_BLD))
#plot(st_geometry(curr.wedge), col="red", add=T)

## graph of stacked bar chart - spread per bearing ----
# organize data for graphingh
sprd.rads.BLD <- sprd_rad %>% gather(key = "bearing", value="sprd_increment", -year)
sprd.rads.BLD$bearing <- as.numeric(substr(sprd.rads.BLD$bearing,7,nchar(sprd.rads.BLD$bearing)))
bearings <- unique(sprd.rads.BLD$bearing) 

# for loop to remove years after which spread reached the end of the bearing
for(b in bearings){
  # b <- bearings[1]
  curr_bearing <- sprd.rads.BLD[which(sprd.rads.BLD$bearing == b),]
  curr_wedge_year_max <- wedge_id[which(wedge_id$wedge_n %in% b), "year"]
  curr_bearing_trunc <- curr_bearing[which(curr_bearing$year <= curr_wedge_year_max),]
  if(b == bearings[1]){
    sprd.rads.BLD_proc <- curr_bearing_trunc} else{
      sprd.rads.BLD_proc <- rbind.data.frame(sprd.rads.BLD_proc,curr_bearing_trunc)}
}
sprd.rads.BLD.g <- sprd.rads.BLD_proc %>% group_by(bearing) %>% summarise(mn_sprd_km = mean(sprd_increment), max_sprd_km = sum(sprd_increment), SE = stderr(sprd_increment))
summary(sprd.rads.BLD.g)
stderr(sprd.rads.BLD.g$mn_sprd_km)

#https://www.r-graph-gallery.com/136-stacked-area-chart
# ggplot(sprd.rads.BLD, aes(x = bearing, y = sprd_increment, fill=factor(year))) +
# geom_area(stat="identity") +
# theme_bw()+
# xlab("Bearing")+
# ylab("Spread distance (km)")+ 
# scale_x_continuous(breaks = seq(0, 360, 22.5*2), limits = c(0,360), expand=c(0,0))+
# scale_y_continuous(breaks = seq(0, 1200, 300), limits = c(0,1200), expand=c(0,0))+
# theme_bw() +
# coord_cartesian(clip="off")+
# scale_fill_manual(values=my_BLD_colors, name="")+
# theme(legend.key.size = unit(0.3, 'cm'), legend.text = element_text(size=7), legend.justification = "left")+
# theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#       panel.background = element_blank(), axis.line = element_line(colour = "black"))

### ### ### ### ### ### ### ###
# Bar chart by bearings/year plot
### ### ### ### ### ### ### ###
col_bearings <- ggplot(sprd.rads.BLD, aes(x = bearing, y = sprd_increment, fill=factor(year))) +
  geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = my_BLD_colors)+ 
  # name="", na.value = "white", na.translate = F)+   theme_bw()+
  xlab("Bearing")+
  ylab("Spread (km)")+
  theme(legend.position = c(0.8,0.6), legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6), legend.justification = "left")+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  guides(fill = guide_legend(title = NULL))




# create an insect circle for clarifying wedge location
wedges_16 = st_wedges(discx,discy,wedge_length*1000,16)
st_crs(wedges_16) <- crs(contigUSCAN_ALBERS_BLD)
wedges_4 = st_wedges(discx,discy,(wedge_length+400)*1000,12)
st_crs(wedges_4) <- crs(contigUSCAN_ALBERS_BLD)
#
sequence_labels_circle <- seq(0,337.5,22.5*2) # needs to reflect numer of wedges
label_x <- (wedge_length+300)*1000*sin((pi*sequence_labels_circle)/180)+discx
label_y <- (wedge_length+300)*1000*cos((pi*sequence_labels_circle)/180)+discy
# 
# my_BLD_colors_spec <- brewer.pal(11,"Spectral")
# inset_circle <- ggplot() +
#   geom_sf(data = wedges_4, fill="transparent", col="transparent") + # polygons filled based on the density value
#   #geom_sf(data = BLD_ALBERS_cropped, aes(fill = BLD_YR), color="light gray") + #polygons filled based on the density value
#   theme_bw()+theme_void()+
#   geom_sf(data=states_provinces_albers_crop, fill="transparent", color="black", lwd=0.1)+
#   scale_fill_gradientn(colors = my_BLD_colors_spec, name="", na.value = "white",   guide = "none")+
#   theme(legend.position="none",legend.key.size = unit(0.5, 'cm'),
#         legend.text = element_text(size=8), legend.justification = "left")+
#   geom_sf(data=pts_LakeOH_ALBERS, color = "black", size = 3)+
#   #geom_sf(data = wedges, fill="transparent", col="light gray") + # polygons filled based on the density value
#   theme_bw()+theme_void()+
#   geom_sf(data =wedges_16, fill="transparent", col="black", size=0.5)+# polygons filled based on the density value
#   annotate("text", label=sequence_labels_circle, x=label_x, y=label_y, size=3)
# 
# resize.win(6.85,4)
# gg_inset_bearings = ggdraw() +
#   draw_plot(col_bearings) +
#   #draw_plot(inset_circle, x = 0.55, y = 65, width = 0.35, height = 0.35)
# gg_inset_bearings
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 




## graph of spread along 22.5° bearing wedges ------------------------------------
# https://stackoverflow.com/questions/43033628/circular-histogram-in-ggplot2-with-even-spacing-of-bars-and-no-extra-lines
sprd.rads.avg <- sprd.rads.BLD_proc %>% group_by(bearing) %>% summarise(mean_sprd = mean(sprd_increment))
summary(sprd.rads.avg)

sprd.rads.avg[order(sprd.rads.avg$mean_sprd, decreasing=T),]
BD_col <- color_dots
resize.win(6.85,6.85)
unit_adj <- 1

bearing_histogram <- ggplot(sprd.rads.avg,
                            aes(x = bearing, y = mean_sprd)) +
  geom_col(width = 22.5, fill = BD_col, color = "black") +
  coord_polar(start = -pi/12+0.05) + # change start value if you want a different orientation
  scale_x_continuous(breaks = seq(0,270, 45)) +
  scale_y_continuous(breaks = seq(0, 100, 20))+ # x axis is getting cut-off
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.ontop = TRUE, # change to FALSE for grid lines below the wind rose
        panel.background = element_blank())+
  annotate("text", x=30, y=seq(20,100,20), label= seq(20,100,20), size=3, hjust=1.2, vjust=1, fontface = 2)+
  #annotate("text", x=seq(0,270, 45), y=120, label= seq(0,270, 45), size=3, hjust=1, vjust=1, fontface = 2)+
  
  theme(panel.grid.minor = element_blank(),
        #panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(colour="black", linetype="dashed"),
        panel.grid.major.y = element_line(colour="black", linetype="dashed"),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text=element_text(size=12,face="bold"),
        axis.title=element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"))


# NOTE that spread distances jump further here (with BD method) because you are
# just measuring jumps along radii, and not jumps from nearest county
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 





# Invaded area and US map with wedges --------------------------------------------
#
#
#---
wedges_16_padj <- wedges_16
invaded_area_discrete <- ggplot() +
  geom_sf(data = contigUSCAN_ALBERS_BLD, aes(fill = (BLD_YR)), color="transparent") + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data=states_provinces_albers, fill="transparent", color="black", lwd=1)+
  scale_fill_viridis(
    name="", na.value = "white", guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", direction = "horizontal"))+
  theme(legend.position=c(0.49,0.25),legend.key.height =  unit(0.2, 'cm'),
        legend.text = element_text(size=7), legend.justification = "left")+
  geom_sf(data=pts_LakeOH_ALBERS, color = "yellow", fill="black", size = 3, shape=21)+
  ggrepel::geom_text_repel(
    data = pts_LakeOH_ALBERS,
    aes(label = NAME, geometry = geometry),
    stat = "sf_coordinates",
    size=2.8,
    min.segment.length = 0,
    nudge_x = 0.5e6,
    nudge_y = 0.8e5,
    colour = "black",
    segment.colour = "black")+
  geom_sf(data = wedges_16_padj, fill="transparent", color="transparent") + #polygons filled based on the density value
  coord_sf(datum=st_crs(BLD_ALBERS_cropped))
invaded_area_discrete


label_x_1 <- (wedge_length+100)*1000*sin((pi*seq(0,337.5,22.5))/180)+discx
label_y_1 <- (wedge_length+100)*1000*cos((pi*seq(0,337.5,22.5))/180)+discy

wedges_p <- ggplot() +
  geom_sf(data = contigUSCAN_ALBERS_BLD, aes(fill = (BLD_YR)), color="transparent") + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data=states_provinces_albers, fill="transparent", color="black", lwd=1)+
  scale_fill_viridis(
    name="", na.value = "white", guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", direction = "horizontal"))+
  theme(legend.position="none")+
  geom_sf(data=states_provinces_albers, fill="transparent", color="light gray", lwd=1)+
  #geom_sf(data = wedges_padj, fill="transparent", color="gray") + #polygons filled based on the density value
  geom_sf(data = wedges_16_padj, fill="transparent", color="black", size=1) + #polygons filled based on the density value
  geom_sf(data=pts_LakeOH_ALBERS, color = "yellow", fill="black", size = 3, shape=21)+
  coord_sf(datum=st_crs(BLD_ALBERS_cropped))+
  annotate("text", label=seq(0,337.5,22.5), x=label_x_1, y=label_y_1, size=3, fontface=2)
#
wedges_p


resize.win(6.85,6.85)
figure2_export <- ggarrange(fig_ERR_SEGMENT, fig_ERR2,
                            fig_DR, fig_DR2, 
                            col_bearings, bearing_histogram, 
                            ncol = 2, nrow = 3,
                            labels = c("A", "B", "C", "D", "E", "F"))
figure2_export
# the WARNING is for the predicted regression line in panel A going below 0


# 2. Export to a file (supports pdf, png, jpeg, etc.)
# pdf("figures/figure2.pdf", width = 3.30709*2, height = 7)
# print(figure2_export)
# dev.off()
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 