#' Calculating distance to road from GBIF data
#'
#' This function sets the paramters in the appropriate manner to be used by the
#' simulation function
#'
#' @param gbif_data, bioclim
#'
#' @return dataframe with distance to road a column
#' @export
distGBIFdata <- function(gbif_data, bioclim = TRUE){

message("Getting the countries from the GBIF data")
countries <- unique(gbif_data$countryCode)

if(bioclim == TRUE){
  message("Extracting the Bioclim data")
alt<- raster::getData('worldclim', var='alt', res=10)
bio <- raster::getData('worldclim', var='bio', res=10)
}

data_df <- list()
for(country.tag in 1: length(countries)){
  message(paste("Extracting states and bioclim data for", countries[country.tag]))
  all_data <- gbif_data%>%
    dplyr::filter(countryCode == countries[country.tag])
  States <- raster::getData("GADM", country = countries[country.tag], level = 2)
 # USborder <- rgeos::gUnaryUnion(States, id = States$ISO)
  #boundary_points <- fortify(USborder)
  sp_df <- sp::SpatialPointsDataFrame(cbind(as.numeric(all_data$decimalLongitude),
                                        as.numeric(all_data$decimalLatitude)),
                                  all_data,proj4string = CRS("+proj=longlat +datum=WGS84"))

  df_data <- raster::intersect(sp_df, States)
  data_with_admin <- df_data@data

  if(bioclim == TRUE){
  altitude <- raster::extract(alt, df_data, df=T)
  bio_data <- raster::extract(bio, df_data, df = T)
  data_df[[country.tag]] <- cbind(data_with_admin, altitude, bio_data)%>%
    data.frame()
  }else{
    data_df[[country.tag]] <- cbind(data_with_admin)%>%
      data.frame()
  }
}

 data_df <- do.call("rbind", data_df)

#save(gull_df, file = "gull_df.RData")

#######################
# Calculate distance to road
#load filtered data

distance_estimate <- function(i, data){

  # message("Obtaining the states from the data")
  states <- unique(data$level2Name)
  states <- gsub("\\s*\\([^\\)]+\\)","",as.character(states))
  #states <- states[states != "Rhein-Neckar-Kreis"]
  # define Belgrade's bounding box based on ggmap's bounding box
  message(paste("Extracting maps for", states[i]))
  bg_map <- ggmap::get_map(getbb(states[i]),
                           maptype = "toner-lite",
                           source = "stamen",
                           color="bw",
                           force=T)

  bg_bbox <- attr(bg_map, 'bb')

  message(paste("map coordinates for", states[i]))
  bbox <- c(xmin=bg_bbox[,2],
            ymin= bg_bbox[,1],
            xmax= bg_bbox[,4],
            ymax=bg_bbox[,3])


  #get states's paved roads
  message(paste("Getting roads for", states[i]))
  bg_way <- osmdata::opq(bbox = bbox, timeout = 1240, memsize = 1004857600) %>%
    osmdata::add_osm_feature(
      key = 'highway') %>%
    osmdata::osmdata_sf(quiet = T)

  message(paste("Extracting the lines", states[i]))

  bg_r <- bg_way$osm_lines

  message(paste("Formatting data from", states[i], " to Polygons"))
  pk <- data %>%
    dplyr::filter(level2Name == states[i]) %>%
    sf::st_as_sf(., coords = c("decimalLongitude", "decimalLatitude"),
             crs = st_crs(bg_r))

  message(paste("Calculating distance for ", states[i]))
  dists_0 <- sf::st_distance(pk,bg_r)
  dist2road <- apply(dists_0,1,min)

  message(paste("Returning distance for ", states[i]))
  ret <- data %>%
    dplyr::filter(level2Name == states[i]) %>%
    dplyr::mutate(distance = dist2road)

  return(ret)

}

states_length <- length(unique(data_df$level2Name))
data_with_distance <- lapply(as.list(seq(1,states_length)), function(x){
  ret <- distance_estimate(x, data_df)
})

return(data_with_distance)

}
