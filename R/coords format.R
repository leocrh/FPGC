############################################################
#                                                          #
#       FUNCTION TO FORMAT LONGITUD AND LATITUD FROM       #
#           UQ123 DATA AND PLOT A MAP IN GGPLOT            #
#                                                          #
############################################################

#' Changes the coordinates to decimal degrees from IWIS raw site files. It adds two columns to the input dataframe: lat and lon
#' @param coords a data frame from IWIS database that contains the columns Long_degrees, Long_minutes, Lat_degrees and Lat_minutes.
#' @examples
#' coords = format.coords(coords)

format.coords = function(coords) {
    coords$Long_degress = as.numeric(coords$Long_degress)
    coords$Lat_degress = as.numeric(coords$Lat_degress)
    coords$Long_minutes = as.numeric(coords$Long_minutes)
    coords$Lat_minutes = as.numeric(coords$Lat_minutes)

    coords$lat = coords$Lat_degress + coords$Lat_minutes/60
    coords$lon = coords$Long_degress + coords$Long_minutes/60
    coords$lat = ifelse(coords$Latitud == "S", coords$lat*(-1), coords$lat)
    coords$lon = ifelse(coords$Longitude == "W", coords$lon*(-1), coords$lon)
    coords
}
