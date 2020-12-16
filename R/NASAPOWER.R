############################################################
#                                                          #
#             DOWNLOADS AND FORMATS NASA POWER             #
#              DATA FROM A LIST OF LOCATIONS               #
#       lonlat argument of get_power function must        #
#      indicate c(x[,],x[x]) the longitud and latitud      #
#                                                          #
############################################################

#' Download NASA POWER data
#'
#' @description  Downloads climate data from the NASA POWER data base. This functions requires the nasapower library to be installed
#' @param t Is a dataframe with longitud and latitude in decimal formats and a "Loc_no" column as location ID.
#' @param lon Is the column index in t that contains the longitud in decimal format
#' @param lat Is the column index in t that contains the latitude in decimal format
#' @param new.names Is a vector of new names for the variables to be extracted. Used for simplicity in naming the variables.
#' @param ... Arguments passed to the get_power function
#' @export

climdata = function(t = NULL, lon = NULL, lat = NULL, new.names = NULL, ... ){

#library(nasapower)
t.split= split(t, t$Loc_no)
t.split = lapply(t.split, droplevels)
#lon & lat must indicate column index where longitud and latitud are in the dataframe
d = lapply(t.split, function(x) nasapower::get_power(lonlat = c(x[,lon], x[,lat]),...))

d.d = reshape2::melt(d, id = c("LON", "LAT"))

tmp <- plyr::ddply(d.d, plyr::.(L1, variable), transform, newid = paste(L1, seq_along(variable)))


df = reshape2::dcast(tmp, L1 + LAT + LON + newid ~ variable, value.var = "value")
df = df[,-which(colnames(df) == "newid")]

#library(zoo)
df$YYYYMMDD = zoo::as.Date(df$YYYYMMDD)

df = df[order(df$L1, df$YEAR, df$MM, df$DD),]
colnames(df)[1] = "Loc_no"
colnames(df)[8] = "Date"
#the new.names vector should be the new names of the climate variables
colnames(df)[9:ncol(df)] = new.names

df = merge(x = df, y = t, by = "Loc_no")
df
}




############################################################
#                                                          #
#     Function to summarize the climate date by groups      #
#      of days for a single location it requires the       #
#            locs.merged.climate object and the            #
#             locs.env.data object dates must              #
#  be in the appropriate format otherwise errors are shown  #
#                                                          #
############################################################

#' Summarize climdata outputs
#'
#' @description Function to summarize the climate date by groups of days for a single location. It requires the locs.merged.climate object which is the climate data.
#' Dates must be in the appropriate format otherwise errors are shown.
#' @param df Is a one row (one location) data frame with the date of sowing (sowing.d), or a date from which the variables will start to be summarized.
#' @param locs.merged.climate Is a list of locations with the climate data
#' @param days Is the number of days after the initial date that will be consider to extract the data
#' @param last.day Is the number of days that will be used to average the data. i.e, las.day must be <= than days and multiple of days
#' @export

periodicextraction = function(df = NULL, locs.merged.climate = NULL,
                               days=120, last.day=120){
    period = days-1 # "last.day" must be multiple of "days"
    df$date.covs = df$sowing.d + period
    locno = as.character(df$Loc_no)
    tmp = locs.merged.climate[[locno]]
    period.covs = subset(tmp, subset = Date >= df$sowing.d
                         & Date <= df$date.covs
                         & Loc_no == df$Loc_no )
    period.covs
    period.covs = droplevels(period.covs)
    group.means = rep(1:(days/last.day), each=last.day)
    period.covs = cbind(period.covs, group.means)
    #library(doBy)
    cov.means = doBy::summaryBy(rad + t.av + t.min + t.max + rel.h + dew.p + cdd
                          + clrsky + wind.speed
                          ~ group.means, data=period.covs, FUN = mean)
    cov.pp.sum = doBy::summaryBy(pp ~ group.means, data = period.covs, FUN = sum)
    f.env.covs = cbind(cov.means, cov.pp.sum$pp.sum)
    colnames(f.env.covs)[11] = "pp"
    #library(reshape2)
    f.conv.melt = reshape2::melt(f.env.covs, id.vars = "group.means")
    conv.cast = reshape2::dcast(f.conv.melt, 1 ~ variable+group.means)
    loc.env.covariables = cbind(df, conv.cast)
    loc.env.covariables = loc.env.covariables[, -c(15)]
    loc.env.covariables


}
