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
#' @description  Downloads climate data from the NASA POWER data base. It is useful for situations when data is required from seeveral sites.
#' @param t Is a dataframe with longitud and latitude in decimal formats and a location ID.
#' @param lon Is the column index in t that contains the longitud in decimal format
#' @param lat Is the column index in t that contains the latitude in decimal format
#' @param loc.id Is column name of the location ids in the dataframe t
#' @param ... Arguments passed to the get_power function
#' @export

climdata = function(t = NULL, lon = NULL, lat = NULL, loc.id = NULL, ... ){

t$Loc_no = t[, loc.id]
t.split= split(t, t$Loc_no)
t.split = lapply(t.split, droplevels)

d = lapply(t.split, function(x) nasapower::get_power(lonlat = c(x[,lon], x[,lat]), ...))

d.d = reshape2::melt(d, id = c("LON", "LAT"))

tmp <- plyr::ddply(d.d, plyr::.(L1, variable), transform, newid = paste(L1, seq_along(variable)))


df = reshape2::dcast(tmp, L1 + LAT + LON + newid ~ variable, value.var = "value")
df = df[,-which(colnames(df) == "newid")]

df$YYYYMMDD = zoo::as.Date(df$YYYYMMDD)

df = df[order(df$L1, df$YEAR, df$MM, df$DD),]
colnames(df)[1] = loc.id
colnames(df)[8] = "Date"


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
#' @param df A dataframe of one row (one location) that contains the date from which meteorological variables (MVs) will be extracted.
#' @param climate.data A list of of locations. Each element of the list must contain the MVs.
#' @param total.days The number of days from which the MVs will be extracted. Is the number of days after the initial date that will be consider to extract the data.
#' @param by.periods The number of days which the MVs will be averaged over. i.e, it must be multiple and <= than total.days.
#' @param start.date The column name in df that contains the date from which the extraction will start
#' @param loc.id The column name in df that contains the location ids
#' @param Date The column name in climate.data that contains the date of the recorded MVs
#' @param mv.names The index of the columns that contain the MVs in climate.data
#' @param pp.name The name of the column that contains the rainfall. This is used since percipitation is summed over the by.periods, instead of averaged as the other MVs
#'
#' @export
#'

periodicextraction <- function(df = NULL, climate.data = NULL,
                               total.days = 120,
                               by.periods = 10,
                               start.date = "sowing.d",
                               loc.id = "Loc_no",
                               Date = "Date",
                               mv.names = 9:18,
                               pp.name = "pp") {


    period = total.days - 1
    df$start.date = df[,start.date]
    df$date.covs = df$start.date + period
    loc.no = as.character(df[, loc.id])
    tmp = climate.data[[loc.no]]
    tmp$Loc_no = df[,loc.id]
    tmp$Date = tmp[, Date]


    period.covs = subset(x = tmp, subset = tmp$Date >= df$start.date & Date <= df$date.covs & tmp$Loc_no == loc.no)
    period.covs = droplevels(period.covs)
    group.means = rep(1:(total.days/by.periods), each = by.periods)
    period.covs = cbind(period.covs, group.means)

    mvs = base::colnames(period.covs)[mv.names]
    mvs = mvs[!mvs %in% pp.name]
    frm = base::paste(mvs, collapse = " + ")
    frm = base::paste(frm, "group.means", sep = "~")
    f.mvs = stats::as.formula(frm)

    cov.means = doBy::summaryBy(formula = f.mvs, data = period.covs, FUN = mean, keep.names = TRUE, na.rm=TRUE)

    p.form = paste(pp.name, "group.means", sep = "~")
    p.form = stats::as.formula(p.form)

    cov.pp.sum = doBy::summaryBy(formula = p.form, data = period.covs,
                                 FUN = sum, keep.names = TRUE)
    f.env.covs = cbind(cov.means, cov.pp.sum[,2])
    f.conv.melt = reshape2::melt(f.env.covs, id.vars = "group.means")
    conv.cast = reshape2::dcast(f.conv.melt, 1 ~ variable + group.means)
    loc.env.covariables = cbind(df, conv.cast[, 2:ncol(conv.cast)])

    covs.final.name = colnames(loc.env.covariables)
    covs.final.name = gsub(pattern = "cov.pp.sum\\[\\, 2\\]", replacement = "PP", x = covs.final.name)

    colnames(loc.env.covariables) = covs.final.name
    loc.env.covariables

}
