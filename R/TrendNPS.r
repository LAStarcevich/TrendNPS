#' TrendNPS: A package for dealing with trending stuff for the NPS.
#'
#' The TrendNPS package does trend stuff.  It's pretty cool.  
#' 
#' @section TrendNPS functions:
#' TrendNPS_Cont
#' LinearizationVar_StRS
#' LinearizationVar
#' PWIGLS_ALL
#'
#' @docType package
#' @name TrendNPS
NULL

#' NPS Lakes and Stuff.  
#'
#' @name SEKIANC
#' @docType data
#' @usage data(SEKIANC)
#' @format A data frame with 80 rows and 16 variables.
#' @keywords datasets
#' @references Somewhere.
#' 
#' A dataset containing the sample observations of something cool. The variables are as follows:
#'
#' \itemize{
#'   \item Site. Lake identifier.
#'   \item Y. Something cool measured. (7.49--189.73)
#'   \item Year. The observation year. (2008--2013)
#'   \item WYear. The observation year minus 2008. (0--5)
#'   \item Panel. Something cool. (1--5)
#'   \item Lat. Northern hemisphere latitude. (36.35831--37.18324)
#'   \item Long. Western hemisphere longitude. (-118.805---118.268)
#'   \item InclProb. Inclusion probability. (0.041--0.357)
#'   \item PanelProb. Panel probability. (0.25--1)
#'   \item AdjInclProb. Adjusted inclusion probability. (0.01025--0.357)
#'   \item utmx_n83. UTM-xxN easting. (339565.2--386563.9)
#'   \item utmy_n83. UTM-xxN northing.  (4025111--4116651)
#'   \item lake_eleva. Lake elevation in meters? (2615--3750)
#'   \item raster_dec. Something cool. (2.547821--18.02151)
#'   \item inver_cost. Super cool. (0.055--0.392)
#'   \item Elev_Stratum. Classification as low- or high-elevation. (Low High)
#' }
NULL
