#' TrendNPS: An R package for trend analysis for data from complex survey designs.
#'
#' @section Description:
#' The TrendNPS package provides trend analysis tools for complex survey designs.  
#' Please refer to the accompanying National Park Service Natural Resource Report (NRR) series report
#' Starcevich et al. (2017) for background on how and under what circumstances to use the trend analysis methods.
#' Example analyses are provided using the tools offered in thte TrendNPS R package.
#' 
#' @section TrendNPS functions:
#' TrendNPS_Cont, 
#' LinearizationVar_StRS, 
#' LinearizationVar, 
#' PWIGLS_ALL 
#'
#' @docType package
#' @name TrendNPS
NULL

#' NPS lake chemistry monitoring data for sites with at least two visits.  Data have been merged with GRTS sampling
#' design factors.
#'
#' @name LakeDataRep
#' @docType data
#' @usage data(LakeDataRep)
#' @format A data frame with 77 rows and 17 variables.
#' @keywords datasets
#' @references Heard, A.M., L.A. H. Starcevich, J.O. Sickman, M. Goldin Rose, and D.W. Schweizer. 2012. 
#' Sierra Nevada Network lake monitoring protocol. Natural Resource Report NPS/SIEN/NRR—2012/551. 
#' National Park Service, Fort Collins, Colorado.
#' 
#' A dataset containing the Sierra Nevada Network lake chemistry data from Sequoia-Kings Canyon and Yosemite National Parks. 
#' See Starcevich et al. (2017) for more information. The variables are as follows:
#'
#' \itemize{
#'   \item Site. Lake identifier.
#'   \item Park. National Park used as stratification variable. (SEKI = Sequoia-Kings Canyon, YOSE = Yosemite)
#'   \item Year. The observation year. (2008--2013)
#'   \item WYear. The observation year minus 2008. (0--5)
#'   \item xcoord. UTM-xxN easting. 
#'   \item ycoord. UTM-xxN northing.  
#'   \item mdcaty. Multi-density category representing cost class from the unequal probability sample draw.
#'   \item panel. Temporal revisit panel, or oversample when an oversample site was used. 
#'   \item wgt. Design weight from original GRTS sample draw.
#'   \item lake_area. Lake area in square meters.
#'   \item lake_perim. Lake perimeter in meters.
#'   \item Elevation_m. Lake elevation in meters.
#'   \item pH. pH measurement of lake. 
#'   \item Cl. Chloride measurement of lake. 
#'   \item EvalStatus. Evaluation status of lake relative to frame and nonresponse error.
#'   \item AdjWgt. Design weight adjusted for frame error and nonresponse error. 
#'   \item PanelWt. Panel weight, or the inverse of the panel inclusion probability. 
#' }
NULL

#' NPS lake chemistry monitoring data for all sites.  
#'
#' @name LakeData
#' @docType data
#' @usage data(LakeData)
#' @format A data frame with 130 rows and 6 variables.
#' @keywords datasets
#' @references Heard, A.M., L.A. H. Starcevich, J.O. Sickman, M. Goldin Rose, and D.W. Schweizer. 2012. 
#' Sierra Nevada Network lake monitoring protocol. Natural Resource Report NPS/SIEN/NRR—2012/551. 
#' National Park Service, Fort Collins, Colorado.
#' 
#' A dataset containing the Sierra Nevada Network lake chemistry data from Sequoia-Kings Canyon and Yosemite National Parks. 
#' See Starcevich et al. (2017) for more information. The variables are as follows:
#'
#' \itemize{
#'   \item Site. Lake identifier.
#'   \item Year. The observation year. (2008--2013)
#'   \item Park. National Park used as stratification variable. (SEKI = Sequoia-Kings Canyon, YOSE = Yosemite)
#'   \item Elevation_m. Lake elevation in meters.
#'   \item pH. pH measurement of lake. 
#'   \item Cl. Chloride measurement of lake. 
#' }
NULL

#' NPS lake chemistry GRTS sample.
#'
#' @name GRTSsample
#' @docType data
#' @usage data(GRTSsample)
#' @format A data frame with 461 rows and 10 variables.
#' @keywords datasets
#' @references Heard, A.M., L.A. H. Starcevich, J.O. Sickman, M. Goldin Rose, and D.W. Schweizer. 2012. 
#' Sierra Nevada Network lake monitoring protocol. Natural Resource Report NPS/SIEN/NRR—2012/551. 
#' National Park Service, Fort Collins, Colorado.
#'
#' Kincaid, T. M. and Olsen, A. R. 2015. spsurvey: Spatial Survey Design and Analysis. R package version 3.1. 
#' URL:http://www.epa.gov/nheerl/arm/.
#' 
#' The design file containing the Sierra Nevada Network lake chemistry GRTS sample draw. Sites were stratified by Park
#' and unequal probability sampling relative to a cost class (mdcaty) were drawn. The sample was drawn using the spsurvey
#' package (Kincaid and Olsen 2015). The variables are as follows:
#'
#' \itemize{
#'   \item siteID. Lake identifier.
#'   \item Park. National Park used as stratification variable. (SEKI = Sequoia-Kings Canyon, YOSE = Yosemite)
#'   \item WYear. The observation year minus 2008. (0--5)
#'   \item xcoord. UTM-xxN easting. 
#'   \item ycoord. UTM-xxN northing.  
#'   \item mdcaty. Multi-density category representing cost class from the unequal probability sample draw.
#'   \item panel. Temporal revisit panel, or oversample when an oversample site was used. 
#'   \item wgt. Design weight from original GRTS sample draw.
#'   \item lake_area. Lake area in square meters.
#'   \item lake_perim. Lake perimeter in meters.
#' }
NULL