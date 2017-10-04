#' TrendNPS: An R package for trend analysis for data from complex survey
#' designs.
#'
#' @section Description: The TrendNPS package provides trend analysis tools for 
#'   complex survey designs. Please refer to the accompanying National Park 
#'   Service Natural Resource Report (NRR) series report Starcevich et al. 
#'   (2017) and Starcevich et al. (2018) for background on how and under what
#'   circumstances to use the trend analysis methods. Example analyses are
#'   provided using the tools offered in the TrendNPS R package.
#' 
#' @section TrendNPS functions:
#' \code{TrendNPS_Cont}, \cr
#' \code{TrendNPS_Count}, \cr
#' \code{TrendNPS_Binary}, \cr
#' \code{LinearizationVar_StRS}, \cr 
#' \code{LinearizationVar}, \cr 
#' \code{PWIGLS_ALL}
#'
#' @docType package
#' @name TrendNPS
NULL

#' NPS lake chemistry monitoring data for sites with at least two visits.  Data
#' have been merged with GRTS sampling design factors.
#'
#' @name LakeDataRep
#' @docType data
#' @usage data(LakeDataRep)
#' @format A data frame with 77 rows and 17 variables.
#' 
#' \describe{ A dataset containing the Sierra Nevada Network lake chemistry data
#' from Sequoia-Kings Canyon and Yosemite National Parks. See Starcevich et al.
#' (2017) for more information. The variables are as follows:
#'
#' \itemize{
#'   \item Site. Lake identifier.
#'   \item Park. National Park used as stratification variable. (\code{SEKI} = Sequoia-Kings Canyon, \code{YOSE} = Yosemite)
#'   \item Year. The observation year. (2008--2013)
#'   \item WYear. The observation year minus 2008. (0--5)
#'   \item xcoord. UTM easting. 
#'   \item ycoord. UTM northing.  
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
#' }
#' @keywords datasets
#' @references Heard, A. M., L. A. H. Starcevich, J. O. Sickman, M. Goldin Rose,
#'   and D. W. Schweizer. 2012. Sierra Nevada Network lake monitoring protocol.
#'   Natural Resource Report NPS/SIEN/NRR-2012/551. National Park Service, Fort
#'   Collins, Colorado.
#' 
NULL

#' NPS lake chemistry monitoring data for all sites.  
#'
#' @name LakeData_orig
#' @docType data
#' @usage data(LakeData_orig)
#' @format A data frame with 130 rows and 6 variables.
#' 
#' \describe{ A dataset containing the Sierra Nevada Network lake chemistry data
#' from Sequoia-Kings Canyon and Yosemite National Parks. See Starcevich et al.
#' (2017) for more information. The variables are as follows:
#'
#' \itemize{
#'   \item Site. Lake identifier.
#'   \item Year. The observation year. (2008--2013)
#'   \item Park. National Park used as stratification variable. (\code{SEKI} = Sequoia-Kings Canyon, \code{YOSE} = Yosemite)
#'   \item Elevation_m. Lake elevation in meters.
#'   \item pH. pH measurement of lake. 
#'   \item Cl. Chloride measurement of lake. 
#' }
#' }
#' @keywords datasets
#' @references Heard, A. M., L. A. H. Starcevich, J. O. Sickman, M. Goldin Rose,
#'   and D. W. Schweizer. 2012. Sierra Nevada Network lake monitoring protocol.
#'   Natural Resource Report NPS/SIEN/NRR-2012/551. National Park Service, Fort
#'   Collins, Colorado.
#' 
NULL

#' NPS lake chemistry GRTS sample.
#'
#' @name GRTSsample
#' @docType data
#' @usage data(GRTSsample)
#' @format A data frame with 461 rows and 10 variables.
#' 
#' \describe{The design file containing the Sierra Nevada Network lake chemistry
#' GRTS sample draw. Sites were stratified by Park and unequal probability
#' sampling relative to a cost class (\code{mdcaty}) were drawn. The sample was
#' drawn using the \code{spsurvey} package (Kincaid and Olsen 2015). The
#' variables are as follows:
#'
#' \itemize{
#'   \item siteID. Lake identifier.
#'   \item Park. National Park used as stratification variable. (\code{SEKI} = Sequoia-Kings Canyon, \code{YOSE} = Yosemite)
#'   \item WYear. The observation year minus 2008. (0--5)
#'   \item xcoord. UTM easting. 
#'   \item ycoord. UTM northing.  
#'   \item mdcaty. Multi-density category representing cost class from the unequal probability sample draw.
#'   \item panel. Temporal revisit panel, or oversample when an oversample site was used. 
#'   \item wgt. Design weight from original GRTS sample draw.
#'   \item lake_area. Lake area in square meters.
#'   \item lake_perim. Lake perimeter in meters.
#' }
#' }
#' @keywords datasets
#' @references Heard, A.M., L.A. H. Starcevich, J.O. Sickman, M. Goldin Rose,
#'   and D.W. Schweizer. 2012. Sierra Nevada Network lake monitoring protocol.
#'   Natural Resource Report NPS/SIEN/NRR-2012/551. National Park Service, Fort
#'   Collins, Colorado.
#'   
#'   Kincaid, T. M. and Olsen, A. R. 2015. spsurvey: Spatial Survey Design and
#'   Analysis. R package version 3.1. \url{http://www.epa.gov/nheerl/arm/}.
#' 
NULL

#' Seastar count data.
#'   
#' @name Seastar
#' @docType data
#' @usage data(Seastar)
#' @format A data frame with 65 rows and 8 variables.
#' 
#' \describe{The file contains counts of \emph{Dermasterias imbricata} (leather
#' seastar) collected between 2009 and 2015 across two parks, Katmai National
#' Park (KATM) and Kenia Fjords National Park (KEFJ) in the Southwest Alaskan
#' Network. The variables are as follows:
#' 
#' \itemize{
#'   \item Site. Sampling unit selected with GRTS sampling.
#'   \item Park. National Park. \code{KATM} = Katmai NP and \code{KEFJ} = Kenai Fjords NP.
#'   \item Year. Survey year formatted as a factor for random effects estimation.
#'   \item WYear. Scalar year covariate for trend analysis.
#'   \item Gradient. Site gradient (degrees).
#'   \item Long. Site longitude.   
#'   \item Lat. Site latitude.
#'   \item Count. Number of sea stars observed.
#' }
#' }
#' @keywords datasets
#' @author Leigh Ann Starcevich of Western EcoSystems Technology, Inc.
#' @references Alaska Ocean Observing System Gulf of Alaska Data Integration
#'   Portal. 2017. Nearshore: Intertidal Systems in Gulf of Alaska data
#'   available at: 
#'   \url{http://portal.aoos.org/gulf-of-alaska.php#metadata/53c052b6-8874-46d1-b40a-acc615a3879a/project/files}.
#'   
NULL

#' \emph{Acrosiphonia} binary cover data.
#' 
#' @name Cover_Acro
#' @docType data
#' @usage data(Cover_Acro)
#' @format A data frame with 1794 rows and 10 variables.
#' 
#' \describe{The file contains binary indicators of \emph{Acrosiphonia} cover
#' collected between 2009 and 2015 across two parks, Katmai National Park (KATM)
#' and Kenia Fjords National Park (KEFJ) in the Southwest Alaskan Network. Data
#' hierarchy includes Elevation class within Quadrat within Site.  The variables
#' are as follows:
#' 
#' \itemize{ 
#'   \item Site. Sampling unit selected with GRTS sampling. 
#'   \item Park. National Park. \code{KATM} = Katmai NP and \code{KEFJ} = Kenai Fjords NP.
#'   \item Quad. Quadrat within Site. 
#'   \item Elev. Elevation class within \code{Quad} Quadrat.
#'   \item Year. Survey year formatted as a factor for random effects estimation.
#'   \item WYear. Scalar year covariate for trend analysis.
#'   \item Gradient. Site gradient (degrees).
#'   \item Long. Site longitude.   
#'   \item Lat. Site latitude.
#' }
#' }
#' @keywords datasets
#' @author Leigh Ann Starcevich of Western EcoSystems Technology, Inc.
#' @references Alaska Ocean Observing System Gulf of Alaska Data Integration
#'   Portal. 2017. Nearshore: Intertidal Systems in Gulf of Alaska data
#'   available at: 
#'   \url{http://portal.aoos.org/gulf-of-alaska.php#metadata/53c052b6-8874-46d1-b40a-acc615a3879a/project/files}.
#'   
NULL 
