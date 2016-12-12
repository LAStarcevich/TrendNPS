#' @export
#' 
#' @title Adjust design weights for nonsampling error.
#'   
#' @description Using the EvalStatus field in the design file, adjust the design
#' weights for frame error and nonsampling error.  
#'   
#' #' @param dat	 Data frame containing columns at least for \code{Site}, \code{WYear}, 
#'   \code{Year}, the outcome of interest, and the survey evaluation status \code{EvalStatus}.
#'   See section "Options for variable \code{EvalStatus}" below.
#' 
#' @param popn	Takes the value of "Target" or "Sampled" to reflect inference to either the original
#' target population or to the subpopulation of sites that would not be subject to nonresponse, respectively. 
#'   
#' @param EvalStatus	 Status of the survey at each site relative to frame error and nonresponse error.
#' Note that this field is calculated separately for each year at a site.
#'   
#' @param wgt	Original design weight. The design weights for a given year should sum to the 
#' size of the population before accounting for any frame error. 
#'   
#' @return Returns a vector of adjusted design weights that should sum to the extent of the 
#' population of interest after accounting for frame error (popn = "Target") or frame and nonresponse
#' error (popn = "Sampled"). 
#'   
#' @details Adjusts the original sampling design weights for frame and nonresponse error.
#' See (Starcevich et al. 2016) for more information on weighting adjustments and populations of inference. 
#' See (Starcevich et al. 2017) for an example on applying the AdjWgt_FrameNR function to a monitoring data set. 
#' 
#' @section Options for variable \code{EvalStatus}:
#'   
#'   The \code{EvalStatus} field requires specific formatting for weighting adjustment.  Valid options include:
#'   
#'   \tabular{ll}{
#'   
#'   \code{"Target - Surveyed"} \tab Site is in target population and was successfully surveyed. \cr
#'   \code{"Target - Not Surveyed"}     \tab Site is in target population but was not successfully surveyed.  \cr
#'   \code{"Non-Target"}    \tab Site is not in target population, not surveyed.  \cr
#'   \code{"Not Evaluated"}     \tab Site was not evaluated as target or not, not surveyed. \cr
#'   
#'   }
#'   
#' @author Leigh Ann Starcevich of Western EcoSystems Technology, Inc.
#'   
#' @references Starcevich, L. A., G. DiDonato, T. McDonald, and J. Mitchell. 2016. A GRTS user’s manual for the 
#' SDrawNPS package: A graphical user interface for generalized random tessellation stratified (GRTS) sampling 
#' and estimation. Natural Resource Report NPS/PWRO/NRR—2016/1233. National Park Service, Fort Collins, Colorado.
#'   
#' Starcevich,L.A., T. McDonald, A. Chung-MacCoubrey, A. Heard, and J. Nesmith. 2017. Trend Estimation for Complex 
#' Survey Designs. Natural Resource Report NPS/xxxx/NRR-2017/xxxx. National Park Service, Fort Collins, Colorado. 
#'
#' @examples 
#' \dontrun{
#' #  ---- Read example data set.
#' 	LakeData$AdjWgt[LakeData$WYear == 0] = AdjWgt(dat=LakeData[LakeData$WYear == 0,],popn="Target",evalstatus = "EvalStatus",wgt = "wgt")
#' }
#' 
#' 
AdjWgt <- function (dat, popn, evalstatus = "EvalStatus", wgt = "wgt") {

# Format EvalStatus field
# remove all spaces and dashes
    dat[, evalstatus] <- gsub(" ", "", dat[, evalstatus], fixed = TRUE)
    dat[, evalstatus] <- gsub("-", "", dat[, evalstatus], fixed = TRUE)
# make evalstatus lower case
    dat[, evalstatus] <- tolower(dat[, evalstatus])
# standardize possible entries
    dat[, evalstatus] <- as.character(dat[, evalstatus])
    dat[, evalstatus][dat[, evalstatus] %in% c("targetsampled","targetsurveyed")] <- "targetsurveyed"
    dat[, evalstatus][dat[, evalstatus] %in% c("targetnotsampled","targetnotsurveyd")] <- "targetnotsurveyed"
    dat[, evalstatus][dat[, evalstatus] %in% c("nontarget","nottarget")] <- "nontarget"
    dat[, evalstatus][dat[, evalstatus] %in% c("noteval", "notevaluated")] <- "notevaluated"

# Error check EvalStatus field
	EvalVector = c("targetsurveyed","targetnotsurveyed","notevaluated","oversamp","nontarget")
	if(any(!dat[, evalstatus]%in%EvalVector)) print(paste("EvalStatus field should take only the following values:","Target - Surveyed,Target - Not Surveyed,Not Evaluated,OverSamp,Non-Target"))

# Calculate sample sizes of evaluated, target, responding, and oversample sites
    n.all <- nrow(dat[dat$panel != "OverSamp", ])
    dat.eval <- dat[dat[, evalstatus] != "notevaluated", ]
    n.eval <- nrow(dat.eval)
    dat.T <- dat.eval[dat.eval$EvalStatus != "nontarget", ]
    n.T <- nrow(dat.T)
    dat.R <- dat.T[dat.T$EvalStatus != "targetnotsurveyed",]
    n.R <- nrow(dat.R)
    n.over <- nrow(dat[(dat$panel == "OverSamp") & (dat[, evalstatus] == "targetsurveyed"), ])

# Calculate adjusted weights
    adj <- ifelse(popn == "Target", n.all * n.T/(n.eval * n.R), n.all/n.eval)
    adj.wt <- (dat[, wgt]) * (dat[, evalstatus] == "targetsurveyed") * adj 
    return(adj.wt)
}
