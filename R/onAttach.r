#' @import  spsurvey stats
#' @importFrom lmerTest lmer
#' @importFrom lme4 VarCorr getME fixef glmer glmerControl
.onAttach <- function(libname, pkgname){
  
  v <- utils::packageVersion(pkgname) 
  
  packageStartupMessage( paste(pkgname, " (vers ", v ,")", sep=""))  
  packageStartupMessage( paste("Developed by WEST, Inc. for use by the National Park Service") )
  
}
