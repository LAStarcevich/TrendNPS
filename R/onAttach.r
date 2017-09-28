#' @import  spsurvey stats
#' @importFrom lmerTest lmer
#' @importFrom lme4 VarCorr getME fixef glmer glmerControl
.onAttach <- function(libname, pkgname){
  
  v <- utils::packageVersion(pkgname) 
  
  packageStartupMessage( paste(pkgname, " (vers ", v ,")", sep=""))  
  packageStartupMessage( paste("Developed by WEST, Inc. for use by the National Park Service") )
  
  #   ---- Install all the packages we could ever possible need.
  list.of.packages <- c("lmerTest", "spsurvey", "stats", "lme4") 
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages) 
  lapply(new.packages, require, character.only=TRUE)
  
}
