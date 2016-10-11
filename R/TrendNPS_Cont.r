#' @export
#' 
#' @title Trend for Continuous Outcomes from Complex Survey Designs
#'   
#' @description \code{TrendNPS_Cont} fits trend models of continuous outcomes 
#'   for four approaches:
#'   
#'   \enumerate{ \item the unreplicated linear mixed model as specified by 
#'   Piepho and Ogutu (2002) that does not incorporate design weights ("PO");
#'   
#'   \item simple linear regression of annual design-based estimates ("SLRDB");
#'   
#'   \item weighted linear regression of annual design-based estimates (WLRDB);
#'   
#'   \item and probability-weighted iterative generalized least squares of a 
#'   linear mixed model (PWIGLS) with 6 types of scaling.
#'   
#'   } Fixed effects structure includes a year term for trend estimation and an 
#'   optional two-level stratification factor. The function assumes random site 
#'   and year intercepts, and an optional random site-level slope effect may be 
#'   included.
#'   
#'   
#' @param dat	Data frame containing columns at least for \code{Site}, 
#'   \code{WYear}, \code{Year}, and the continuous outcome of interest \code{Y}.
#'   See section "Data frame \code{dat}" below.
#'   
#' @param method	Method for trend estimation entered as a string. Valid values 
#'   include \code{"PO"} for the Piepho and Ogutu (2002) unreplicated linear 
#'   mixed model, \code{"SLRDB"} for simple linear regression of design-based 
#'   estimates, \code{"WLRDB"} for weighted linear regression of design-based 
#'   estimates, or \code{"PWIGLS"} for probability-weighted iterative 
#'   generalized least squares (Pfeffermann et al. 1998, Asparouhov 2002).
#'   
#' @param slope	 Logical value indicating inclusion of a random site-level slope
#'   effect in the variance components structure used for the \code{"PO"} and 
#'   \code{"PWIGLS"} trend methods.  Site- and year-level random intercept terms
#'   included as default.
#'   
#' @param type	Scaling type when \code{method="PWIGLS"}. Valid values include 
#'   \code{"Aonly"}, \code{"A"}, \code{"AI"}, \code{"B"}, \code{"BI"} 
#'   \code{"C"}. See section "Options for variable \code{type}" below.
#'   
#' @param stratum	 Text string identifying an optional two-level stratification 
#'   factor in \code{dat}.
#'   
#' @param Y	 Text string indicating the outcome variable in \code{dat}.
#'   
#' @param lat	 Site latitudes using an equal-area projection, e.g., UTM or 
#'   Albers.
#'   
#' @param long	Site longitudes using an equal-area projection, e.g., UTM or 
#'   Albers.
#'   
#' @param stage1wt	Design weights from the original sample draw without 
#'   accounting for temporal revisit designs.
#'   
#' @param stage2wt	Panel inclusion weights for each site each year.
#'   
#' @param str1prop	Proportion of the first stratum in the population in a 
#'   two-level stratification variable.
#'   
#' @param nbhd	Logical indictor for using the neighborhood variance estimator 
#'   when the \code{method="WLRDB"}.  If \code{FALSE}, \code{"WLRDB"} calculates
#'   the design-based mean estimate's standard error assuming independent random
#'   sampling.
#'   
#' @return A displayed vector containing the trend estimate, along with its standard
#'   error, estimated variance components, and degrees of freedom used to test
#'   for trend and calculate confidence intervals.
#'   
#'   Additionally, a data frame containing each of 
#'   
#'   \tabular{rl}{
#'   
#'   \code{mu}      \tab	Estimated intercept from the trend model/ \cr
#'   \code{trend}   \tab	Estimated trend of the logged outcome \code{Y} from the trend model.  \cr
#'   \code{SEtrend} \tab	Estimated standard error of the trend estimate for the logged outcome \code{Y}. \cr
#'   \code{sig2a}   \tab	Estimated site-to-site variation. \cr
#'   \code{sig2b}   \tab	Estimated year-to-year variation. \cr
#'   \code{sig2t}   \tab	Estimated variation among site-level slopes. \cr
#'   \code{sig2e}   \tab	Estimated residual variation. \cr
#'   \code{eta}     \tab	Degrees of freedom used to test for trend and calculate confidence intervals. \cr
#'   
#'   }
#'   
#' @details For long-term monitoring, \code{TrendNPS_Cont} assumes sampling 
#'   occasions occur no more than once per unique record of \code{WYear} and 
#'   \code{Site} in data frame \code{dat}.  [Tell people what to do to get 
#'   around this?]
#'   
#'   Sampling (\code{stage1wt}) and temporal-revisit (\code{stage2wt}) design 
#'   weights must be adjusted for any missing data or frame error prior to trend
#'   analysis. See Oh and Scheuren (1983), Little and Rubin (2002), and 
#'   Starcevich et al. (2016).
#'   
#'   For background on the unweighted linear mixed model, see Piepho & Ogutu
#'   (2002).  For information on the probability-weighted generalized least
#'   squares approach, see Pfeffermann et al. (1998), Asparouhov (2006), and
#'   Starcevich et al. (2017).
#'   
#' @section Method requirements:
#'   
#'   Spatially balanced samples utilizing \code{method="SLRDB"} or 
#'   \code{method="WLRDB"} require specification of \code{lat} and \code{long} 
#'   on setting the neighborhood variance estimator \code{nbhd=TRUE}.  Set 
#'   option \code{nbhd=FALSE} for non-spatially balanced designs so that the 
#'   \code{"WLRDB"} method [and \code{"SLRDB"}?] assumes independent random 
#'   sampling for standard errors.  These methods also require specification of 
#'   \code{stratum}.
#'   
#'   The stratification variable \code{stratum} is required for 
#'   \code{method="PO"} and \code{method="PWIGLS"}.  Note that 
#'   \code{method="PWIGLS"} further requires a value for \code{type}.  See 
#'   section "Options for variable \code{type}" below.
#'   
#'   Variables \code{state1wt} and \code{stage2wt} must be specified when 
#'   \code{method="PWIGLS"}.
#'   
#' @section Data frame \code{dat}:
#'   
#'   Data frame \code{dat} requires at least variables \code{Site}, 
#'   \code{WYear}, \code{Year}, and \code{Y}.  Variable \code{Site} serves as a 
#'   statistical-unit text descriptior.
#'   
#'   Underlying statistical functions called by \code{TrendNPS_Cont} require two
#'   separate temporal variables: \code{WYear} is a scalar-year value, while
#'   \code{Year} is a factor.
#'   
#'   Ultimately, \code{WYear} enters the model as a trend effect, and designates
#'   the survey occasions used in the fixed-effects portion of the trend model. 
#'   The regression coefficient associated with \code{WYear} provides the basis 
#'   for trend estimation and testing. In line with Piepho & Ogutu (2002), 
#'   function \code{TrendNPS_Cont} reassigns the recorded temporal unit 
#'   \code{WYear} in which the outcome \code{Y} is least variable as zero, with 
#'   all other values \code{WYear} calculated with respect to this redefined 
#'   zero point.
#'   
#'   Variable \code{Year} enters the model as a random effect, and identifies 
#'   unique values for each survey occasion.  As a random effect, it provides 
#'   the basis for estimation of random year-to-year variation among survey 
#'   occasions.
#'   
#'   Variable \code{Y} must contain only positive non-zero values.  Analyses 
#'   utilizing \code{method="PO"} or \code{method="PWIGLS"} assumes variable 
#'   \code{Y} has been previously transformed to a logarithmic scale.
#'   
#' @section Options for variable \code{type}:
#'   
#'   Selection of \code{method="PWIGLS"} requires further specification of 
#'   argument \code{type}.  Valid options include
#'   
#'   \tabular{rl}{
#'   
#'   \code{"Aonly"} \tab Probability weighting but no scaling at either stage. \cr
#'   \code{"A"}     \tab Scaling for panel weights with the mean site-level design weight.  \cr
#'   \code{"AI"}    \tab Scaling for panel weights with the mean site-level design weight, but no site-level scaling.  \cr
#'   \code{"B"}     \tab Scaling for panel weights with the effective mean site-level design weight.  \cr
#'   \code{"BI"}    \tab Scaling for panel weights with the effective mean site-level design weight, but no site-level scaling.  \cr
#'   \code{"C"}     \tab Scaling only at the year level with the inverse of the average year-level weight. \cr
#'   
#'   }
#'   
#' @author Leigh Ann Starcevich and Jason Mitchell of Western EcoSystems 
#'   Technology, Inc.
#'   
#' @references Asparouhov, T. (2006).  General multi-level modeling with 
#'   sampling weights. Communications in Statistics - Theory and Methods 35: 
#'   439-460.
#'   
#'   Little, J.A.R., and D.B. Rubin. 2002. Statistical analysis with missing 
#'   data, 2nd edition. John Wiley and Sons, Inc., New Jersey.
#'   
#'   Oh, H.L., and F.J. Scheuren. 1983. Weighting adjustment for unit 
#'   nonresponse. Pages 143-184 in W.G. Madow, I. Olkin, and D.B. Rubin, 
#'   editors. Incomplete data and in sample surveys. (Vol. 2). Academic Press, 
#'   New York.
#'   
#'   Pfeffermann, D., C.J. Skinner, D.J. Holmes, H. Goldstein, and J. Rasbash 
#'   (1998).  Weighting for unequal selection probabilities in multilevel 
#'   models. Journal of the Royal Statistical Society, Series B 60(1): 23-40.
#'   
#'   Piepho, H.P., & Ogutu, J.O. (2002).  A simple mixed model for trend 
#'   analysis in wildlife populations.  Journal of Agricultural, Biological, and
#'   Environmental Statistics, 7(3), 350-360.
#'   
#'   L.A. Starcevich, G. DiDonato, T. McDonald, and J. Mitchell. 2016. A GRTS 
#'   user's manual for the SDrawNPS package: A graphical user interface for 
#'   generalized random tessellation stratified (GRTS) sampling and estimation. 
#'   Natural Resource Report NPS/PWRO/NRR-2016/1233. National Park Service, Fort
#'   Collins, Colorado.
#'   
#'   L.A. Starcevich, T. McDonald, A. Chung-MacCoubrey, A. Heard, and J. Nesmith
#'   (2017). Trend Estimation for Complex Survey Designs. Natural Resource 
#'   Report NPS/xxxx/NRR-2017/xxxx. National Park Service, Fort Collins, 
#'   Colorado.
#'   
#' @seealso \code{lme4}, \code{lmerTest}, \code{spsurvey}
#'   
#' @examples 
#' 
#' #  ---- Read example data set.
#' data(SEKIANC)
#' 
#' #  ---- Load dependent packages. 
#' pkgList <- c("lme4","lmerTest","spsurvey")
#' inst <- pkgList %in% installed.packages()
#' if (length(pkgList[!inst]) > 0) install.packages(pkgList[!inst])
#' lapply(pkgList, library, character.only = TRUE)
#' 
#' ###########################
#' #  ---- A1. PO method: full random effects model.
#' PO_ests <- TrendNPS_Cont(dat=SEKIANC, method="PO", type=NA, slope=TRUE, 
#'                          stratum=NA, Y="Y", lat="Lat", long="Long",
#'                          stage1wt="wgt", stage2wt="Panelwgt", str1prop=NA)
#' # Error: number of observations (=80) <= number of random effects (=92) for 
#' # term (1 + WYear | Site); the random-effects parameters and the residual 
#' # variance (or scale parameter) are probably unidentifiable
#' 
#' #  ---- Conclusion: the random effects are not estimable due to too few
#' #  ---- observations.  Reduce the random effects model by removing the 
#' #  ---- random slope term; i.e., set slope=FALSE.
#' 
#' #  ---- A2. PO method: reduced random effects model.
#' PO_ests <- TrendNPS_Cont(dat=SEKIANC, method="PO", type=NA, slope=FALSE,
#'                          stratum=NA, Y="Y", lat="Lat", long="Long",
#'                          stage1wt="wgt", stage2wt="Panelwgt", str1prop=NA)
#' round(PO_ests,4)
#' 
#' #     mu  trend SEtrend  sig2a sig2t sigat sig2b  sig2e    eta
#' #  3.964 0.0194  0.0322 0.3422     0     0     0 0.0161 3.6202
#' 
#' #  ---- Construct a 95%-confidence interval on the trend estimate of 
#' #  ---- the logged outcome.
#' CI <- c(PO_ests$trend - qt(1-(0.05/2), PO_ests$eta) * PO_ests$SEtrend, 
#'         PO_ests$trend + qt(1-(0.05/2), PO_ests$eta) * PO_ests$SEtrend)
#' round(c(PO_ests$trend, CI),4)
#' 
#' # 0.0194 -0.0737  0.1125
#' 
#' #  ---- Construct a 95%-confidence interval on the trend estimate of 
#' #  ---- the original outcome.  
#' round(exp(c(PO_ests$trend, CI)),4)
#' 
#' # 1.0196 0.9290 1.1191
#' 
#' #  ---- Conduct a two-sided test of trend for an alternative hypothesis 
#' #  ---- of trend in either direction.
#' t_stat <- PO_ests$trend/ PO_ests$SEtrend
#' p_value <- 2*(1-pt(abs(t_stat), PO_ests$eta, lower=TRUE))
#' round(data.frame(t_stat,p_value) ,4)
#' 
#' # 0.6038 0.5817
#' 
#' 
#' ###########################
#' #  ---- B1. WLRDB method.
#' WLRDB_ests <- TrendNPS_Cont(dat=SEKIANC, method="WLRDB", type=NA, 
#'                             slope=FALSE, stratum=NA, Y="Y", lat="Lat",
#'                             long="Long",stage1wt="wgt", stage2wt="Panelwgt",
#'                             str1prop=NA)
#' round(WLRDB_ests,4)
#' 
#' #      mu  trend SEtrend  sig2a sig2t sigat sig2b  sig2e   eta
#' #  4.1914 -0.0177  0.0492     0     0     0     0 0.0002   4
#' 
#' 
#' ###########################
#' #  ---- C1. PWIGLS with type="Aonly".
#' PWIGLS_Aonly_ests <- TrendNPS_Cont(dat=SEKIANC, method="PWIGLS", type="Aonly",
#'                                    slope=FALSE, stratum=NA, Y="Y", lat="Lat",
#'                                    long="Long",stage1wt="wgt", 
#'                                    stage2wt="Panelwgt", str1prop=NA)
#' round(PWIGLS_Aonly_ests,4)
#' 
#' #      mu  trend SEtrend  sig2a sig2t sigat  sig2b  sig2e    eta
#' #  3.8812 0.0383   0.031 6.2650     0     0 0.0732 0.0167 3.6202
#' 
TrendNPS_Cont<-function(dat,method,slope=TRUE,type=NA,stratum=NA,Y,lat=NA,long=NA,stage1wt=NA,stage2wt=NA,str1prop=NA,nbhd=TRUE) {

#   dat <- SEKIANC
#   method <- "PO"
#   slope <- FALSE
#   type <- NA
#   stratum <- NA
#   Y <- "Y"
#   lat <- "Lat"
#   long <- "Long"
#   stage1wt <- "wgt"
#   stage2wt <- "Panelwgt"
#   str1prop <- NA
#   nbhd <- TRUE
  
  # Calculate sample sizes
  Sites = sort(unique(dat$Site))
  ma = length(Sites)
  Years = as.character(sort(unique(dat$Year)))
  WYears = sort(unique(dat$WYear))
  mb = length(WYears)

  # Assign stratum
  if(!is.na(stratum)) dat$Stratum <- dat[,stratum]

  ############################################################################################
  # METHOD PO: Piepho & Ogutu (2002) linear mixed model 
  # Unreplicated model - one observation per Year and Site
  ############################################################################################
  if(method %in% c("PO","PWIGLS")) {	
  
    # create a logged outcome of interest
    dat$LogY = log(dat[,Y])
	  if(is.na(stratum)) {
	  	if(slope) fit<-lmer(LogY ~ WYear + (1|Year) +(1+WYear|Site), data=dat)
	  	if(!slope) fit<-lmer(LogY ~ WYear + (1|Year) +(1|Site), data=dat)
	  } 
    
	  if(!is.na(stratum)) {
	  	if(slope) fit<-lmer(LogY ~ WYear + Stratum+ (1|Year) +(1+WYear|Site), data=dat)
	  	if(!slope) fit<-lmer(LogY ~ WYear + Stratum+ (1|Year) +(1|Site), data=dat)
	  }

    if(method=="PWIGLS") {
    	fit_PO = fit
    	rm(fit)
    }

    if(method=="PO") {
	    mu<- summary(fit)$coef[1,1]
	    trend<- summary(fit)$coef[2,1]
	    SEtrend<- summary(fit)$coef[2,2]
	    if(slope) sig2a<- VarCorr(fit)$Site[1,1]		# var(ai)
	    if(!slope) sig2a<- VarCorr(fit)$Site			# var(ai)
	    if(slope) {
	 	    sig2t<- VarCorr(fit)$Site[2,2]		# var(ti)
	 	    sigat<- VarCorr(fit)$Site[1,2]		# cov(ai,ti)
	    }
	    if(!slope) sig2t<- sigat<- 0		
 	    sig2e<- attr(VarCorr(fit),"sc")^2
	    eta<- summary(fit)$coef[2,3]
	   #q.t<-qt(1-(alfa/2),eta)
	  }
  }		# End PO   ... jason & PWIGLS i think
  
  ###############################################################################	
  # Methods SLRDB and WLRDB: Simple linear regression on annual design-based estimates 
  # SLRDB: unweighted 
  # WLRDB: weighted by inverse of variance of each annual estimate
  ###############################################################################

  if(method %in% c("SLRDB","WLRDB")) {
    MeanEsts=data.frame(matrix(NA,mb,3))
    MeanEsts[,1] <- Years

    # Calculate annual design-based status estimates of the mean
	  for (g in 1:mb) {
      dat.g<-dat[dat$WYear==WYears[g],]
       
      if(is.na(stratum)) {
	      ests <-cont.analysis(
		      sites=data.frame(siteID=dat.g$Site, rep(TRUE,nrow(dat.g))), 
		      subpop= data.frame(siteID=dat.g$Site, Popn1=rep(1, nrow(dat.g))), 
		      design= data.frame(siteID=dat.g$Site, 
		        wgt= dat.g$wgt, xcoord = dat.g[,long], ycoord = dat.g[,lat]), 
		      data.cont= data.frame(siteID=dat.g$Site, Y=dat.g[,Y]), conf=90, vartype=ifelse(nbhd,"Local","SRS"))
	      MeanEsts[g,2:3] <-ests$Pct[ests$Pct$Statistic=="Mean",6:7]
      }
      
      if(!is.na(stratum)) {
	      ests <-cont.analysis(
		      sites=data.frame(siteID=dat.g$Site, rep(TRUE,nrow(dat.g))), 
		      subpop= data.frame(siteID=dat.g$Site, Popn1=rep(1, nrow(dat.g))), 
		      design= data.frame(siteID=dat.g$Site, 
		        wgt= dat.g$wgt, xcoord = dat.g[,long], ycoord = dat.g[,lat], stratum=dat.g[,stratum]), 
		      data.cont= data.frame(siteID=dat.g$Site, Y=dat.g[,Y]), conf=90, vartype=ifelse(nbhd,"Local","SRS"))
	      MeanEsts[g,2:3] <-ests$Pct[ests$Pct$Statistic=="Mean",6:7]
		  }
	  }
    
	  MeanEsts[,3] <- as.numeric(MeanEsts[,3])
	  names(MeanEsts) <- c("Year","Est.Mean","SE")
	  MeanEsts$WYear = WYears
	  
    # calculate trend in logged mean over time
	  if(method=="SLRDB")  fit<-lm(log(Est.Mean) ~ WYear, data=MeanEsts)
	  if(method=="WLRDB")  fit<-lm(log(Est.Mean) ~ WYear, weights=1/(SE^2), data=MeanEsts)

  	mu <- coef(fit)[1]				# mu
	  trend<- coef(fit)[2]				# slope
	  SEtrend<-sqrt(vcov(fit)[2,2])			# SE
	  sig2a<- sig2t<- sigat<- sig2b<- 0
	  sig2e<- (summary(fit)$sigma)^2
	  eta<-summary(fit)$df[2]
	 #q.t <- qt(1-(alfa/2),eta)
  }  # End SLRDB/WLRDB

  ############################################################################################
  # PWIGLS methods
  ############################################################################################

  if(method=="PWIGLS") {
    if(!type %in% c("Aonly","A","AI","B","BI","C")) return("PWIGLS scaling type not recognized.")

	  fit<-PWIGLS_ALL(Z=getME(fit_PO,"Z"),dat=dat,stage1wt=stage1wt,stage2wt=stage2wt,type=type,stratum=stratum,slope=slope)
  #	fit<-fit[[1]]
	  mu<- fixef(fit)[1]				# intercept
    if(slope) {
	    sig2a<- VarCorr(fit)$Site[1,1]		# var(ai)
	    sig2t<- VarCorr(fit)$Site[2,2]		# var(ti)
	    sigat<- VarCorr(fit)$Site[1,2]		# cov(ai,ti)
    	sig2b<- VarCorr(fit)$Year[1]			# var(bj)
	    sig2e<- attr(VarCorr(fit),"sc")^2		# var(eij)
	  } 
    if(!slope) {
	    sig2a<- VarCorr(fit)$Site[1]			# var(ai)
	    sig2t<- sigat<- 0
	    sig2b<- VarCorr(fit)$Year[1]			# var(bj)
	    sig2e<- attr(VarCorr(fit),"sc")^2		# var(eij)
	  }
	  varij = sig2a + sig2b + dat$WYear*sigat +(dat$WYear^2)*sig2t +sig2e

    if(is.na(stratum)) {		# No strata
	    trend<- fixef(fit)[2]				# slope
	    eta<- summary(fit_PO)$coef[2,3]		# Use PO df
    # Linearization variance -- Pfeffermann 1988
    	SEtrend<-sqrt(LinearizationVar(Site=dat$Site,
		    wij=dat[,stage1wt]*dat[,stage2wt],
		    xij=dat$WYear,
		    eij=dat$LogY-mu-(trend*dat$WYear),
		    varYij=varij)) 
	  }	# END No strata

    if(!is.na(stratum))  {		# Stratification with two levels
	    trend1<- fixef(fit)[2]				# slope of stratum 1
	    trend2<- fixef(fit)[4]				# slope of stratum 2 - slope of stratum 1
	    trend=trend1 + (1-str1prop)*trend2		# population estimate of trend
	    eta<- sum(summary(fit_PO)$coef[c(2,4),3])	# Use PO df - sum of df for two trend coefs
	    Ind1 = as.numeric(dat$Stratum==levels(dat$Stratum)[1])
	    SEtrend<-sqrt(LinearizationVar_StRS(Site=dat$Site,
		    wij=dat$wgt*dat$Stage2wt,
		    xij=data.frame(dat$WYear,Ind1),
		    eij=dat$LogY-mu-(trend1*dat$WYear)-(fixef(fit)[3]*Ind1)-(trend2*dat$WYear*Ind1),
		    varYij=varij,
		    str1prop=str1prop))
	  }   # END Stratification with two levels
  }			# End PWIGLS

  if(slope==TRUE){
    ans=data.frame(mu,trend,SEtrend,sig2a,sig2t,sigat,sig2b,sig2e,eta)	# ,q.t
    names(ans)=c('mu','trend','SEtrend','sig2a','sig2t','sigat','sig2b','sig2e','eta')
  } else {
    ans=data.frame(mu,trend,SEtrend,sig2a,sig2t,sigat,sig2e,eta)	
    names(ans)=c('mu','trend','SEtrend','sig2a','sig2t','sigat','sig2e','eta')
  }
  rownames(ans) = NULL
  return(ans)
}
