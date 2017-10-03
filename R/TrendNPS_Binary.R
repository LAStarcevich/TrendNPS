#' @export
#' 
#' @title Trend for Binary Outcomes from Complex Survey Designs
#'   
#' @description \code{TrendNPS_Binary} fits trend models of binary outcomes for
#'   four approaches:
#'   
#'   \enumerate{ \item the unreplicated generalized linear mixed model with 
#'   variance components as specified by Piepho and Ogutu (2002) that does 
#'   not incorporate design weights ("PO");
#'   
#'   \item simple linear regression of annual design-based estimates ("SLRDB");
#'   
#'   \item weighted linear regression of annual design-based estimates ("WLRDB");
#'   
#'   \item and probability-weighted iterative generalized least squares ("PWIGLS") 
#    for 6 types of scaling.
#'   
#'   } Fixed effects structure includes a year term for trend estimation and an 
#'   optional two-level stratification factor. The function assumes random site 
#'   and year intercepts, and an optional random site-level slope effect may be 
#'   included.
#'   
#' @param alpha The probability of a Type I error. 
#'   
#' @param dat	Data frame containing columns at least for \code{Site}, 
#'   \code{WYear}, \code{Year}, and the binary outcome of interest \code{Y}.
#'   See section "Data frame \code{dat}" below.
#'   
#' @param method	Method for trend estimation entered as a string. Valid values 
#'   include \code{"PO"} for the generalized linear mixed model extension of the 
#'   Piepho and Ogutu (2002) unreplicated linear mixed model, \code{"SLRDB"} for 
#'   simple linear regression of design-based estimates, \code{"WLRDB"} for 
#'   weighted linear regression of design-based estimates, or \code{"PWIGLS"} for 
#'   probability-weighted iterative generalized least squares (Pfeffermann et al. 
#'   1998, Asparouhov 2002).
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
#' @param Y	 Text string indicating the binary outcome variable in \code{dat}.
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
#' @return A list containing three or four elements, depending on the specified trend method. 
#'   The first element of the list, ModelEstimates, contains a data frame with the trend model output:
#'   
#'   \tabular{rl}{
#'   
#'   \code{mu}      \tab	Estimated intercept from the trend model. \cr
#'   \code{trend}   \tab	Estimated trend of the outcome \code{Y} on the link function scale from the trend model.  \cr
#'   \code{SEtrend} \tab	Estimated standard error of the trend estimate. \cr
#'   \code{sig2a}   \tab	Estimated site-to-site variation. \cr
#'   \code{sig2b}   \tab	Estimated year-to-year variation. \cr
#'   \code{sig2t}   \tab	Estimated variation among site-level slopes. \cr
#'   \code{sig2e}   \tab	Estimated residual variation. \cr
#'   \code{eta}     \tab	Degrees of freedom used to test for trend and calculate confidence intervals. \cr
#'   
#'   }
#'   
#'   The contents of the second list element, TrendTest, include the trend coefficient on the link scale, 
#'   the trend standard error, the t-statistic, degrees of freedom, and p-value for a two-sided test of no trend.
#'   
#'   The contents of the third list element, TrendCI, include the back-transformed trend on the original scale 
#'   and the lower and upper (1-\code{alpha})*100\%-confidence interval bounds. 
#'   
#'   A fourth list element, DBests, contains the annual design-based estimates and is returned only when the 
#'   SLRDB and WLRDB trend approaches are used.  The design-based estimates may be used to assess the assumption
#'   of independence for the linear regression models.  The output of the DBests list element include:
#'   
#'   \tabular{rl}{
#'   
#'   \code{Year}      \tab	Year of survey. \cr
#'   \code{Est.Mean}   \tab	Design-based estimate of the mean of the outcome of interest. \cr
#'   \code{SE} \tab	   Estimated standard error of design-based mean estimate via neighborhood 
#'   variace estimator if \code{nbhd = TRUE}. \cr
#'   \code{WYear}   \tab	WYear variable. \cr
#'   \code{Resid}     \tab	Residuals from the linear model fit. \cr
#'   
#'   }
#'   
#' @details For long-term monitoring, \code{TrendNPS_Binary} assumes a single annual sampling 
#'   occasion per unique record of \code{WYear} and \code{Site} in data frame \code{dat}.  
#'   
#'   Sampling (\code{stage1wt}) and temporal-revisit (\code{stage2wt}) design 
#'   weights must be adjusted for any missing data or frame error prior to trend
#'   analysis. See Oh and Scheuren (1983), Little and Rubin (2002), and 
#'   Starcevich et al. (2016).
#'   
#'   For background on the variance structure, see Piepho & Ogutu (2002).  For information 
#'   on the probability-weighted generalized least squares approach, see Pfeffermann 
#'   et al. (1998), Asparouhov (2006), and Starcevich et al. (2017).
#'   
#' @section Method requirements:
#'   
#'   Spatially balanced samples utilizing \code{method="SLRDB"} or 
#'   \code{method="WLRDB"} require specification of \code{lat} and \code{long} 
#'   on setting the neighborhood variance estimator \code{nbhd=TRUE}.  If 
#'   spatial balance cannot be assumed, set option \code{nbhd=FALSE} for 
#'   non-spatially balanced designs so that the \code{"SLRDB"} and 
#'   \code{"WLRDB"} methods assume independent random sampling for standard 
#'   errors.
#'   
#'   The stratification variable \code{stratum} is required for all methods when
#'   design strata are used.  Note that \code{method="PWIGLS"} further requires
#'   a value for \code{type}.  See section "Options for variable \code{type}"
#'   below.
#'   
#'   Variables \code{stage1wt} and \code{stage2wt} must be specified for all 
#'   methods except the unweighted \code{method="PO"}.
#'   
#' @section Data frame \code{dat}:
#'   
#'   Data frame \code{dat} requires at least variables \code{Site}, 
#'   \code{WYear}, \code{Year}, and \code{Y}.  Variable \code{Site} defines the 
#'   sampling unit.
#'   
#'   Underlying statistical functions called by \code{TrendNPS_Cont} require two
#'   separate temporal variables: \code{Year} is a factor and \code{WYear} is a 
#'   scalar-year value.  \code{WYear} is shifted so that the year for which \code{Y}
#'   is least variable is centered on 0. See Piepho & Ogutu (2002) for more 
#'   information.
#'   
#'   Ultimately, \code{WYear} represents time in the fixed-effects portion of the  
#'   trend model. The regression coefficient associated with \code{WYear} provides
#'   the basis for trend estimation and testing. 
#'   
#'   Variable \code{Year} enters the model as a random effect and provides 
#'   the basis for estimation of random year-to-year variation among survey 
#'   occasions.
#'   
#'   Variable \code{Y} must contain only binary values coded as 0 or 1.  
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
#'   Little, J. A. R., and D. B. Rubin. (2002). Statistical analysis with missing 
#'   data, 2nd edition. John Wiley and Sons, Inc., New Jersey.
#'   
#'   Oh, H. L., and F. J. Scheuren. (1983). Weighting adjustment for unit 
#'   nonresponse. Pages 143-184 in W.G. Madow, I. Olkin, and D.B. Rubin, 
#'   editors. Incomplete data and in sample surveys. (Vol. 2). Academic Press, 
#'   New York.
#'   
#'   Pfeffermann, D., C. J. Skinner, D. J. Holmes, H. Goldstein, and J. Rasbash 
#'   (1998).  Weighting for unequal selection probabilities in multilevel 
#'   models. Journal of the Royal Statistical Society, Series B 60(1): 23-40.
#'   
#'   Piepho, H. P., & Ogutu, J. O. (2002).  A simple mixed model for trend 
#'   analysis in wildlife populations.  Journal of Agricultural, Biological, and
#'   Environmental Statistics, 7(3), 350-360.
#'   
#'   Starcevich, L. A., G. DiDonato, T. McDonald, and J. Mitchell. (2016). A GRTS 
#'   user's manual for the SDrawNPS package: A graphical user interface for 
#'   generalized random tessellation stratified (GRTS) sampling and estimation. 
#'   Natural Resource Report NPS/PWRO/NRR-2016/1233. National Park Service, Fort
#'   Collins, Colorado.
#'   
#'   Starcevich, L. A., T. McDonald, A. Chung-MacCoubrey, A. Heard, and J. Nesmith.
#'   (2017). Trend Estimation for Complex Survey Designs. Natural Resource 
#'   Report NPS/xxxx/NRR-2017/xxxx. National Park Service, Fort Collins, 
#'   Colorado.

#'   R. Wolfinger. (1993). Laplace’s approximation for nonlinear mixed models. 
#'   Biometrika 80(4): 791-795. 
#'   
#' @seealso \code{lme4}, \code{lmerTest}, \code{spsurvey}
#'   
#' @examples 
#' 
#' #  ---- Read example data set.
#' data(Cover_Acro)
#' 
#' #  ---- Load dependent packages. 
#' pkgList <- c("lme4","lmerTest","spsurvey")
#' inst <- pkgList %in% installed.packages()
#' if (length(pkgList[!inst]) > 0) install.packages(pkgList[!inst])
#' lapply(pkgList, library, character.only = TRUE)
#' 
#' ###########################
#' # Example 1: Trend analysis of a binary indicator of Acrosophonia presence
#' # with the PO approach: stratification and full random effects model.
#' TrendAcro_PO_StRS = TrendNPS_Binary(alpha=0.1,
#' dat=Cover_Acro,method="PO",slope=TRUE,type=NA,stratum="Park",Y="Y",
#' stage1wt="wgt",stage2wt="PanelWt",str1prop=0.13227) 

#' # $ModelEstimates
#' #           mu    trend    SEtrend     sig2a      sig2t       sigat     sig2b
#' #   -0.2992415 0.221928 0.09341594 0.2867116 0.01394005 -0.04732509 0.2411157
#' 
#' # $TrendTest
#' #      trend    SEtrend   z-stat    pvalue
#' #   0.221928 0.09341594 2.375698 0.0175158
#' 
#' # $TrendCIofOddsRatio
#' #   Annual Pct Change     CI low   CI high
#' #           0.2484815 0.07065702 0.4558408
#'
#' 
#' # Example 2: Trend analysis of a binary indicator of Acrosophonia presence
#' # with the SLRDB approach for stratification.
#' TrendAcro_SLRDB_StRS = TrendNPS_Binary(alpha=0.1,
#' dat=Cover_Acro,method="SLRDB",slope=TRUE,type=NA,stratum="Park",Y="Y",
#' lat="Lat",long="Long",stage1wt="wgt",stage2wt="PanelWt",
#' str1prop=0.13227) 
#' TrendAcro_SLRDB_StRS
#' # $ModelEstimates
#' #           mu    trend    SEtrend sig2a sig2t sigat sig2b
#' #   -0.9909923 0.211127 0.08966504     0     0     0     0
#' # 
#' # $TrendTest
#' #      trend    SEtrend   z-stat     pvalue
#' #   0.211127 0.08966504 2.354619 0.01854169
#' # 
#' # $TrendCIofOddsRatio
#' #  Annual Pct Change     CI low   CI high
#' #           0.2350692 0.06570992 0.4313426
#' # 
#' # $DBests
#' #   Year  Est.Mean         SE WYear      Resid
#' #   2008 0.1604049 0.06934590     0 -0.6642261
#' #   2009 0.2652986 0.05963881     1 -0.2387430
#' #   2010 0.4121248 0.10439841     2  0.2135496
#' #   2011 0.6416667 0.07341420     3  0.9402165
#' #   2012 0.5439923 0.11698602     4  0.3229097
#' #   2013 0.5875758 0.03719284     5  0.2893100
#' #   2014 0.4954618 0.06389670     6 -0.2939231
#' #   2015 0.4794624 0.04936033     7 -0.5690937
#' 
#' num <- TrendAcro_SLRDB_StRS$DBests$Est.Mean
#' denom <- 1-TrendAcro_SLRDB_StRS$DBests$Est.Mean 
#' RrendAcro_SLRDB_StRS$DBests$LogOdds = num / denom
#' plot(TrendAcro_SLRDB_StRS$DBests$Year,
#'      TrendAcro_SLRDB_StRS$DBests$LogOdds,
#'      xlab="Year", ylab="Odds ratio of cover proportion")
#' lines(2008:2015, exp(-0.9909923 + 0.211127*(0:7)), col=2)
#'
#' # Plot trend on proportional cover scale
#' plot(TrendAcro_SLRDB_StRS$DBests$Year,
#' 		TrendAcro_SLRDB_StRS$DBests$Est.Mean, 
#' 		xlab="Year", ylab="Cover proportion")
#' lines(2008:2015, expit(-0.9909923 + 0.211127*(0:7)), col=2)


TrendNPS_Binary<-function(alpha,dat,method,slope=TRUE,type=NA,stratum=NA,Y,lat=NA,long=NA,stage1wt=NA,stage2wt=NA,str1prop=NA,nbhd=TRUE) {

# dat = data set for trend analysis
# dat must contain fields for Site, WYear (scalar year value), and Year (year factor)
# other covariate fields in dat may be used in the trend model, but WYear must be used to estimate trend

# Calculates a trend model specified by method:
# 'PO' = Piepho and Ogutu (2002) unreplicated linear mixed model 
# 'SLRDB' = Simple linear regression of design-based estimates
# 'WLRDB' = Weighted linear regression of design-based estimates using the inverse of the variance of the estimate as the weight
# 'PWIGLS' = probability-weighted iterative generalized least squares (Pfeffermann et al. 1998, Asparouhov 2002)
#    PWIGLS requires the input, type =
#    'Aonly' = probability weighting but no scaling at either stage
#    'A' = scales panel weights with the mean site-level design weight
#    'AI' = scales panel weights with the mean site-level design weight, but removes site-level scaling
#    'B' = scales panel weights with the effective mean site-level design weight
#    'BI' = scales panel weights with the effective mean site-level design weight, but removes site-level scaling
#    'C' = scales only at the year level with the inverse of the average year-level weight 

# Inputs needed for SLRDB and WLRDB methods
# stratum = column name of stratification variable in dat
#    stratum is used directly in the SLRDB and WLRDB approaches
#    The stratification variable must be included in the fe input for methods PO and PWIGLS
# Y = column name of the outcome of interest for design-based estimates
# lat = column name of the latitude 
# long = column name of the longitude

# Inputs needed for PWIGLS method
# stage1wt = column name of the weight from the sampling design, adjusted for frame and/or nonresponse error if needed
# stage2wt = column name of the panel weight from temporal revisit design

# srt1prop = proportion of the population represented by stratum 1 in a two-level stratification variable 
#   str1 will be the stratum listed first when stratum levels are ordered alphabetically.

# requires package lme4, lmerTest, spsurvey

# Calculate sample sizes
Sites = sort(unique(dat$Site))
ma = length(Sites)
Years = as.character(sort(unique(dat$Year)))
WYears = sort(unique(dat$WYear))
mb = length(WYears)

# Assign stratum
if(!is.na(stratum)) dat$Stratum <- as.factor(dat[,stratum])
dat$Y <- dat[,Y]

############################################################################################
# METHOD PO: Piepho & Ogutu (2002) linear mixed model 
# Unreplicated model - one observation per Year and Site
############################################################################################
if(method %in% c("PO","PWIGLS")) {	
	if(is.na(stratum)) {
		if(slope) fit<-glmer(Y ~ WYear + (1|Year) +(1+WYear|Site), family=binomial(link=logit), data=dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6))) 
		if(!slope) fit<-glmer(Y ~ WYear + (1|Year) +(1|Site), family=binomial(link=logit), data=dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6))) 
	}
	if(!is.na(stratum)) {
		if(slope) fit<-glmer(Y ~ WYear * Stratum + (1|Year) +(1+WYear|Site),family=binomial(link=logit), data=dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
		if(!slope) fit<-glmer(Y ~ WYear * Stratum + (1|Year) +(1|Site),family=binomial(link=logit), data=dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
	}

if(method=="PWIGLS") {
	fit_PO = fit
	rm(fit)
	}

if(method=="PO") {
	if(!slope) {
		sig2a<- VarCorr(fit)$Site		# var(ai)
		sig2t<- sigat<- 0
	}
	if(slope) {
		sig2a<- VarCorr(fit)$Site[1,1]		# var(ai)
		sig2t<- VarCorr(fit)$Site[2,2]		# var(ti)
		sigat<- VarCorr(fit)$Site[1,2]		# cov(ai,ti)
	}
	sig2b<- VarCorr(fit)$Year
	if(is.na(stratum))  {		# No stratification
		mu<- summary(fit)$coef[1,1]
		trend<- summary(fit)$coef[2,1]
		SEtrend<- summary(fit)$coef[2,2]
			}
	if(!is.na(stratum))  {		# Stratification with two levels
		mu<- summary(fit)$coef[1,1]
		c.vec = matrix(c(0,1,0,1-str1prop),1,4)
		trend = c.vec %*% fixef(fit)
		SEtrend<- sqrt(c.vec%*% as.matrix(vcov(fit))%*%t(c.vec))
	}

}		# End PO
} 		# End PO or PWIGLS
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
		# summarize proportion of positive trials (hits) at the site level
		if(is.na(stratum)) {
		dat.g = aggregate(list(MeanY=dat.g[,Y]), list(Site=dat.g$Site, long=dat.g[,long], lat=dat.g[,lat],wgt= dat.g[,stage1wt]), mean)
		ests <-cont.analysis(
		sites=data.frame(siteID=dat.g$Site, rep(TRUE,nrow(dat.g))), 
		subpop= data.frame(siteID=dat.g$Site, Popn1=rep(1, nrow(dat.g))), 
		design= data.frame(siteID=dat.g$Site, 
		wgt= dat.g$wgt, xcoord = dat.g$long, ycoord = dat.g$lat), 
		data.cont= data.frame(siteID=dat.g$Site, Y=dat.g$MeanY), conf=90, vartype=ifelse(nbhd,"Local","SRS"))
			MeanEsts[g,2:3] <-ests$Pct[ests$Pct$Statistic=="Mean",6:7]
		}
		if(!is.na(stratum)) {
		dat.g = aggregate(list(MeanY=dat.g[,Y]), list(Site=dat.g$Site, long=dat.g[,long], lat=dat.g[,lat],wgt= dat.g[,stage1wt],stratum=dat.g[,stratum]), mean)
		ests <-cont.analysis(
		sites=data.frame(siteID=dat.g$Site, rep(TRUE,nrow(dat.g))), 
		subpop= data.frame(siteID=dat.g$Site, Popn1=rep(1, nrow(dat.g))), 
		design= data.frame(siteID=dat.g$Site, 
		wgt= dat.g$wgt, xcoord = dat.g$long, ycoord = dat.g$lat, stratum=dat.g$stratum), 
		data.cont= data.frame(siteID=dat.g$Site, Y=dat.g$MeanY), conf=90, vartype=ifelse(nbhd,"Local","SRS"))
			MeanEsts[g,2:3] <-ests$Pct[ests$Pct$Statistic=="Mean",6:7]
		}
	}
	MeanEsts[,3] <- as.numeric(MeanEsts[,3])
	names(MeanEsts) <- c("Year","Est.Mean","SE")
	MeanEsts$WYear = WYears

# calculate trend in logged mean over time
	if(method=="SLRDB")  fit<-lm(logit(Est.Mean) ~ WYear, data=MeanEsts)
	if(method=="WLRDB")  fit<-lm(logit(Est.Mean) ~ WYear, weights=1/(MeanEsts$SE^2), data=MeanEsts)

  	mu <- coef(fit)[1]				# mu
	trend<- coef(fit)[2]				# slope
	SEtrend<-sqrt(vcov(fit)[2,2])			# SE
	sig2a<- sig2t<- sigat<- sig2b<- 0
	sig2e<- (summary(fit)$sigma)^2
}  # End SLRDB/WLRDB

############################################################################################
# PWIGLS methods
############################################################################################

if(method == "PWIGLS") {
if(!type %in% c("Aonly","A","AI","B","BI","C")) return("PWIGLS scaling type not recognized.")

	fit<-PWIGLS_Binary(Z=getME(fit_PO,"Z"),dat=dat,stage1wt=stage1wt,stage2wt=stage2wt,type=type,stratum=stratum,slope=slope,PO.fitted=predict(fit_PO,type="response"),PO.resid=resid(fit_PO,type="response"))
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
	

if(is.na(stratum)) {		# No strata
	trend<- fixef(fit)[2]				# slope
	eta<- summary(fit)$coef[2,3]		# Use PO df
# Linearization variance -- Pfeffermann 1988
	SEtrend <- SEbeta<- sqrt(vcov(fit_PO)[2,2])
	}	# END No strata


if(!is.na(stratum))  {		# Stratification with two levels
	trend1<- fixef(fit)[2]				# slope of stratum 1
	trend2<- fixef(fit)[4]				# slope of stratum 2 - slope of stratum 1
	trend=trend1 + (1-str1prop)*trend2		# population estimate of trend
	Cont = matrix(c(0,1,0,str1prop),1,4)
	SEtrend<-as.numeric(sqrt(Cont%*%vcov(fit)%*%t(Cont)))	}   # END Stratification with two levels
}			# End PWIGLS

ans.model = data.frame(mu,trend,SEtrend,sig2a,sig2t,sigat,sig2b)	
names(ans.model)=c('mu','trend','SEtrend','sig2a','sig2t','sigat','sig2b')
rownames(ans.model) = NULL

# z test
ans.trendtest = data.frame(trend = ans.model[2], SE = ans.model[3], z.stat=ans.model[2]/ans.model[3], pvalue=2*pnorm(as.numeric(abs(ans.model[2]/ans.model[3])),lower.tail =FALSE))
names(ans.trendtest)[3] = 'z-stat'

# CI on trend:
CI = data.frame(ans.model[2],ans.model[2]-qnorm(1-(alpha/2))*ans.model[3],ans.model[2]+qnorm(1-(alpha/2))*ans.model[3])
ans.trendCI = exp(CI)-1		# return CI of odds ratio
names(ans.trendCI) = c("Annual Pct Change", "CI low", "CI high")

if(!method %in% c("SLRDB","WLRDB")) return(list(ModelEstimates = ans.model, TrendTest = ans.trendtest, TrendCIofOddsRatio = ans.trendCI))
if(method %in% c("SLRDB","WLRDB")) return(list(ModelEstimates = ans.model, TrendTest = ans.trendtest, TrendCIofOddsRatio = ans.trendCI, DBests = data.frame(MeanEsts, Resid = resid(fit))))

}


