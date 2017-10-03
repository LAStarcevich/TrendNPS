#' @export
#' 
#' @title Trend for Count Outcomes from Complex Survey Designs
#'   
#' @description \code{TrendNPS_Count} fits trend models of count outcomes 
#  for four approaches:
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
#    of predictions from a LaPlace approximation for 6 types of scaling.
#'   
#'   } Fixed effects structure includes a year term for trend estimation and an 
#'   optional two-level stratification factor. The function assumes random site 
#'   and year intercepts, and an optional random site-level slope effect may be 
#'   included.
#'   
#' @param alpha The probability of a Type I error. 
#'   
#' @param dat	Data frame containing columns at least for \code{Site}, 
#'   \code{WYear}, \code{Year}, and the count outcome of interest \code{Y}.
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
#' @details For long-term monitoring, \code{TrendNPS_Count} assumes a single annual sampling 
#'   occasion per unique record of \code{WYear} and \code{Site} in data frame \code{dat}.  
#'   
#'   Sampling (\code{stage1wt}) and temporal-revisit (\code{stage2wt}) design 
#'   weights must be adjusted for any missing data or frame error prior to trend
#'   analysis. See Oh and Scheuren (1983), Little and Rubin (2002), and 
#'   Starcevich et al. (2016).
#'   
#'   For background on the variance structure, see Piepho & Ogutu (2002).  For more information 
#'   on the LaPlace approximation, see Wolfinger (1993) and Raudenbush et al. (2000). For 
#'   information on the probability-weighted generalized least squares approach, see Pfeffermann 
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
#'   Variable \code{Y} must contain only integer count values coded as numeric values 0, 1, 2, 3, ...  
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
#'   Raudenbush, S.W., M-L Yang, and M. Yosef. (2000). Maximum likelihood for 
#'   generalized linear models with nested random effects via a high-order, 
#'   multivariate Laplace approximation. Journal of Computational and Graphical 
#'   Statistics 9(1): 141-157. 

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
#' data(Seastar)
#' 
#' #  ---- Load dependent packages. 
#' pkgList <- c("lme4","lmerTest","spsurvey")
#' inst <- pkgList %in% installed.packages()
#' if (length(pkgList[!inst]) > 0) install.packages(pkgList[!inst])
#' lapply(pkgList, library, character.only = TRUE)
#' 
#' ###########################
#' # Example 1: Trend analysis of annual site-level counts of leather sea stars
#' # with the PO approach: stratification and full random effects model.
#' TrendSeaster_PO_StRS = TrendNPS_Count(alpha=0.1,
#' dat=Seastar,method="PO",slope=TRUE,type=NA,stratum="Park",Y="Y",
#' stage1wt="wgt",stage2wt="PanelWt",str1prop=0.13227) 

#' # $ModelEstimates
#' #         mu      trend    SEtrend    sig2a       sig2t        sigat      sig2b
#' #   1.470868 0.04604476 0.04518698 0.248528 0.001259151 -0.007764177 0.03724626
#' # 
#' # $TrendTest
#' #        trend    SEtrend   z-stat    pvalue
#' #   0.04604476 0.04518698 1.018983 0.3082111
#' # 
#' # $TrendCI
#' #  Annual Pct Change      CI low   CI high
#' #         0.04712128 -0.02788504 0.1279149
#' 
#' 
#' # Example 2: Trend analysis of annual site-level counts of leather sea stars
#' # with the SLRDB approach for stratification.
#' TrendSeastar_SLRDB_StRS = TrendNPS_Count(alpha=0.1,
#' dat=Seastar,method="SLRDB",slope=TRUE,type=NA,stratum="Park",Y="Count",
#' lat="Lat",long="Long", stage1wt="wgt",stage2wt="PanelWt",str1prop=0.13227) 
#' 
#' # TrendSeastar_SLRDB_StRS
#' # $ModelEstimates
#' #        mu     trend    SEtrend sig2a sig2t sigat sig2b
#' #  2.687804 0.0317377 0.04920822     0     0     0     0
#' # 
#' # $TrendTest
#' #       trend    SEtrend    z-stat    pvalue
#' #   0.0317377 0.04920822 0.6449674 0.5189483
#' # 
#' # $TrendCI
#' #   Annual Pct Change      CI low   CI high
#' #          0.03224671 -0.04801179 0.1192715
#' # 
#' # $DBests
#' #   Year Est.Mean       SE WYear       Resid
#' #   2009 11.66248 2.268814     0 -0.23142756
#' #   2010 14.66567 4.057595     1 -0.03403269
#' #   2011 18.60000 6.206449     2  0.17188195
#' #   2012 17.17467 5.045078     3  0.06041852
#' #   2013 21.71013 5.658295     4  0.26302412
#' #   2014 20.29531 4.128595     5  0.16389719
#' #   2015 11.99475 2.917410     6 -0.39376154
#' 
#' # Plot trend on original scale
#' plot(TrendSeastar_SLRDB_StRS$DBests$Year,TrendSeastar_SLRDB_StRS$DBests$Est.Mean, xlab="Year", ylab="Mean Sea Star Counts")
#' lines(2009:2015, exp(2.687804 + 0.0317377*(0:6)), col=2)
#' 
#' 
#' 
TrendNPS_Count<-function(alpha,dat,method,slope=TRUE,type=NA,stratum=NA,Y,lat=NA,long=NA,stage1wt=NA,stage2wt=NA,str1prop=NA,nbhd=TRUE) {

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
		if(slope) fit<-glmer(Y ~ WYear + (1|Year) +(1+WYear|Site), family=poisson, data=dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6))) 
		if(!slope) fit<-glmer(Y ~ WYear + (1|Year) +(1|Site), family=poisson, data=dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6))) 
	}
	if(!is.na(stratum)) {
		if(slope) fit<-glmer(Y ~ WYear * Stratum + (1|Year) +(1+WYear|Site),family=poisson, data=dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
		if(!slope) fit<-glmer(Y ~ WYear * Stratum + (1|Year) +(1|Site),family=poisson, data=dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
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
		if(is.na(stratum)) {
			ests <-cont.analysis(
		sites=data.frame(siteID=dat.g$Site, rep(TRUE,nrow(dat.g))), 
		subpop= data.frame(siteID=dat.g$Site, Popn1=rep(1, nrow(dat.g))), 
		design= data.frame(siteID=dat.g$Site, 
		wgt= dat.g$wgt, xcoord = dat.g[,long], ycoord = dat.g[,lat]), 
		data.cont= data.frame(siteID=dat.g$Site, Y=dat.g$Y), conf=90, vartype=ifelse(nbhd,"Local","SRS"))
			MeanEsts[g,2:3] <-ests$Pct[ests$Pct$Statistic=="Mean",6:7]
		}
		if(!is.na(stratum)) {
			ests <-cont.analysis(
		sites=data.frame(siteID=dat.g$Site, rep(TRUE,nrow(dat.g))), 
		subpop= data.frame(siteID=dat.g$Site, Popn1=rep(1, nrow(dat.g))), 
		design= data.frame(siteID=dat.g$Site, 
		wgt= dat.g$wgt, xcoord = dat.g[,long], ycoord = dat.g[,lat], stratum=dat.g[,stratum]), 
		data.cont= data.frame(siteID=dat.g$Site, Y=dat.g$Y), conf=90, vartype=ifelse(nbhd,"Local","SRS"))
			MeanEsts[g,2:3] <-ests$Pct[ests$Pct$Statistic=="Mean",6:7]
		}
	}
	MeanEsts[,3] <- as.numeric(MeanEsts[,3])
	names(MeanEsts) <- c("Year","Est.Mean","SE")
	MeanEsts$WYear = WYears
# calculate trend in logged mean over time
	if(method=="SLRDB")  fit<-lm(log(Est.Mean) ~ WYear, data=MeanEsts)
	if(method=="WLRDB")  fit<-lm(log(Est.Mean) ~ WYear, weights=1/(MeanEsts$SE^2), data=MeanEsts)

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

	fit<-PWIGLS_Count(Z=getME(fit_PO,"Z"),dat=dat,stage1wt=stage1wt,stage2wt=stage2wt,type=type,stratum=stratum,slope=slope,PO.fitted=predict(fit_PO,type="response"),PO.resid=resid(fit_PO,type="response"))
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
ans.trendCI = exp(CI)-1
names(ans.trendCI) = c("Annual Pct Change", "CI low", "CI high")

if(!method %in% c("SLRDB","WLRDB")) return(list(ModelEstimates = ans.model, TrendTest = ans.trendtest, TrendCI = ans.trendCI))
if(method %in% c("SLRDB","WLRDB")) return(list(ModelEstimates = ans.model, TrendTest = ans.trendtest, TrendCI = ans.trendCI, DBests = data.frame(MeanEsts, Resid = resid(fit))))

}


