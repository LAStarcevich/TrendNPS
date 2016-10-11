TrendNPS_Cont<-function(dat,method,slope=TRUE,type=NA,stratum=NA,Y,lat=NA,long=NA,stage1wt=NA,stage2wt=NA,str1prop=NA,nbhd=TRUE) {

# dat = data set for trend analysis
# dat must contain fields for Site, WYear (scalar year value), and Year (year factor)
# other covariate fields in dat may be used in the trend model, but WYear must be used to estimate trend

# Calculates a trend model specified by method:
# 'PO' = Piepho and Ogutu (2002) unreplicated linear mixed model 
# 'SLRDB' = 
# 'WLRDB' = 
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
# stage1wt = weight from the sampling design, adjusted for frame and/or nonresponse error if needed
# stage2wt = panel weight from temporal revisit design

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
}		# End PO
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

ans=data.frame(mu,trend,SEtrend,sig2a,sig2t,sigat,sig2b,sig2e,eta)	# ,q.t
names(ans)=c('mu','trend','SEtrend','sig2a','sig2t','sigat','sig2b','sig2e','eta')
rownames(ans) = NULL
return(ans)
}
