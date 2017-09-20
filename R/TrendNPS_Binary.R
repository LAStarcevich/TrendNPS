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
	SE <- NULL  # Jason adds to fool devtools::check wrt "SE" on next line. 
	if(method=="WLRDB")  fit<-lm(log(Est.Mean) ~ WYear, weights=1/(SE^2), data=MeanEsts)

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
ans.trendCI = exp(CI)-1
names(ans.trendCI) = c("Annual Pct Change", "CI low", "CI high")

if(!method %in% c("SLRDB","WLRDB")) return(list(ModelEstimates = ans.model, TrendTest = ans.trendtest, TrendCI = ans.trendCI))
if(method %in% c("SLRDB","WLRDB")) return(list(ModelEstimates = ans.model, TrendTest = ans.trendtest, TrendCI = ans.trendCI, DBests = data.frame(MeanEsts, Resid = resid(fit))))

}


