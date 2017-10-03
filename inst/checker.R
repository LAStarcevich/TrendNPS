


#   ---- Add in new datasets.  
# load("C:/Users/jmitchell/Downloads/TrendNPS_BinaryCount_DataFrames.RData")
# 
# devtools::use_data(Cover_Acro)
# 
# devtools::use_data(Seastar)
# 
# devtools::document()
# 
# devtools::check(manual=TRUE)


#   ---- Jason testing after package compilation.  

#   ---- Check lazy loading of R dataframes. 

rm(list=ls())
sessionInfo()

detach("package:TrendNPS", unload=TRUE)
require(TrendNPS)

data(Cover_Acro)
data(Seastar)
data(GRTSsample)
data(LakeData_orig)
data(LakeDataRep)

rm(list=c("Cover_Acro","Seastar","GRTSsample","LakeData_orig","LakeDataRep"))


#   ---- Check for csv loading.  Use function system.file to speed up the problem of 
#   ---- finding the user's local installation of the TrendNPS R package.  These 
#   ---- files are stored in the "extdata" folder of the compiled package.  
Cover_Acro <- read.csv(system.file("extdata","Cover_Acro.csv",package="TrendNPS",mustWork=TRUE))
Seastar <- read.csv(system.file("extdata","Seastar.csv",package="TrendNPS",mustWork=TRUE))
GRTSsample <- read.csv(system.file("extdata","GRTSsample.csv",package="TrendNPS",mustWork=TRUE))
LakeData_orig <- read.csv(system.file("extdata","LakeData_orig.csv",package="TrendNPS",mustWork=TRUE))
LakeDataRep <- read.csv(system.file("extdata","LakeDataRep.csv",package="TrendNPS",mustWork=TRUE))

rm(list=c("Cover_Acro","Seastar","GRTSsample","LakeData_orig","LakeDataRep"))












# Trend NPS Binary Examples

# ---- Read example data set.
data(Cover_Acro)
# ---- Load dependent packages.
pkgList <- c("lme4","lmerTest","spsurvey")
TrendNPS_Binary 15
inst <- pkgList %in% installed.packages()
if (length(pkgList[!inst]) > 0) install.packages(pkgList[!inst])
lapply(pkgList, library, character.only = TRUE)
###########################
# Example 1: Trend analysis of a binary indicator of Acrosophonia presence
# with the PO approach: stratification and full random effects model.
TrendAcro_PO_StRS = TrendNPS_Binary(alpha=0.1,
                                    dat=Cover_Acro,method="PO",slope=TRUE,type=NA,stratum="Park",Y="Y",
                                    stage1wt="wgt",stage2wt="PanelWt",str1prop=0.13227)
# $ModelEstimates
# mu trend SEtrend sig2a sig2t sigat sig2b
# -0.2992415 0.221928 0.09341594 0.2867116 0.01394005 -0.04732509 0.2411157
# $TrendTest
# trend SEtrend z-stat pvalue
# 0.221928 0.09341594 2.375698 0.0175158
# $TrendCIofOddsRatio
# Annual Pct Change CI low CI high
# 0.2484815 0.07065702 0.4558408
# Example 2: Trend analysis of a binary indicator of Acrosophonia presence
# with the SLRDB approach for stratification.
TrendAcro_SLRDB_StRS = TrendNPS_Binary(alpha=0.1,
                                       dat=Cover_Acro,method="SLRDB",slope=TRUE,type=NA,stratum="Park",Y="Y",
                                       lat="Lat",long="Long",stage1wt="wgt",stage2wt="PanelWt",
                                       str1prop=0.13227)
# TrendAcro_SLRDB_StRS
# $ModelEstimates
# mu trend SEtrend sig2a sig2t sigat sig2b
# -0.9909923 0.211127 0.08966504 0 0 0 0
#
# $TrendTest
# trend SEtrend z-stat pvalue
# 0.211127 0.08966504 2.354619 0.01854169
#
# $TrendCIofOddsRatio
# Annual Pct Change CI low CI high
# 0.2350692 0.06570992 0.4313426
#
# $DBests
# Year Est.Mean SE WYear Resid
# 2008 0.1604049 0.06934590 0 -0.6642261
# 2009 0.2652986 0.05963881 1 -0.2387430
# 2010 0.4121248 0.10439841 2 0.2135496
# 2011 0.6416667 0.07341420 3 0.9402165
# 2012 0.5439923 0.11698602 4 0.3229097
# 2013 0.5875758 0.03719284 5 0.2893100
# 2014 0.4954618 0.06389670 6 -0.2939231
# 2015 0.4794624 0.04936033 7 -0.5690937
# Plot trend on log-odds scale

num <- TrendAcro_SLRDB_StRS$DBests$Est.Mean
denom <- 1-TrendAcro_SLRDB_StRS$DBests$Est.Mean 
TrendAcro_SLRDB_StRS$DBests$LogOdds = num / denom
plot(TrendAcro_SLRDB_StRS$DBests$Year,
     TrendAcro_SLRDB_StRS$DBests$LogOdds,
     xlab="Year", ylab="Odds ratio of cover proportion")
lines(2008:2015, exp(-0.9909923 + 0.211127*(0:7)), col=2)

# Plot trend on proportional cover scale
plot(TrendAcro_SLRDB_StRS$DBests$Year,
     TrendAcro_SLRDB_StRS$DBests$Est.Mean, 
     xlab="Year", ylab="Cover proportion")
lines(2008:2015, TrendNPS::expit(-0.9909923 + 0.211127*(0:7)), col=2)









# Count (not Cont)

# ---- Read example data set.
data(LakeDataRep)
# ---- Load dependent packages.
pkgList <- c("lme4","lmerTest","spsurvey")
inst <- pkgList %in% installed.packages()
if (length(pkgList[!inst]) > 0) install.packages(pkgList[!inst])
lapply(pkgList, library, character.only = TRUE)
###########################
# Example 1: Trend analysis with the PO approach: stratification and full random
# effects model.
str1prop <- 808.6936/(808.6936+362.8214)
PO_ests <- TrendNPS_Cont(alpha=0.1,dat=LakeDataRep, method="PO", slope=TRUE,
                         stratum="Park", Y="Cl", str1prop=str1prop)
PO_ests

# $ModelEstimates
# mu trend SEtrend sig2a sig2t sigat sig2b
# 1.269532 -0.07196304 0.05624014 0.2846501 0.009663692 -0.05208829 0.04079863
# sig2e eta
# 0.08133577 24.94572
# $TrendTest
# trend SEtrend t-stat eta pvalue
# -0.07196304 0.05624014 -1.279567 24.94572 0.2124709
# $TrendCI
# Annual Pct Change CI low CI high
# -0.06943471 -0.1546776 0.0244041
###########################
# Example 2: Trend analysis with the WLRDB approach: stratification and full random effects model.
str1prop <- 808.6936/(808.6936+362.8214)
WLRDB_ests <- TrendNPS_Cont(alpha=0.1,dat=LakeDataRep, method="WLRDB",
                            slope=TRUE, stratum="Park", Y="Cl", lat="ycoord",
                            long="xcoord",stage1wt="AdjWgt", stage2wt="PanelWt",
                            str1prop=str1prop)
WLRDB_ests
# $ModelEstimates
# mu trend SEtrend sig2a sig2t sigat sig2b sig2e eta
# 0.4218735 0.06037018 0.07802896 0 0 0 0 3.402272 4
# $TrendTest
# trend SEtrend t-stat eta pvalue
# 0.06037018 0.07802896 0.7736894 4 0.4822974
# $TrendCI
# Annual Pct Change CI low CI high
# 0.06222969 -0.1005534 0.2544735
#$DBests
# Year Est.Mean SE WYear Resid
# 2008 4.505724 0.55610529 0 1.083475074
# 2009 2.162006 0.20158290 1 0.288793000
# 2010 1.712219 0.03578556 2 -0.004823561
# 2011 1.592197 0.07584161 3 -0.137869385
# 2012 2.611047 0.14932238 4 0.296397135
# 2013 2.348186 0.18412120 5 0.129918636
###########################
# Example 3: Trend analysis with the PWIGLS method with type="Aonly".
str1prop <- 808.6936/(808.6936+362.8214)
PWIGLS_Aonly_ests <- TrendNPS_Cont(alpha=0.1,dat=LakeDataRep, method="PWIGLS",
                                   type="Aonly", slope=TRUE, stratum="Park", Y="Cl", lat="ycoord",
                                   long="xcoord", stage1wt="AdjWgt", stage2wt="PanelWt",
                                   str1prop=str1prop)
PWIGLS_Aonly_ests
# $ModelEstimates
# mu trend SEtrend sig2a sig2t sigat sig2b
# 1.441448 -0.08109772 0.09649006 10.66269 0.3102913 -1.818939 0.4271289
# sig2e eta

# 0.09145612 24.94572
# $TrendTest
# trend SEtrend t-stat eta pvalue
# -0.08109772 0.09649006 -0.8404774 24.94572 0.4086247
# $TrendCI
# Annual Pct Change CI low CI high
# -0.07789642 -0.2180231 0.08734037








# Trend Cont (not Count)

# ---- Read example data set.
data(Seastar)
# ---- Load dependent packages.
pkgList <- c("lme4","lmerTest","spsurvey")
inst <- pkgList %in% installed.packages()
if (length(pkgList[!inst]) > 0) install.packages(pkgList[!inst])
lapply(pkgList, library, character.only = TRUE)
###########################
# Example 1: Trend analysis of annual site-level counts of leather sea stars
# with the PO approach: stratification and full random effects model.
TrendSeaster_PO_StRS = TrendNPS_Count(alpha=0.1,
                                      dat=Seastar,method="PO",slope=TRUE,type=NA,stratum="Park",Y="Count",
                                      stage1wt="wgt",stage2wt="PanelWt",str1prop=0.13227)
TrendSeaster_PO_StRS
# $ModelEstimates
# mu trend SEtrend sig2a sig2t sigat sig2b

# -0.2992415 0.221928 0.09341594 0.2867116 0.01394005 -0.04732509 0.2411157
# $TrendTest
# trend SEtrend z-stat pvalue
# 0.221928 0.09341594 2.375698 0.0175158
# $TrendCIofOddsRatio
# Annual Pct Change CI low CI high
# 0.2484815 0.07065702 0.4558408
# Example 2: Trend analysis of annual site-level counts of leather sea stars
# with the SLRDB approach for stratification.
TrendSeastar_SLRDB_StRS = TrendNPS_Count(alpha=0.1,
                                         dat=Seastar,method="SLRDB",slope=TRUE,type=NA,stratum="Park",Y="Count",
                                         lat="Lat",long="Long", stage1wt="wgt",stage2wt="PanelWt",str1prop=0.13227)
# TrendSeastar_SLRDB_StRS
# $ModelEstimates
# mu trend SEtrend sig2a sig2t sigat sig2b
# 2.687804 0.0317377 0.04920822 0 0 0 0
#
# $TrendTest
# trend SEtrend z-stat pvalue
# 0.0317377 0.04920822 0.6449674 0.5189483
#
# $TrendCI
# Annual Pct Change CI low CI high
# 0.03224671 -0.04801179 0.1192715
#
# $DBests
# Year Est.Mean SE WYear Resid
# 2009 11.66248 2.268814 0 -0.23142756
# 2010 14.66567 4.057595 1 -0.03403269
# 2011 18.60000 6.206449 2 0.17188195
# 2012 17.17467 5.045078 3 0.06041852
# 2013 21.71013 5.658295 4 0.26302412
# 2014 20.29531 4.128595 5 0.16389719
# 2015 11.99475 2.917410 6 -0.39376154
# Plot trend on original scale
plot(TrendSeastar_SLRDB_StRS$DBests$Year,
     TrendSeastar_SLRDB_StRS$DBests$Est.Mean,
     xlab="Year", ylab="Mean Sea Star Counts")
lines(2009:2015, exp(2.687804 + 0.0317377*(0:6)), col=2)


