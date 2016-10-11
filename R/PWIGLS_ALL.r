#' @export
#' 
#' @title Calculate cool stuff.
#'   
#' @description Calculate cool stuff.  
#'   
#' @param Z  Something.
#' 
#' @param dat	 Data frame containing columns at least for \code{Site}, 
#'   \code{WYear}, \code{Year}, and the continuous outcome of interest \code{Y}.
#'   See section "Data frame \code{dat}" below.
#' 
#' @param type	Scaling type when \code{method="PWIGLS"}. Valid values include 
#'   \code{"Aonly"}, \code{"A"}, \code{"AI"}, \code{"B"}, \code{"BI"} 
#'   \code{"C"}. See section "Options for variable \code{type}" below.
#'   
#' @param slope	 Logical value indicating inclusion of a random site-level slope
#'   effect in the variance components structure used for the \code{"PO"} and 
#'   \code{"PWIGLS"} trend methods.  Site- and year-level random intercept terms
#'   included as default.
#'   
#' @param stratum	 Text string identifying an optional two-level stratification 
#'   factor in \code{dat}.
#'   
#' @param stage1wt	Design weights from the original sample draw without 
#'   accounting for temporal revisit designs.
#'   
#' @param stage2wt	Panel inclusion weights for each site each year.
#'   
#' @return Something cool.
#'   
#' @details Something cool.
#' 
#' @section Options for variable \code{type}:
#'   
#'   Selection of \code{method="PWIGLS"} requires further specification of 
#'   argument \code{type}.  Valid options include
#'   
#'   \tabular{ll}{
#'   
#'   \code{"Aonly"} \tab probability weighting but no scaling at either stage \cr
#'   \code{"A"}     \tab scaling for panel weights with the mean site-level design weight  \cr
#'   \code{"AI"}    \tab scaling for panel weights with the mean site-level design weight, but no site-level scaling  \cr
#'   \code{"B"}     \tab scaling for panel weights with the effective mean site-level design weight  \cr
#'   \code{"BI"}    \tab scaling for panel weights with the effective mean site-level design weight, but no site-level scaling  \cr
#'   \code{"C"}     \tab scaling only at the year level with the inverse of the average year-level weight \cr
#'   
#'   }
#'   
#' @author Leigh Ann Starcevich of Western EcoSystems Technology, Inc.
#'   
#' @references Asparouhov, T. (2006).  General multi-level modeling with 
#'   sampling weights. Communications in Statistics - Theory and Methods 35: 
#'   439-460.
#'
#'   Pfeffermann, D., C.J. Skinner, D.J. Holmes, H. Goldstein, and J. Rasbash 
#'   (1998).  Weighting for unequal selection probabilities in multilevel 
#'   models. Journal of the Royal Statistical Society, Series B 60(1): 23-40.
#'   
#' @seealso \code{TrendNPS_Cont}, \code{LinearizationVar}
#'   
#' @examples 
#' \dontrun{
#' #  ---- Read example data set.
#' SEKIANC_orig = LinearizationVar(Site,wij,xij,eij,varYij)
#' }
#' 
#' 
PWIGLS_ALL<-function(Z,dat,slope,type,stratum,stage1wt,stage2wt) {


SitesTables<-table(dat$Site)
Sites<-unique(as.character(dat$Site))
Years<-sort(unique(dat$Year))
ma=length(Sites)
mb=length(Years)

if(slope) {
	sitecol<-length(unique(dat$Site))*2	    # assuming random slope model
	sitecolindex<-sitecol/2
}
if(!slope) {
	sitecol<-length(unique(dat$Site))	    # assuming no random slope model
	sitecolindex<-sitecol
}
d<-ncol(Z)

# Adjust RE design matrix with sqrt of site design weights
for (j in 1:sitecolindex) {
	Site.j<-Sites[j]
	#rows.j<-which(Z[,colnames(Z) %in% Site.j][,1]==1)
	rows.j<-which(Z[,colnames(Z) %in% Site.j]==1)

if(type=="Aonly") {
# Pfeffermann Step A only from Pfeffermann et al. 1998, no scaling 
	s2j<-s1j<-1
	}

if(type=="AI") {
	n.j<-length(rows.j)
	wj.bar<-sum(dat[rows.j,stage2wt])/n.j
	s2j<-wj.bar
	s1j<-1/wj.bar
	}

if(type=="BI") {
	n0.j<-(sum(dat[rows.j,stage2wt])^2)/sum(dat[rows.j,stage2wt]^2)
	wj.tilde<-sum(dat[rows.j,stage2wt])/n0.j
	s2j<-wj.tilde
	s1j<-1/wj.tilde
	}

if(type=="A") {   # Asperouhov type A = Pfeffermann Step A, scaling 2 = type in Kovacevic
	n.j<-length(rows.j)
	wj.bar<-sum(dat[rows.j,stage2wt])/n.j
	s2j<-1
	s1j<-1/wj.bar
	}

if(type=="B") {   # Asperouhov type B = Pfeffermann Step A, scaling 1
	n0.j<-(sum(dat[rows.j,stage2wt])^2)/sum(dat[rows.j,stage2wt]^2)
	wj.tilde<-sum(dat[rows.j,stage2wt])/n0.j
	s2j<-1
	s1j<-1/wj.tilde
	}

if(type=="C") {   
# Asperouhov type C 
	s2j<-1
	n.j.sum<-sum(SitesTables)	
	wij.sum<-sum(dat[,stage2wt])
	s1j<-n.j.sum/wij.sum
	}

# PSU weight adj
Z[rows.j,]<- Z[rows.j,]/sqrt(dat[rows.j,stage1wt]*s2j)	# multiply sqrt design weight by all cols
# SSU weight adj
if(ma>mb) Z[rows.j,(sitecol+1):d]<-Z[rows.j,(sitecol+1):d]/(sqrt(dat[rows.j,stage2wt]*s1j))  # multiply panel wgt for year effects
if(ma<=mb) Z[rows.j,1:mb]<-Z[rows.j,1:mb]/(sqrt(dat[rows.j,stage2wt]*s1j))  # multiply panel wgt for year effects

if(slope) {
	if(ma>mb) {
		A<-Z[,seq(1,sitecol-1,2)]
		T<-Z[,seq(2,sitecol,2)]
		B<-Z[,(sitecol+1):d]
	}

	if(ma<=mb) {
		A<-Z[,mb+seq(1,sitecol-1,2)]
		T<-Z[,mb+seq(2,sitecol,2)]
		B<-Z[,1:mb]
	}
# Sum across rows to pick up each weight 
	dat$SiteWt<-rowSums(A)			
	dat$SlopeWt<-rowSums(T)
	dat$YearWt<-rowSums(B)

	if(is.na(stratum)) fit.PWIGLS<-lmer(LogY ~ WYear +(-1+YearWt|Year) +(-1+SiteWt+SlopeWt|Site), data=dat, REML=FALSE)  
	if(!is.na(stratum)) fit.PWIGLS<-lmer(LogY ~ WYear*Stratum +(-1+YearWt|Year) +(-1+SiteWt+SlopeWt|Site), data=dat, REML=FALSE)  

	#VarCor<-VarCorr(fit.PWIGLS)
	#varcor<- c(VarCor$Year[1], VarCor$Site[1,1], VarCor$Site[2,2], attr(VarCor,"sc")^2, VarCor$Site[2,1])
}	# End Random Slope model

if(!slope) {
	if(ma>mb) {
		A<-Z[,1:sitecol]
		#T<-Z[,seq(2,sitecol,2)]
		B<-Z[,(sitecol+1):d]
	}

	if(ma<=mb) {
		A<-Z[,(mb+1):d]
		#T<-Z[,mb+seq(2,sitecol,2)]
		B<-Z[,1:mb]
	}
# Sum across rows to pick up each weight 
	dat$SiteWt<-rowSums(A)			
	#dat$SlopeWt<-rowSums(T)
	dat$YearWt<-rowSums(B)
	if(is.na(stratum)) fit.PWIGLS<-lmer(LogY ~ WYear +(-1+YearWt|Year) +(-1+SiteWt|Site), data=dat, REML=FALSE)  
	if(!is.na(stratum)) fit.PWIGLS<-lmer(LogY ~ WYear*Stratum +(-1+YearWt|Year) +(-1+SiteWt|Site), data=dat, REML=FALSE)  

	#VarCor<-VarCorr(fit.PWIGLS)
	#varcor<- c(VarCor$Year[1], VarCor$Site, 0, attr(VarCor,"sc")^2, 0)

}	# End No Random Slope model

}  # end site loop

return(fit.PWIGLS)
} 



