PWIGLS_ALL<-function(Z,dat,stage1wt,stage2wt,type,stratum,slope) {
#' @export
#' 
#' @title Probability-weighted iterative generalized least squares (PWIGLS).
#'   
#' @description Calculate trend from data collected with complex survey designs by 
#' incorporating weights with a linear mixed model. This is an internal function 
#' called from TrendNPS_Cont.  
#'   
#' @param Z  Random-effects design matrix from the unweighted (PO) model.
#' 
#' @param dat	 Data frame containing columns at least for \code{Site}, 
#'   \code{WYear}, \code{Year}, and the continuous outcome of interest \code{Y}.
#'   See section "Data frame \code{dat}" below.
#' 
#' @param stage1wt	Design weights from the original sample draw without 
#'   accounting for temporal revisit designs.
#'   
#' @param stage2wt	Panel inclusion weights for each site each year.
#'   
#' @param type	Scaling type when \code{method="PWIGLS"}. Valid values include 
#'   \code{"Aonly"}, \code{"A"}, \code{"AI"}, \code{"B"}, \code{"BI"} 
#'   \code{"C"}. See section "Options for variable \code{type}" below.
#'   
#' @param stratum	 Text string identifying an optional two-level stratification 
#'   factor in \code{dat}. Use stratum = NA to indicate no stratification used.
#'   
#' @param slope	 Logical value indicating inclusion of a random site-level slope
#'  effect in the variance components structure in addition to the Site- and Year-
#'  level random intercept terms. Default = TRUE.
#'
#' @return Returns a vector of regression coefficient estimates for the trend model. 
#'   
#' @details Calculates the probability-weighted iterative generalized least squares
#' (PWIGLS) trend model (Pfeffermann et al. 1998; Asparouhov 2006).
#' 
#' @section Options for variable \code{type}:
#'   
#'   Selection of \code{method="PWIGLS"} requires further specification of 
#'   argument \code{type}.  Valid options include
#'   
#'   \tabular{ll}{
#'   
#'   \code{"Aonly"} \tab Probability weighting but no scaling at either stage \cr
#'   \code{"A"}     \tab Panel-weights scaling with mean site-level design weight  \cr
#'   \code{"AI"}    \tab Panel-weights scaling with mean site-level design weight, but no site-level scaling  \cr
#'   \code{"B"}     \tab Panel-weights scaling with effective mean site-level design weight  \cr
#'   \code{"BI"}    \tab Panel-weights scaling with effective mean site-level design weight, but no site-level scaling  \cr
#'   \code{"C"}     \tab Year-level scaling only with inverse of the average year-level weight \cr
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
#' 	fit<-PWIGLS_ALL(Z=getME(fit_PO,"Z"),dat=dat,stage1wt=stage1wt,stage2wt=stage2wt,type=type,
#'       stratum=stratum,slope=slope)
#' }
#' 
#' 

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
	if(slope) rows.j<-which(Z[,colnames(Z) %in% Site.j][,1]==1)
	if(!slope) rows.j<-which(Z[,colnames(Z) %in% Site.j][,1]==1)

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
}  # end site loop


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
	dat$SiteWt<-rowSums(as.matrix(A))			
	dat$SlopeWt<-rowSums(as.matrix(T))
	dat$YearWt<-rowSums(as.matrix(B))

	if(is.na(stratum)) fit.PWIGLS<-lmer(LogY ~ WYear +(-1+YearWt|Year) +(-1+SiteWt+SlopeWt|Site), data=dat, REML=FALSE)  
	if(!is.na(stratum)) fit.PWIGLS<-lmer(LogY ~ WYear*Stratum +(-1+YearWt|Year) +(-1+SiteWt+SlopeWt|Site), data=dat, REML=FALSE)  

}	# End Random Slope model

if(!slope) {

	if(ma>mb) {
		A<-Z[,1:sitecol]
		B<-Z[,(sitecol+1):d]
	}

	if(ma<=mb) {
		A<-Z[,(mb+1):d]
		B<-Z[,1:mb]
	}

# Sum across rows to pick up each weight 
	dat$SiteWt<-rowSums(as.matrix(A))			
	#dat$SlopeWt<-rowSums(T)
	dat$YearWt<-rowSums(as.matrix(B))
	if(is.na(stratum)) fit.PWIGLS<-lmer(LogY ~ WYear +(-1+YearWt|Year) +(-1+SiteWt|Site), data=dat, REML=FALSE)  
	if(!is.na(stratum)) fit.PWIGLS<-lmer(LogY ~ WYear*Stratum +(-1+YearWt|Year) +(-1+SiteWt|Site), data=dat, REML=FALSE)  

}	# End No Random Slope model

return(fit.PWIGLS)
} 



