#' @export
#' 
#' @title Calculate linearization variance of trend (Skinner, et al. 1989)
#'   
#' @description Calculate the linearization variance of the trend coefficient
#'   
#' @param Site  vector of PSUs.
#' @param wij  inclusion weight = PSU wt * SSU wt.
#' @param xij  year var for trend estimation.
#' @param eij  residual.
#' @param varYij  total variance for yij, estimated from PWIGLS estimates of
#'   variance from P&O.
#' @param str1prop 	Proportion of the first stratum in the population in a 
#'   two-level stratification variable.
#'   
#' @return Something cool.
#'   
#' @details Something cool.
#'   
#' @author Leigh Ann Starcevich of Western EcoSystems Technology, Inc.
#'   
#' @references Skinner 1989., p. 82-83.
#'   
#' @seealso \code{TrendNPS_Cont}, \code{LinearizationVar}
#'   
#' @examples 
#' \dontrun{
#' #  ---- Read example data set.
#' SEKIANC_orig = LinearizationVar_StRS(Site,wij,xij,eij,varYij,str1prop)
#' }
#' 
#' 
LinearizationVar_StRS <- function (Site,wij,xij,eij,varYij,str1prop) {

ILij<-xij[,2]
xij<-xij[,1]

dTdB0 = matrix(c(sum(wij/varYij),
sum(xij*wij/varYij),
sum(eij*wij*(xij^2)/(varYij^2)),
sum(xij*wij*eij/(varYij^2)),
sum(eij*wij/(varYij^2)),
sum(wij*ILij/varYij),
sum(xij*wij*ILij/varYij)),1,7)

dTdB1 = matrix(c(sum(xij*wij/varYij),
sum((xij^2)*wij/varYij),
sum(eij*wij*(xij^3)/(varYij^2)),
sum((xij^2)*wij*eij/(varYij^2)),
sum(xij*eij*wij/(varYij^2)),
sum(xij*wij*ILij/varYij),
sum((xij^2)*wij*ILij/varYij)),1,7)

dTdsig2t = matrix(c(sum(eij*wij*(xij^2)/(varYij^2)),
sum(wij*(xij^3)*eij/(varYij^2)),
sum((-(xij^4)*wij/(2*(varYij^2)))+(wij*(xij^4)*(eij^2)/(varYij^3))),
sum((-(xij^3)*wij/(2*(varYij^2)))+(wij*(xij^3)*(eij^2)/(varYij^3))),
sum((-(xij^2)*wij/(2*(varYij^2)))+((xij^2)*wij*(eij^2)/(varYij^3))),
sum(eij*wij*(xij^2)*ILij/(varYij^2)),
sum(wij*(xij^3)*eij*ILij/(varYij^2))),1,7) 

dTdsigat = matrix(c(sum(eij*xij*wij/(varYij^2)),
sum(wij*(xij^2)*eij/(varYij^2)),
sum((-(xij^3)*wij/(2*(varYij^2)))+((xij^3)*wij*(eij^2)/(varYij^3))),
sum((-(xij^2)*wij/(2*(varYij^2)))+((xij^2)*wij*(eij^2)/(varYij^3))),
sum((-wij*xij/(2*(varYij^2)))+(wij*xij*(eij^2)/(varYij^3))),
sum(eij*xij*wij*ILij/(varYij^2)),
sum(wij*(xij^2)*eij*ILij/(varYij^2))),1,7)

dTdsig2abe = matrix(c(sum(eij*wij/(varYij^2)),
sum(xij*wij*eij/(varYij^2)),
sum((-(xij^2)*wij/(2*(varYij^2)))+(wij*(xij^2)*(eij^2)/(varYij^3))),
sum((-wij*xij/(2*(varYij^2)))+(wij*xij*(eij^2)/(varYij^3))),
sum((-wij/(2*(varYij^2)))+(wij*(eij^2)/(varYij^3))),
sum(eij*wij*ILij/(varYij^2)),
sum(xij*wij*eij*ILij/(varYij^2))) ,1,7)

dTdB2 = matrix(c(sum(wij*ILij/varYij),
sum(xij*wij*ILij/varYij),
sum(eij*wij*(xij^2)*ILij/(varYij^2)),
sum(xij*wij*eij*ILij/(varYij^2)),
sum(eij*wij*ILij/(varYij^2)),
sum(wij*ILij*ILij/varYij),
sum(xij*wij*ILij*ILij/varYij)),1,7)

dTdB3 = matrix(c(sum(xij*wij*ILij/varYij),
sum((xij^2)*wij*ILij/varYij),
sum(eij*wij*(xij^3)*ILij/(varYij^2)),
sum((xij^2)*wij*eij*ILij/(varYij^2)),
sum(xij*eij*wij*ILij/(varYij^2)),
sum(xij*wij*ILij*ILij/varYij),
sum((xij^2)*wij*ILij*ILij/varYij)),1,7)

I.mat = rbind(dTdB0,dTdB1,dTdsig2t,dTdsigat,dTdsig2abe,dTdB2,dTdB3)

#print(dim(I.mat))
I.mat.inv = solve(I.mat)

g = aggregate(list(g1 = eij*wij/varYij, 
g2 = xij*eij*wij/varYij, 
g3 = (((-wij*(xij^2)/(2*varYij))) + (wij*(xij^2)*(eij^2))/(2*(varYij^2))),
g4 = (((-wij*xij/(2*varYij))) + (wij*xij*(eij^2))/(2*(varYij^2))),
g5 = (((-wij/(2*varYij))) + (wij*(eij^2))/(2*(varYij^2))),
g6 = eij*wij*ILij/varYij, 
g7 = xij*eij*wij*ILij/varYij), 
list(Site=Site), sum)

#print('nrow(g)')
#print(nrow(g))

VarL.T = nrow(g)* var(g[,2:8])

#print('I.mat.inv')
#print(I.mat.inv)

VarL.Beta = I.mat.inv%*%VarL.T%*%I.mat.inv

# Use a contrast to get the population level trend variance across two strata
C<-matrix(c(0,1,0,(1-str1prop)),1,4)	# weight to get popn trend (note beta4 is a difference in slopes)
VarC<-C%*%VarL.Beta[c(1:2,6:7),c(1:2,6:7)]%*%t(C)
return(VarC)

}


