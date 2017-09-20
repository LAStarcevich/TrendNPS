PWIGLS_Count<-function(Z, dat, stage1wt, stage2wt, type=type,stratum=stratum,slope=slope,PO.fitted,PO.resid) {

# Calculate probability-weighted estimated of GLMM by Laplace approximation
# Calls PWIGLS_ALL
# las 5.2.17

#dat$Stratum = as.factor(dat$Stratum)

# Laplace Approximation for Linearization Var
# Calc nu, on link scale
dat$LogY = log(PO.fitted) + (1/PO.fitted)* PO.resid   

if(is.na(stratum)) {
# Using PWIGLS Poisson predictions
fit_PO_Laplace = lmer(LogY ~ WYear  + (1|Year) +(1+WYear|Site), data = dat)  

Z_PO_L = getME(fit_PO_Laplace,"Z")
fit_PWIGLS_Laplace = PWIGLS_ALL(Z=Z_PO_L,dat=dat,stage1wt=stage1wt,stage2wt=stage2wt,type=type,stratum=NA,slope=slope)
}	# End !StRS

if(!is.na(stratum)) {

# Using PWIGLS Poisson predictions
fit_PO_Laplace = lmer(LogY ~ WYear*Stratum  + (1|Year) +(1+WYear|Site), data = dat)  

Z_PO_L = getME(fit_PO_Laplace,"Z")
fit_PWIGLS_Laplace = PWIGLS_ALL(Z=Z_PO_L,dat=dat,stage1wt=stage1wt,stage2wt=stage2wt,type=type,stratum="Stratum",slope=slope)
}	# End StRS

return(fit_PWIGLS_Laplace)
} 

