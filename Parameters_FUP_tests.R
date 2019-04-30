########################################## AKI DIAGNOSTICS MODEL: FOLLOW UP PERIOD MODEL TEST PARAMETERS  ###############################################################

#WARNING: This code requires that the hospital period model has already been run. 
#This is required for the follow-up model health state starting distribution parameters, which are derived from the hospital model trace. 

parameter_draw_FUP  <- function() { 
  
  set.seed(seed)
  
#Define Hospital trace dependent on test

trace_hosp <- array(NA, c(Nsim,14))
for (i in 1:Nsim){
  for (t in 1:14){
    trace_hosp[i,t]<-ifelse(test==1, trace.D90.Nephro[i,t], ifelse(test==21, trace.D90.NGALp[i,t], ifelse(test==22, trace.D90.NGALu[i,t], ifelse(test==23, trace.D90.NGALs[i,t], ifelse(test==31, trace.D90.CystCp[i,t], ifelse(test==32, trace.D90.CystCu[i,t],ifelse(test==33, trace.D90.CystCs[i,t],NA)))))))
}}


start.FU            <- rep(NA,times=Nsim)
start.mort          <- rep(NA,times=Nsim)
start.CKD           <- rep(NA,times=Nsim)
start.mort_AKI      <- rep(NA,times=Nsim)
start.FU_AKI        <- rep(NA,times=Nsim)
start.ESRD          <- rep(NA,times=Nsim)
start.dialysis      <- rep(NA,times=Nsim) 
start.transplant    <- rep(NA,times=Nsim) 

for (i in 1:Nsim) {
  start.FU[i]            <- sum(trace_hosp[i,1], trace_hosp[i,2], trace_hosp[i,3])                                              #sum of hospital model no AKI cohort states
  start.mort[i]          <- sum(trace_hosp[i,4]) 
  start.CKD[i]           <- sum(trace_hosp[i,9], trace_hosp[i,11], trace_hosp[i,13])                                            #sum of hospital states: in ICU+RRT, in ward+RRT, discharged+RRT 
  start.ESRD[i]          <- 0
  start.dialysis[i]      <- 0
  start.transplant[i]    <- 0
  start.mort_AKI[i]      <- trace_hosp[i,14]                                                                            #hospital model mortality state
  start.FU_AKI[i]        <- sum(trace_hosp[i,5], trace_hosp[i,6], trace_hosp[i,7], trace_hosp[i,8], trace_hosp[i,10],trace_hosp[i,12])      #remaining states
}

start.FU            <<- start.FU
start.mort          <<- start.mort
start.CKD           <<- start.CKD
start.mort_AKI      <<- start.mort_AKI
start.FU_AKI        <<- start.FU_AKI
start.ESRD          <<- start.ESRD
start.dialysis      <<- start.dialysis 
start.transplant    <<- start.transplant 


}