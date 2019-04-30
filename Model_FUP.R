#======================= AKI DIAGNOSTICS ACUTE FOLLOW UP PERIOD MODEL ==================================================#

#This model tracks post-hospital stay outcomes, including CKD and general mortality. 

# Model health states:        1= Outpatient follow-up (No AKI cohort) 
#                             2= Dead (No AKI cohort)
#                             3= Outpatient follow-up (AKI cohort)
#                             4= CKD (Stages 1-4)
#                             5= ESRD (Stage 5) no dialysis
#                             6= ESRD + Maintenance Dialysis
#                             7= Transplant
#                             8= Dead


# Global parameters of model 
S  <- 8                              # Number of health states in model = 8 
H  <- 100- round(mean(startage))     # Time Horizon = 100 years - startage = 41 (in base case)
#H <- 5                              # Different time horizons tested in sensitivity analyses
#H <- 20
T  <- H + 1                          # Number of model cycles. Add one for half-cycle correction.                   


#create empty objects
tps   <- array(NA,c(S,S,T))             # Model Transition Probability Matrix (defines probability of moving between health states over time)
trace <- array(NA, dim=c(T,S,Nsim))     # Trace: to record proportion of cohort in each health state over time, for each simulation
qtime <- rep(NA,len=T)                  # Record QALYs over time 
cost  <- rep(NA,len=T)                  # Record costs over time 
# note the end of the first cycle at the end of T=2


# Model test intervention parameter 'test':

        # test = 0 for standard care (base) arm
        # test = 1 for include test


model <- function(i, test) {                                  

  
#======================= TRANSITION PROBABILITY MATRIX ================================================================================================================
  
  for (t in 1:T){      
  
#Transitions from ICU followup (no AKI cohort) state. 
# Note: excess mortality due to ICU stay applied to background mortality in first 2 years post discharge. 
    tps[1,1,t] <- 1-ifelse(t<2, fupmort_yr1[i], ifelse(t<=5, fupmort_yr2[i], mort[round(startage+t),3]))                                          
    tps[1,2,t] <- ifelse(t<2, fupmort_yr1[i], ifelse(t<=5, fupmort_yr2[i], mort[round(startage+t),3]))         
    tps[1,3,t] <- 0
    tps[1,4,t] <- 0
    tps[1,5,t] <- 0
    tps[1,6,t] <- 0
    tps[1,7,t] <- 0
    tps[1,8,t] <- 0

#Transitions from mortality (no AKI cohort). 
    tps[2,1,t] <- 0                                          
    tps[2,2,t] <- 1         
    tps[2,3,t] <- 0
    tps[2,4,t] <- 0
    tps[2,5,t] <- 0
    tps[2,6,t] <- 0
    tps[2,7,t] <- 0                          
    tps[2,8,t] <- 0                   
                         
#Transitions from ICU followup (AKI cohort) state. 
# Assumption: excess mortality due to ICU stay and AKI applied over first 5 years only     
# Assumption: patients have to develop CKD (stages 1-4) before progressing to ESRD states. 
    tps[3,1,t] <- 0
    tps[3,2,t] <- 0
    tps[3,3,t] <- 1- ifelse(t<2, fupmort_yr1[i], ifelse(t<=5, fupmort_yr2[i], mort[round(startage+t),3])) - tp_FUPtoCKD[i]                                                                                
    tps[3,4,t] <- tp_FUPtoCKD[i]                                                   
    tps[3,5,t] <- 0                                        
    tps[3,6,t] <- 0
    tps[3,7,t] <- 0                 
    tps[3,8,t] <- ifelse(t<2, fupmort_yr1[i], ifelse(t<=5, fupmort_yr2[i], mort[round(startage+t),3]))  
                         
#Transitions from CKD stages 1-4 state.
# Assumption: no recovery from CKD state 
    tps[4,1,t] <- 0       
    tps[4,2,t] <- 0       
    tps[4,3,t] <- 0                     
    tps[4,4,t] <- 1 - (tp_CKDtoESRD[i] + tp_CKDtodial[i] + max(tp_CKDtomort[i], mort[round(startage+t),3]))
    tps[4,5,t] <- tp_CKDtoESRD[i]                  
    tps[4,6,t] <- tp_CKDtodial[i]
    tps[4,7,t] <- 0 
    tps[4,8,t] <- max(tp_CKDtomort[i], mort[round(startage+t),3])     
  
#Transitions from ESRD (SKD stage 5), no dialysis
# Assumption: no recovery from ESRD (either remain in one of ESRD states or get transplant, no return to CKD/fup state)
    tps[5,1,t] <- 0       
    tps[5,2,t] <- 0     
    tps[5,3,t] <- 0                   
    tps[5,4,t] <- 0    
    tps[5,5,t] <- 1 - (tp_ESRDtodial[i] + tp_ESRDtotrans[i] + max(tp_ESRDtomort[i], mort[round(startage+t),3]))
    tps[5,6,t] <- tp_ESRDtodial[i]
    tps[5,7,t] <- tp_ESRDtotrans[i]       
    tps[5,8,t] <- max(tp_ESRDtomort[i], mort[round(startage+t),3])     
  
#Transitions from ESRD on maintenance dialysis 
    tps[6,1,t] <- 0
    tps[6,2,t] <- 0
    tps[6,3,t] <- 0
    tps[6,4,t] <- 0
    tps[6,5,t] <- 0                                                 
    tps[6,6,t] <- tp_dialysis[i]     
    tps[6,7,t] <- tp_dialtotrans[i] 
    tps[6,8,t] <- tp_dialtomort[i]   
      
#Transitions from successful transplant
    tps[7,1,t] <- 0
    tps[7,2,t] <- 0
    tps[7,3,t] <- 0
    tps[7,4,t] <- 0
    tps[7,5,t] <- 0
    tps[7,6,t] <- tp_transtodial[i]
    tps[7,7,t] <- tp_transplant[i]   
    tps[7,8,t] <- tp_transtomort[i]    
    
#Death= absorbing state (i.e. no movement out of death). 
    tps[8,1,t] <- 0
    tps[8,2,t] <- 0
    tps[8,3,t] <- 0
    tps[8,4,t] <- 0
    tps[8,5,t] <- 0
    tps[8,6,t] <- 0
    tps[8,7,t] <- 0
    tps[8,8,t] <- 1  
    
  }   
  


#==================== MARKOV trace ========================================================================================================================
  
# trace <- matrix(nrow = T, ncol = S)


# Define what proportion of the patient cohort start in each health state 
# Dependent on testing strategy

# trace for t = 1 (beginning of year 1)
  trace [1,1,i]  <- start.FU[i]    
  trace [1,2,i]  <- start.mort[i]
  trace [1,3,i]  <- start.FU_AKI[i]
  trace [1,4,i]  <- start.CKD[i]
  trace [1,5,i]  <- start.ESRD[i]                
  trace [1,6,i]  <- start.dialysis[i]
  trace [1,7,i]  <- start.transplant[i]
  trace [1,8,i]  <- start.mort_AKI[i]
  

# Define movement of cohort after first cycle according to transition probability matrix    

# trace for t >= 2

  for (t in 2:T) {
    trace[t,,i] <- trace[t-1,,i] %*% tps[,,t]            #Matrix multiplication, which starts with multiplication of trace at t=1 by the transition matrix, and loops over time until T is reached                           
  }


# Store trace
trace<<-trace
  
  

#======================= Outputs ============================================================================================================================


#### QALYs ####
  
# Calculate QALYs: trace*utility

  for (t in 1:T) {
                                                           
    qtime[t] <-   (trace[t,1,i]*ifelse(t<=1,u_FU_1[i], ifelse(t<=4, u_FU_2[i], u_FU_5[i]))     #Assuming here that ICU stay does not have impact on utility once patients have CKD etc. 
                  + trace[t,3,i]*ifelse(t<=1,u_FU_1[i], ifelse(t<=4, u_FU_2[i], u_FU_5[i]))  
                  + trace[t,4,i]*u_CKD[i]
                  + trace[t,5,i]*u_ESRD[i]
                  + trace[t,6,i]*u_dialysis[i]
                  + trace[t,7,i]*u_transplant[i])/((1+disc)^(t-1))  
  }
  
#QALYs <- qtime[1]/2+qtime[T]/2       #USE FOR SA TIME HORIZON=1 ONLY
QALYs <- qtime[1]/2+sum(qtime[2:(T-1)])+qtime[T]/2   #Half-cycle correction. 


# Store qtime
qtime<<-qtime
  


#### COSTs #### 
  
# Calculate costs: trace*cost
  for (t in 1:T) {
    
    cost[t] <-   (trace[t,1,i]*ifelse(t<=10,  c_fup[i,t],                   c_fup_yr11[i])    
                 + trace[t,3,i]*ifelse(t<=10, c_fupAKI[i,t],                c_fup_yr11[i])
                 + trace[t,4,i]*ifelse(t<=10, c_fupAKI[i,t]+c_CKD[i],       c_fup_yr11[i]+ c_CKD[i])
                 + trace[t,5,i]*ifelse(t<=10, c_fupAKI[i,t]+c_ESRD[i],      c_fup_yr11[i]+ c_ESRD[i])
                 
                 + trace[t,6,i]*ifelse(t<=10, c_fupAKI[i,t]+c_dialysis[i],  c_fup_yr11[i]+ c_dialysis[i])
                 + ifelse(t==1, 0, trace[t-1,4,i]*tp_CKDtodial[i]*(c_dialysis_yr1[i]-c_dialysis[i])) 
                 + ifelse(t==1, 0, trace[t-1,5,i]*tp_ESRDtodial[i]*(c_dialysis_yr1[i]-c_dialysis[i]))    
                 + ifelse(t==1, 0, trace[t-1,7,i]*tp_transtodial[i]*(c_dialysis_yr1[i]-c_dialysis[i])) 
                 
                 + trace[t,7,i]*ifelse(t<=10, c_fupAKI[i,t]+c_transplant[i],c_fup_yr11[i]+ c_transplant[i])
                 + ifelse(t==1, 0, trace[t-1,5,i]*tp_ESRDtotrans[i]*(c_transplant_yr1[i]-c_transplant[i]))    
                 + ifelse(t==1, 0, trace[t-1,6,i]*tp_dialtotrans[i]*(c_transplant_yr1[i]-c_transplant[i])))/  ((1+disc)^(t-1))
  }

#Cost  <- cost[1]/2 + cost[T]/2                            
Cost  <- cost[1]/2 + sum(cost[2:(T-1)])+cost[T]/2       #half cycle correction 


# Store  costs
Cost<<-Cost 
  

#Model return:
return(c(QALYs,Cost))
  
}
