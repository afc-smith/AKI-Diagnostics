########################################## AKI DIAGNOSTICS MODEL: HOSPITAL PERIOD  ##############################################################

#This model is for patients receiving the new test who are identified as having AKI S2+ and therefore cannot benefit from early intervention

# Model health states:        1=  In ICU with normal kidney function  (No AKI cohort)
#                             2=  General hospital ward (No AKI)
#                             3=  Discharged home (No AKI)
#                             4=  Dead (no AKI cohort)
#                             5=  In ICU with normal kidney function (AKI cohort) 
#                             6=  In ICU with AKI Severity S1 (KDIGO grading)                                   
#                             7=  In ICU with AKI Severity S2                                                   
#                             8=  In ICU with AKI Severity S3                                                   
#                             9=  In ICU with AKI Severity S3 and receiving Renal Replacement Therapy (RRT)     
#                             10= General Ward post AKI in ICU 
#                             11= General Ward +RRT post AKI in ICU                            
#                             12= Discharged home                                                               
#                             13= Discharged home +RRT                                                          
#                             14= Dead (AKI cohort)                                                                         

# Global parameters of model 
S  <- 14            # Number of health states in hospital model 
T  <- 90            # Time Horizon = 90 days (3 months; cycle lengths of 1 day)

# Create empty objects
tps   <- array(NA,c(S,S,T))             # Model Transition Probability Matrix (defines probability of moving between health states over time)
trace <- array(NA, dim=c(T,S,Nsim))     # Trace: to record proportion of cohort in each health state over time, for each simulation
qtime <- rep(NA,len=T)                  # Record QALYs over time 
cost  <- rep(NA,len=T)                  # Record costs over time 


model <- function(i,test) {                                   
  
#======================= TRANSITION PROBABILITY MATRIX ================================================================================================================

  for (t in 1:T){               
    
    #Transitions from ICU with normal kidney function (no AKI cohort)      
    tps[1,1,t]   <-  1                                                                 
    tps[1,2,t]   <-  0                       
    tps[1,3,t]   <-  0             
    tps[1,4,t]   <-  0 
    tps[1,5,t]   <-  0
    tps[1,6,t]   <-  0                                                 
    tps[1,7,t]   <-  0     
    tps[1,8,t]   <-  0    
    tps[1,9,t]   <-  0 
    tps[1,10,t]  <-  0  
    tps[1,11,t]  <-  0 
    tps[1,12,t]  <-  0 
    tps[1,13,t]  <-  0
    tps[1,14,t]  <-  0
    
    #Transitions from general ward (no AKI cohort)
    # Note: Assumption of no returns to ICU, and hospital aquired AKI is not considered (focus is on ICU-acquired AKI)
    tps[2,1,t]   <-  0
    tps[2,2,t]   <-  1                                      
    tps[2,3,t]   <-  0                                                      
    tps[2,4,t]   <-  0                                                                
    tps[2,5,t]   <-  0
    tps[2,6,t]   <-  0                             
    tps[2,7,t]   <-  0
    tps[2,8,t]   <-  0
    tps[2,9,t]   <-  0 
    tps[2,10,t]  <-  0  
    tps[2,11,t]  <-  0
    tps[2,12,t]  <-  0
    tps[2,13,t]  <-  0
    tps[2,14,t]  <-  0
    
    #Transitions from discharged home (no AKI cohort)
    tps[3,1,t]   <-  0
    tps[3,2,t]   <-  0
    tps[3,3,t]   <-  1        
    tps[3,4,t]   <-  0
    tps[3,5,t]   <-  0
    tps[3,6,t]   <-  0                        
    tps[3,7,t]   <-  0
    tps[3,8,t]   <-  0
    tps[3,9,t]   <-  0 
    tps[3,10,t]  <-  0 
    tps[3,11,t]  <-  0
    tps[3,12,t]  <-  0 
    tps[3,13,t]  <-  0
    tps[3,14,t]  <-  0
    
    #Transitions from mortality (no AKI cohort)
    tps[4,1,t]   <-  0
    tps[4,2,t]   <-  0
    tps[4,3,t]   <-  0
    tps[4,4,t]   <-  1  
    tps[4,5,t]   <-  0                    
    tps[4,6,t]   <-  0                   
    tps[4,7,t]   <-  0                   
    tps[4,8,t]   <-  0         
    tps[4,9,t]   <-  0                   
    tps[4,10,t]  <-  0                  
    tps[4,11,t]  <-  0                               
    tps[4,12,t]  <-  0                                      
    tps[4,13,t]  <-  0                           
    tps[4,14,t]  <-  0
    
    #Transitions from ICU with normal kidney function (AKI cohort)
    # Assumption: All AKI cohort must spend at least one day in an AKI health state (so no transitions to non-ICU/non-AKI health states here)
    tps[5,1,t]   <-  0 
    tps[5,2,t]   <-  0                                        
    tps[5,3,t]   <-  0                                         
    tps[5,4,t]   <-  0
    tps[5,5,t]   <-  tp_norm[i,"S0",t]
    tps[5,6,t]   <-  tp_norm[i,"S1",t]           
    tps[5,7,t]   <-  tp_norm[i,"S2",t]        
    tps[5,8,t]   <-  tp_norm[i,"S3",t]        
    tps[5,9,t]   <-  tp_norm[i,"RRT",t]     
    tps[5,10,t]  <-  tp_norm[i,"ward",t]                     
    tps[5,11,t]  <-  tp_norm[i,"ward_RRT",t]  
    tps[5,12,t]  <-  tp_norm[i,"disch",t]                                     
    tps[5,13,t]  <-  tp_norm[i,"disch_RRT",t]                               
    tps[5,14,t]  <-  tp_norm[i,"mort",t]    
    
    #Transitions from ICU AKI S1
    tps[6,1,t]   <-  0
    tps[6,2,t]   <-  0
    tps[6,3,t]   <-  0
    tps[6,4,t]   <-  0
    tps[6,5,t]   <-  tp_S1[i,"S0",t]                      
    tps[6,6,t]   <-  tp_S1[i,"S1",t]      
    tps[6,7,t]   <-  tp_S1[i,"S2",t]                      
    tps[6,8,t]   <-  tp_S1[i,"S3",t]          
    tps[6,9,t]   <-  tp_S1[i,"RRT",t]           
    tps[6,10,t]  <-  tp_S1[i,"ward",t]                    
    tps[6,11,t]  <-  tp_S1[i,"ward_RRT",t]    
    tps[6,12,t]  <-  tp_S1[i,"disch",t]       
    tps[6,13,t]  <-  tp_S1[i,"disch_RRT",t]   
    tps[6,14,t]  <-  tp_S1[i,"mort",t]        
    
    #Transitions from ICU AKI S2 
    tps[7,1,t]   <-  0      
    tps[7,2,t]   <-  0
    tps[7,3,t]   <-  0
    tps[7,4,t]   <-  0
    tps[7,5,t]   <-  tp_S2[i,"S0",t]                   
    tps[7,6,t]   <-  tp_S2[i,"S1",t]                   
    tps[7,7,t]   <-  tp_S2[i,"S2",t]
    tps[7,8,t]   <-  tp_S2[i,"S3",t]           
    tps[7,9,t]   <-  tp_S2[i,"RRT",t]          
    tps[7,10,t]  <-  tp_S2[i,"ward",t]                    
    tps[7,11,t]  <-  tp_S2[i,"ward_RRT",t]    
    tps[7,12,t]  <-  tp_S2[i,"disch",t]
    tps[7,13,t]  <-  tp_S2[i,"disch_RRT",t]   
    tps[7,14,t]  <-  tp_S2[i,"mort",t]        
    
    #Transitions from ICU AKI S3
    tps[8,1,t]   <-  0     
    tps[8,2,t]   <-  0
    tps[8,3,t]   <-  0
    tps[8,4,t]   <-  0
    tps[8,5,t]   <-  tp_S3[i,"S0",t]                   
    tps[8,6,t]   <-  tp_S3[i,"S1",t]                  
    tps[8,7,t]   <-  tp_S3[i,"S2",t]                  
    tps[8,8,t]   <-  tp_S3[i,"S3",t]
    tps[8,9,t]   <-  tp_S3[i,"RRT",t]          
    tps[8,10,t]  <-  tp_S3[i,"ward",t]                   
    tps[8,11,t]  <-  tp_S3[i,"ward_RRT",t]                  
    tps[8,12,t]  <-  tp_S3[i,"disch",t]  
    tps[8,13,t]  <-  tp_S3[i,"disch_RRT",t]             
    tps[8,14,t]  <-  tp_S3[i,"mort",t]        
    
    #Transitions from ICU RRT
    tps[9,1,t]   <-  0
    tps[9,2,t]   <-  0
    tps[9,3,t]   <-  0
    tps[9,4,t]   <-  0
    tps[9,5,t]   <-  tp_RRT[i,"S0",t]
    tps[9,6,t]   <-  tp_RRT[i,"S1",t]
    tps[9,7,t]   <-  tp_RRT[i,"S2",t]
    tps[9,8,t]   <-  tp_RRT[i,"S3",t]
    tps[9,9,t]   <-  tp_RRT[i,"RRT",t]
    tps[9,10,t]  <-  tp_RRT[i,"ward",t] 
    tps[9,11,t]  <-  tp_RRT[i,"ward_RRT",t]               
    tps[9,12,t]  <-  tp_RRT[i,"disch",t] 
    tps[9,13,t]  <-  tp_RRT[i,"disch_RRT",t]         
    tps[9,14,t]  <-  tp_RRT[i,"mort",t]       

#Transitions from General ward
tps[10,1,t]   <-  0
tps[10,2,t]   <-  0
tps[10,3,t]   <-  0                                     
tps[10,4,t]   <-  0
tps[10,5,t]   <-  tp_ward[i,"S0",t]
tps[10,6,t]   <-  tp_ward[i,"S1",t]
tps[10,7,t]   <-  tp_ward[i,"S2",t]
tps[10,8,t]   <-  tp_ward[i,"S3",t]
tps[10,9,t]   <-  tp_ward[i,"RRT",t]
tps[10,10,t]  <-  tp_ward[i,"ward",t]      
tps[10,11,t]  <-  tp_ward[i,"ward_RRT",t]                  
tps[10,12,t]  <-  tp_ward[i,"disch",t]                    
tps[10,13,t]  <-  tp_ward[i,"disch_RRT",t]    
tps[10,14,t]  <-  tp_ward[i,"mort",t]    

#Transitions from General Ward +RRT
tps[11,1,t]  <-  0
tps[11,2,t]  <-  0
tps[11,3,t]  <-  0
tps[11,4,t]  <-  0
tps[11,5,t]  <-  tp_ward_RRT[i,"S0",t]         
tps[11,6,t]  <-  tp_ward_RRT[i,"S1",t]
tps[11,7,t]  <-  tp_ward_RRT[i,"S2",t]
tps[11,8,t]  <-  tp_ward_RRT[i,"S3",t]
tps[11,9,t]  <-  tp_ward_RRT[i,"RRT",t]
tps[11,10,t] <-  tp_ward_RRT[i,"ward",t] 
tps[11,11,t] <-  tp_ward_RRT[i,"ward_RRT",t]  
tps[11,12,t] <-  tp_ward_RRT[i,"disch",t]                
tps[11,13,t] <-  tp_ward_RRT[i,"disch_RRT",t]
tps[11,14,t] <-  tp_ward_RRT[i,"mort",t]

#Transitions from dishcarged home
tps[12,1,t]  <-  0
tps[12,2,t]  <-  0
tps[12,3,t]  <-  0
tps[12,4,t]  <-  0
tps[12,5,t]  <-  0     
tps[12,6,t]  <-  0    
tps[12,7,t]  <-  0
tps[12,8,t]  <-  0
tps[12,9,t]  <-  0
tps[12,10,t] <-  0 
tps[12,11,t] <-  0 
tps[12,12,t] <-  1- tp_dischmort[i]    
tps[12,13,t] <-  0
tps[12,14,t] <-  tp_dischmort[i]

#Transitions from discharged home + RRT
tps[13,1,t]  <-  0                        
tps[13,2,t]  <-  0
tps[13,3,t]  <-  0
tps[13,4,t]  <-  0
tps[13,5,t]  <-  0
tps[13,6,t]  <-  0
tps[13,7,t]  <-  0
tps[13,8,t]  <-  0
tps[13,9,t]  <-  0
tps[13,10,t] <-  0
tps[13,11,t] <-  0
tps[13,12,t] <-  0
tps[13,13,t] <-  1-tp_dischmort_RRT[i]
tps[13,14,t] <-  tp_dischmort_RRT[i]

#Transitions from death (AKI cohort)
tps[14,1,t]  <-  0
tps[14,2,t]  <-  0
tps[14,3,t]  <-  0
tps[14,4,t]  <-  0
tps[14,5,t]  <-  0
tps[14,6,t]  <-  0
tps[14,7,t]  <-  0
tps[14,8,t]  <-  0
tps[14,9,t]  <-  0 
tps[14,10,t] <-  0 
tps[14,11,t] <-  0 
tps[14,12,t] <-  0
tps[14,13,t] <-  0
tps[14,14,t] <-  1

  }    



#==================== MARKOV trace ========================================================================================================================

# trace <- matrix(nrow = T, ncol = S)

# Define what proportion of the patient cohort start in each health state 

# trace for t = 1

trace [1,1,i]   <-  0                  
trace [1,2,i]   <-  0
trace [1,3,i]   <-  0
trace [1,4,i]   <-  0
trace [1,5,i]   <-  0
trace [1,6,i]   <-  0
trace [1,7,i]   <-  ifelse(test==1, start.S2.Nephro_TP[i], ifelse(test==21, start.S2.NGALp_TP[i], ifelse(test==22, start.S2.NGALu_TP[i],ifelse(test==23, start.S2.NGALs_TP[i], ifelse(test==31,start.S2.CystCp_TP[i],ifelse(test==32,start.S2.CystCu_TP[i],ifelse(test==33,start.S2.CystCs_TP[i],NA)))))))  
trace [1,8,i]   <-  ifelse(test==1, start.S3.Nephro_TP[i], ifelse(test==21, start.S3.NGALp_TP[i], ifelse(test==22, start.S3.NGALu_TP[i],ifelse(test==23, start.S3.NGALs_TP[i], ifelse(test==31,start.S3.CystCp_TP[i],ifelse(test==32,start.S3.CystCu_TP[i],ifelse(test==33,start.S3.CystCs_TP[i],NA))))))) 
trace [1,9,i]   <-  ifelse(test==1, start.RRT.Nephro_TP[i], ifelse(test==21, start.RRT.NGALp_TP[i], ifelse(test==22, start.RRT.NGALu_TP[i],ifelse(test==23, start.RRT.NGALs_TP[i], ifelse(test==31,start.RRT.CystCp_TP[i],ifelse(test==32,start.RRT.CystCu_TP[i],ifelse(test==33,start.RRT.CystCs_TP[i],NA)))))))  
trace [1,10,i]  <-  0  
trace [1,11,i]  <-  0    
trace [1,12,i]  <-  0 
trace [1,13,i]  <-  0
trace [1,14,i]  <-  0



# Define movement of cohort after first cycle according to transition probability matrix 

# trace for t >= 2

for (t in 2:T) {
  trace[t,,i] <- trace[t-1,,i] %*% tps[,,t]             #Matrix multiplication, which starts with multiplication of trace at t=1 by the transition matrix, and loops over time until T is reached                       
}

# Store trace
trace<<-trace      



#======================= Outputs ============================================================================================================================

# Model health states:  1=  ICU with normal kidney function (No AKI); 2= ward (No AKI); 3= Discharged (No AKI); 4= Dead (no AKI cohort) 5= ICU normal kidney function (AKI cohort); 
# 6= AKI S1; 7= AKI S2; 8=  AKI S3; 9= AKI RRT; 10= General Ward; 11= General Ward +RRT; 12= Discharged; 13= Discharged home +RRT; 14= Dead (AKI cohort)  


#### QALYs ####

# Calculate QALYs: trace*utility
# Note: conversion to daily utility (i.e. QALYs) by dividing by 365. 
# Note: Half-cycle correction not applied in this model                

for (t in 1:T) {
  
  qtime[t] <-   ((trace[t,1,i]*u_ICU[i]/365 
                  + trace[t,2,i]*u_ward[i]/365 
                  + trace[t,3,i]*u_disch[i]/365
                  + trace[t,5,i]*u_ICU[i]/365
                  + trace[t,6,i]*u_ICU[i]/365
                  + trace[t,7,i]*u_ICU[i]/365 
                  + trace[t,8,i]*u_ICU[i]/365
                  + trace[t,9,i]*u_ICU[i]/365
                  + trace[t,10,i]*u_ward[i]/365
                  + trace[t,11,i]*u_ward_dialysis[i]/365
                  + trace[t,12,i]*u_disch[i]/365
                  + trace[t,13,i]*u_disch_dialysis[i]/365))  
}

QALYs <- sum(qtime[1:T])               

# Store QALYs
qtime<<-qtime



#### COSTs ####     

#Note- no discounting or half-cycle correction applied in HOSP model. 

for (t in 1:T) {
  
  cost[t] <-   ((trace[t,1,i]*c_ICU[i] 
                 + trace[t,2,i]*c_ward[i]
                 + trace[t,3,i]*c_disch[i]
                 + trace[t,5,i]*c_ICU[i] 
                 + trace[t,6,i]*c_ICU_AKI[i]
                 + trace[t,7,i]*c_ICU_AKI[i]
                 + trace[t,8,i]*c_ICU_AKI[i]
                 + trace[t,9,i]*c_ICU_dialysis[i]
                 + trace[t,10,i]*c_ward[i]
                 + trace[t,11,i]*c_ward_dialysis[i]
                 + trace[t,12,i]*c_disch[i]
                 + trace[t,13,i]*c_disch_dialysis[i]))        
}

cost <- sum(cost[1:T]) +                                                                     
  (trace [1,7,i]+trace [1,8,i]+trace [1,9,i])*   #i.e proportion in this model
  ifelse(test==1, c_Nephro[i], ifelse(test==21 | test==22 | test==23, c_NGAL[i], ifelse(test==31 | test==32 | test==33,c_CystC[i], 0))) #Cost of doing test
     

# Store costs
cost<<-cost 

#Model return:
return(c(QALYs,cost,trace[T,1,i],trace[T,2,i],trace[T,3,i],trace[T,4,i],trace[T,5,i],trace[T,6,i],trace[T,7,i],trace[T,8,i],trace[T,9,i],trace[T,10,i],trace[T,11,i],trace[T,12,i],trace[T,13,i],trace[T,14,i]))


}
