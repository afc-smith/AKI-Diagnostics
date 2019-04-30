########################################## AKI DIAGNOSTICS MODEL: FOLLOW UP PERIOD MODEL PARAMETERS  ###############################################################

#WARNING: This code requires that the hospital period model has already been run. 
#This is required for the follow-up model health state starting distribution parameters, which are derived from the hospital model trace. 

parameter_draw_FUP  <- function() { 
  
  set.seed(seed)
  
#========================================================================================================================================#
#====================================================== BASELINE PARAMETERS =============================================================#
#========================================================================================================================================#


#====================== CODE FOR FUNCTIONS USED ========================# 

beta.draw <- function(n, m, s) {
  v <- s ^ 2
  rbeta(n, (m ^ 2 - m ^ 3 - m * v) / v,
        (m - 2 * m ^ 2 + m ^ 3 - v + m * v) /v)
}

gamma.draw <- function(n, m, s) {
  v <- s ^ 2
  rgamma(n, m ^ 2 / v, m / v)
}

lnorm.draw  <- function(n, m, s){    
  v <- s ^ 2
  rlnorm(n, log(m) - 0.5 * (log (1 + v / (m ^ 2))), 
         sqrt (log (1 + v / (m ^ 2))))
}


#=========================== PATIENT CHARACTERISTICS ==============================#
  
#Model cohort starting age
# Base case: from Leeds AKI Registry data. Mean age of allcomers to ICU= 63 
# Add 90 days to account for hospital model period
startage <<- 63 + 90/365
  

  
#======================= MODEL STRUCTURAL PARAMETERS ============================#
  
#Time horizon of model
H <- 100- round(mean(startage))     # Time horizon of followup model= 100 - startage (annual cycles) i.e. lifetime horizon capped at 100 years
T <- H+2                            # Add 2 to allow for half-cycle correction       
  
  
  
#============================= PROBABILITIES ======================================# 

#STARTING PROBABILITIES   
#Starting health state probabilities derived from distribution of patients at the end of the hospital period model
#Assumptions: any patient still recieving RRT at end of hospital model is classified as having CKD. 
#All other patients assumed to move to relevant AKI/No AKI follow-up state. 

# Hosptial period Model health states:        1=  In ICU with normal kidney function  (No AKI cohort)
#                                             2=  Transferred to general ward (No AKI)
#                                             3=  Discharged home (No AKI)
#                                             4=  Dead (no AKI cohort)
#                                             5=  In ICU with normal kidney function (AKI cohort)
#                                             6=  In ICU with AKI Severity S1 (KDIGO grading)
#                                             7=  In ICU with AKI Severity S2
#                                             8=  In ICU with AKI Severity S3 
#                                             9=  In ICU with AKI Severity S3 and receiving Renal Replacement Therapy (RRT)
#                                             10=  Transferred to General Ward                          
#                                             11= Transferred to General Ward +RRT
#                                             12= Discharged home 
#                                             13= Discharged home +RRT
#                                             14= Dead (AKI cohort) 

start.FU            <- rep(NA,times=Nsim)
start.mort          <- rep(NA,times=Nsim)
start.CKD           <- rep(NA,times=Nsim)
start.mort_AKI      <- rep(NA,times=Nsim)
start.FU_AKI        <- rep(NA,times=Nsim)
start.ESRD          <- rep(NA,times=Nsim)
start.dialysis      <- rep(NA,times=Nsim) 
start.transplant    <- rep(NA,times=Nsim) 

for (i in 1:Nsim) {
start.FU[i]            <- sum(trace_hosp.D90[i,1], trace_hosp.D90[i,2], trace_hosp.D90[i,3])               #TRACE HOSP IS GENERATED FROM HOSPITAL MODEL RUN                                  #sum of hospital model no AKI cohort states
start.mort[i]          <- sum(trace_hosp.D90[i,4]) 
start.CKD[i]           <- sum(trace_hosp.D90[i,9], trace_hosp.D90[i,11], trace_hosp.D90[i,13])                                            #sum of hospital states: in ICU+RRT, in ward+RRT, discharged+RRT 
start.ESRD[i]          <- 0
start.dialysis[i]      <- 0
start.transplant[i]    <- 0
start.mort_AKI[i]      <- trace_hosp.D90[i,14]                                                                            #hospital model mortality state
start.FU_AKI[i]        <- sum(trace_hosp.D90[i,5], trace_hosp.D90[i,6], trace_hosp.D90[i,7], trace_hosp.D90[i,8], trace_hosp.D90[i,10],trace_hosp.D90[i,12])      #remaining states
}

start.FU            <<- start.FU
start.mort          <<- start.mort
start.CKD           <<- start.CKD
start.mort_AKI      <<- start.mort_AKI
start.FU_AKI        <<- start.FU_AKI
start.ESRD          <<- start.ESRD
start.dialysis      <<- start.dialysis 
start.transplant    <<- start.transplant 

#Validity check (should all sum to 1)
#sum(start.FU, start.CKD,start.ESRD, start.dialysis,start.transplant, start.mort,start.FU_AKI, start.mort_AKI)



#============================ TRANSITION PROBABILITIES ============================#

#====== Follow-up Mortality ======#

#Follow-up mortality for ICU cohort 
#Base case: from Lone et al. Reported 1-year mortality 10.9%, 95% CI: 10.0 to 11.7. 
#(0.117-0.10)/(3.92)  #sd= 0.004336735
#estBetaParams(0.109,0.004336735)  #alpha=0.006706571, beta=0.0548216
fupmort_yr1 <<- rbeta(Nsim,2.332 , 19.0625)
#hist(fupmort_yr1)

#Lone et al. Reported 5 year mortality= 32.3% 95% CI: 31 to 33.6. So at end of year 1 89.1% are still alive, 
# and by year 5 out of those still alive at yr 1, an extra (32.3-10.9)/0.891 = 24.01796% died. (95% CI: 0.2255892 to 0.2547699)
#convert 4 year mortality in to 1 year mortality. 
#1-exp(-(-log(1-0.2401796)/4)*1) #=0.06636368
#1-exp(-(-log(1-0.2255892)/4)*1) #=0.06191358
#1-exp(-(-log(1-0.2547699)/4)*1) #=0.07087831
##sd 
#(0.07087831-0.06191358)/3.92  #0.002286921
#estBetaParams(0.06636368,0.002286921)  #0.4860051, 6.837354
fupmort_yr2<<- rbeta(Nsim, 1.731627,24.36137)
#hist(fupmort_yr2, breaks=100)

#Background mortality 
#Mortality for normal population, from ONS 2014      
mort <<- read.csv("Background_mortality.csv")    
#nb- daily mort in second column, annual mort in 3rd column. 




#====== CKD ======#

#Risk of CKD from follow-up post AKI 
#fROM Rimes-Stigare et al 2015. Annual rate CKD in no AKI ICU cohort=0.44% (0.0044). Adjusted Incident Rate Ratio for AKI cohort of 7.6 (95% CI:  5.5 to 10.4)
RR_CKD <<- lnorm.draw(Nsim,7.6,(10.4-5.5)/3.92)

tp_FUPtoCKD  <<- 1-exp(log(1-0.0044*RR_CKD))-0.0044         
#hist(tp_FUPtoCKD)
#hist(RR_CKD)


#Mortality from CKD (1-4)    
#Base case: from Kent 2015. 
#Number of patient years with events: For CKD 1-3 (n=1494): vascular deaths 36 (1.3%), non-vasc 97 (1.6%).  For CKD 4 (n=2228)= 92 (1%), 219 (2.5%). 
tp_CKDtomort  <<- (1494/(1494+2228))*(rbeta(Nsim,36,(36/0.013)-36) +rbeta(Nsim,97, (97/0.016)-97)) +(2228/(1494+2228))*(rbeta(Nsim,92,(92/0.01)-92) + rbeta(Nsim,219,(219/0.025)-219))     
#hist(tp_CKDtomort)
#summary(tp_CKDtomort)
#sd(tp_CKDtomort)

#CKD to ESRD on maintanence dialysis 
#Base case: from Kent et al. 2015. Number of patient years with events: For CKD 1-3: 77 (1.3%). For CKD 4=548 (6.2%)  
tp_CKDtodial <<-   (1494/(1494+2228))*rbeta(Nsim,77, (77/0.013)-77) + (2228/(1494+2228))*rbeta(Nsim,548, (548/0.062)-548)    
#hist(tp_CKDtodial)
# summary(tp_CKDtodial)
# sd(tp_CKDtodial)

#CKD to ESRD no dialysis      #check assumption here
#From kent we have 1,017 ESRD patients not on dialysis, vs. 2498 on dialysis. Assume this proportion remains the same over time. 
#Then, using tp_CKDtodial we can make an assumption about the additional number of patients who would transition to ESRD no dialysis. 
tp_CKDtoESRD <<- tp_CKDtodial*rbeta(Nsim,1017,2498)
#1494/2498
# hist(tp_CKDtoESRD)
# summary(tp_CKDtoESRD)
# sd(tp_CKDtoESRD)


#ESRD no dialysis to ESRD+dialysis      
#Base case Kent et al 2015. Total number of patient years with event= 718= 18.2%. Total number patient years= 718/0.182=3945. 
tp_ESRDtodial <<- rbeta(Nsim,718,3945-718)
#hist(tp_ESRDtodial)


#Transitions from ESRD +dialysis
#Base case: 18th Annual Renal Registry Report 
temp <- rdirichlet(Nsim, c(3265,258,641))
tp_dialysis     <<- temp[,1]
tp_dialtotrans  <<- temp[,2]
tp_dialtomort   <<- temp[,3]

#ESRD no dialysis mortality            
#Base case: Kent et al 2015. Total number of patient years with event= 86+151= 237= 6%. Total patient years= 237/0.06=3950. 
tp_ESRDtomort <<- tp_dialtomort*(6/7.5)
#tp_ESRDtomort <<- rbeta(Nsim,237,3950-273)
#hist(tp_ESRDtomort)

#ESRD no dialysis to Transplant 
#Base case: Kent et al 2015. 335 patient years with event = 8.5%. Total patient years=335/0.085= 3941.     
tp_ESRDtotrans <<-  rbeta(Nsim,335,3941-335)
#hist(tp_ESRDtotrans)
  
#Transitions from Transplant 
#Base case: 18th Annual Renal Registry Report 
temp <- rdirichlet(Nsim, c(329,5,2))
tp_transplant   <<- temp[,1]
tp_transtodial  <<- temp[,2]
tp_transtomort  <<- temp[,3]



#=============================== UTILITIES ========================================# 
  
#=====Post-discharge recovery utility====#
# Base case: Cuthbertson et al 2011.           
# Utility up to year 1 post discharge= 0.67, sd= 0.28; up to year 4 = 0.70, sd=0.281; year 5 onwards=0.68. sd=0.301. 

# YEAR 1
# Cannot fit beta distribution here, as sd is too large and function collapses. 
# Set to disutility from 1. mean=1-0.666, sd=0.28
u_FU_1 <<- 1-gamma.draw(Nsim, 1-0.66, 0.28)   #u_FU_1_dec
#hist(u_FU_1, breaks=200)        

# YEAR 2-4
u_FU_2 <<- 1-gamma.draw(Nsim, 1-0.701,0.281)
#hist(u_FU_2, breaks=100)                    

# YEAR 5 onwards
u_FU_5  <<- 1-gamma.draw(Nsim, 1-0.677, 0.301)
#hist(u_FU_5, breaks=100)                    

#_______________________________________________

#Utility successful transplant
# Base case: assumed equivalent to 'recovered' value
u_transplant     <<- u_FU_5
#hist(u_transplant, breaks=100)    

#Utility for CKD 1-4
# Base case: decrement from successful transplant taken from Wyld et al. 2012
u_CKD_decrement  <- gamma.draw(Nsim, 0.02 ,0.03)
#u_CKD_decrement  <- rnorm(Nsim, 0.02 ,0.03)
u_CKD <<- u_transplant - u_CKD_decrement  
#hist(u_CKD_decrement, breaks=1000)         

#Utility ESRD no dialysis 
# Base case: decrement from successful transplant taken from Wyld et al. 2012
u_ESRD_decrement  <- gamma.draw(Nsim, 0.2,0.09)
u_ESRD <<- u_transplant - u_ESRD_decrement  
#hist(u_ESRD, breaks=100)

#Utility ESRD on maintenance dialysis 
# Base case: decrement from successful transplant taken from Wyld et al. 2012
u_dialysis_decrement  <- rnorm(Nsim, 0.11,0.02)  
u_dialysis <<- u_transplant - u_dialysis_decrement  
# hist(u_dialysis, breaks=100)


  
#=================================== COSTS ========================================#
  
#NB: All costs in 2015 prices. 

#Follow-up cost post ICU stay    
# Base case: Lone et al. 2016. Reported 1-year post discharge cost mean 2014 US$=8863 (8332-9393). N=6317. 
#Converted to 2015 GBP using EPPI converter: mean 6230.46 (5857.18 to 6603.03). 
c_fup_yr1  <<- lnorm.draw(Nsim, 6230.46, (6603.03-5857.18)/(3.92))  
#hist(c_fup_yr1, breaks=100)  

#Year 2. Mean 5704 (5269 to 6138). N=4116. 
#Converted to GBP using EPPI converter: mean 4009.76 (3703.97 to 4314.85). 
c_fup_yr2  <<- lnorm.draw(Nsim, 4009.76, (4314.85-3703.97)/(3.92))  
#hist(c_fup_yr2, breaks=100)

#Year 3. Mean 5421 (4949 to 5892). N=4104.
#Converted to GBP using EPPI converter: mean 3810.82 (3479.02 to 4141.92). 
c_fup_yr3  <<- lnorm.draw(Nsim, 3810.82, (4141.92-3479.02)/(3.92))
#hist(c_fup_yr3, breaks=100)

#Year 4. Mean 5146 (4639 to 5652). N=4078.
#Converted to GBP using EPPI converter: mean 3617.50 (3261.10 to 3973.21).
c_fup_yr4  <<- lnorm.draw(Nsim, 3617.50, (3973.21-3261.10)/(3.92))
#hist(c_fup_yr4, breaks=100)

#Year 5. Mean $4521 (4061 to 4982). N=3593. 
#Converted to GBP using EPPI converter: mean 3178.14 (2854.78 to 3502.21).
#(3502.21-2854.78)/(3.92)    #sd= 165.1607
c_fup_yr5  <<- lnorm.draw(Nsim, 3178.14, (3502.21-2854.78)/(3.92) )
#hist(c_fup_yr5, breaks=100)

#Year 6-10 and year 11 onwards
#Assume same downward trend and variance as observed in last year reported in Lone et al. (165.1607) 
# 3617.50-3178.14 #=439.36 drop in costs
c_fup_yr6 <<-lnorm.draw(Nsim, 3178.14-439.36, 165.1607)
c_fup_yr7 <<-lnorm.draw(Nsim, 3178.14-2*439.36, 165.1607)
c_fup_yr8 <<-lnorm.draw(Nsim, 3178.14-3*439.36, 165.1607)
c_fup_yr9 <<-lnorm.draw(Nsim, 3178.14-4*439.36, 165.1607)
c_fup_yr10 <<-lnorm.draw(Nsim,3178.14-5*439.36, 165.1607)
c_fup_yr11 <<-lnorm.draw(Nsim,3178.14-6*439.36, 165.1607)       #Use for year 11 and onwards
#hist(c_fup_yr11)

c_fup<<- data.frame(c_fup_yr1, c_fup_yr2,c_fup_yr3,c_fup_yr4,c_fup_yr5,c_fup_yr6,c_fup_yr7,c_fup_yr8,c_fup_yr9,c_fup_yr10)


#Impact of AKI in ICU on follow-up costs
#Base case: Lone et al. 2016. Reported Admission rate ratio of 1.15 (95% CI: 1.02 to 1.31) for use of RRT in ICU stay on 5-year follow-up hospital admissions. 
#Use as proxy for impact on costs. 
#(1.31-1.02)/3.92
ratio_fupcost <- rnorm(Nsim, 1.15, 0.07397959)
#hist(ratio_fupcost)

#Apply ratio for first 5 years of follow up costs. 
c_fupAKI <<- data.frame(c_fup_yr1*ratio_fupcost, c_fup_yr2*ratio_fupcost,c_fup_yr3*ratio_fupcost,c_fup_yr4*ratio_fupcost,c_fup_yr5*ratio_fupcost,c_fup_yr6,c_fup_yr7,c_fup_yr8,c_fup_yr9,c_fup_yr10)



#Additional Cost of CKD                 
# Nb: assume that c_fup captures additional costs due to stay in ICU; cost of CKD, ESRD etc. (derived from papers not looking specifically at post-ICU populations) is therefore added on in addition to the fup costs in the model (see model cost code)
# Base case: Kent 2015. 
# Costs based on reported baseline cost of CKD, plus additional costs of pre-existing diabetes and cardiac diseases, plus cost of non-fatal MVE
# Mean=£578.74, 95% CI: 470.4 - 689.55. 
c_CKD <<-  lnorm.draw(Nsim, 578.74,(689.55-470.4 )/(3.92))             
#hist(c_CKD, breaks=200)


#Cost of ESRD no dialysis
# Base case: Kent 2015. 
# Mean=£759.55, 95% CI: 625.42 - 897. 
c_ESRD <<-  lnorm.draw(Nsim, 759.55,(897-625.42)/(3.92)) 
#hist(c_ESRD, breaks=100)

#Cost of ESRD + dialysis in year 1         
# Base case: Kent 2015.
# Mean=£20440.01, 95% CI: 19904.07 - 20820.92.
c_dialysis_yr1 <<-   lnorm.draw(Nsim, 20440.01,(20820.92-19904.07 )/(3.92))  
#hist(c_dialysis_yr1)

#Cost of ESRD + dialysis in year 2+
# Base case: Kent 2015. 
# Mean=£25035.34, 95% CI: 24786.34 - 25129.31
c_dialysis <<- lnorm.draw(Nsim, 25035.34,(25129.31-24786.34)/(3.92))  
#hist(c_dialysis, breaks=100)

#Cost of successful transplant in year 1    
# Base case: Kent 2015. 
# Mean=£26301.22, 95% CI: 25578.29 - 26915.51
c_transplant_yr1 <<-  lnorm.draw(Nsim, 26301.22,(26915.51-25578.29)/(3.92))  
#hist(c_transplant_yr1)

#Cost of successful transplant year 2+
# Base case: Kent 2015. Reported costs (£1,148 [978 - 1,318]; 2011 prices) inflated to 2013/14 prices using EPPI-Centre Cost Converter (http://eppi.ioe.ac.uk/costconversion/default.aspx) 
# Mean=£1467.39, 95% CI: 1173.29- 1651.8
c_transplant <<- lnorm.draw(Nsim, 1467.39,(1651.8-1173.29)/(3.92)  ) 
#hist(c_transplant)


}
