########################################## AKI DIAGNOSTICS MODEL: HOSPITAL PERIOD MODEL PARAMETERS  ###############################################################

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

parameter_draw_HOSP  <- function() { 
  
  set.seed(seed)
  
# Global parameters of model 
  S  <- 14            # Number of health states in model 
  T  <- 90            # Time Horizon = 90 days (3 months; cycle lengths of 1 day)

  
#========================================================================================================================================#
#====================================================== BASELINE PARAMETERS =============================================================#
#========================================================================================================================================#
  

#====================== CODE FOR FUNCTIONS USED ========================# 

#Code to generate beta/gamma/lnorm distribution parameters based on aggregate reported data (e.g. number [n], mean [m], standard error or deviation [s])  
  
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
# Base case: from Leeds AKI Registry data. Mean age of allcomers to ICU= 61.3351, sd=15.37216.
startage <<- 61.3351  

#Proportion of cohort male 
# Base case: from Leeds AKI Registry data. 263 out of 376 allcomers to ICU were male (113 female)
# Used within background mort calculation only (within excel file), not an explicit parameter.  



#============================= PROBABILITIES ======================================# 

#Proportion of allcomers to ICU who arrive or develop AKI in ICU
# Base case: Susantitaphong 2013. Pooled rate in Critical Care setting. Mean= 0.317. 95% CI: 0.28 to 0.35. N=888,604. 
# Assume normality to derive variance from mean and 95% CI. Large N, so assume standard error= standard deviation. 
# Sensitivity analysis (Cardiac Population): Susantitaphong 2013. Mean= 0.243. 95% CI: 0.204 to 0.288. N=164,333.  
p_AKI  <<- if(Model.pop=="Normal.pop") beta.draw(Nsim, 0.317, (0.35-0.28)/(2*1.96)) else (if(Model.pop=="CVsurg.pop") beta.draw(Nsim, 0.243, (0.288-0.204)/(2*1.96)) else (NA))
#hist(p_AKI, breaks=100)

#Health state starting distributions for AKI cohort, baseline 
# Base case: from analysis of Leeds AKI Registry data (AKI ICU cohort)
# Health states: c(AKI S0 (currently normal kidney function but destined to get AKI), S1, S2, S3, S3+RRT)
day_0  <- c(18,25,8,7,2)  #raw data
p_day0  <<- data.frame(rdirichlet(Nsim, day_0))
start.noAKI        <<- 1-p_AKI                 
start.S0           <<- p_AKI*p_day0[,1]    
start.S1           <<- p_AKI*p_day0[,2]    
start.S2           <<- p_AKI*p_day0[,3]    
start.S3           <<- p_AKI*p_day0[,4]    
start.RRT          <<- p_AKI*p_day0[,5]    
# Validity check: starting proportions should add to 1     
# Nb: For testing arms, patients with existing AKI may not be tested, and therefore need to be separated (see testing parameters at end)



#============================== TRANSITION PROBABILITIES ==============================#

#============== AKI cohort in-hospital transitions ==============#                      

#All AKI cohort transition probabilities are derived from AKI Registry daily data for cohort of 60 Leeds patients who had or developed AKI in ICU in 7 days from admission. 
#All transition probability (TP) matrices have dimension [Nsim, AKI cohort health states, T] i.e. [10000, 10, 90]

#Note all transition probability array columns are labbeled according to names above e.g. "mean(tp_norm[,"mort",])" gives mean probability of moving from S0 to death over time in the model
#Nb column numbers in these tp matrices correspond to the health states as follows:
  #                             #1=  "S0" - In ICU with normal kidney function (AKI cohort)
  #                             #2=  "S1" - In ICU with AKI Severity S1 (KDIGO grading)
  #                             #3=  "S2" - In ICU with AKI Severity S2
  #                             #4=  "S3" - In ICU with AKI Severity S3
  #                             #5=  "RRT" - In ICU with AKI Severity S3 + RRT 
  #                             #6=  "ward" - Transferred to General Ward 
  #                             #7=  "ward_RRT" - Transferred to General Ward +RRT
  #                             #8=  "mort" - Dead
  #                             #9=  "disch" - Discharged home 
  #                             #10=  "disch_RRT" Discharged home +RRT

#Load transition probability .RData files
#NEW USERS - Will need to update the file path names here

load("Model Transition Probabilities\\TPs for model- 10000 simulations 0.01 added\\tp_norm.RData")  
tp_norm <<- tp_norm
# plot(colMeans(tp_norm[,"S0",1:90]))
#sum(tp_norm[,1:10,1])

load("Model Transition Probabilities\\TPs for model- 10000 simulations 0.01 added\\tp_S1.RData")
tp_S1 <<- tp_S1

load("Model Transition Probabilities\\TPs for model- 10000 simulations 0.01 added\\tp_S2.RData")
tp_S2 <<- tp_S2

load("Model Transition Probabilities\\TPs for model- 10000 simulations 0.01 added\\tp_S3.RData")
tp_S3 <<- tp_S3

load("Model Transition Probabilities\\TPs for model- 10000 simulations 0.01 added\\tp_RRT.RData")
tp_RRT <<- tp_RRT

load("Model Transition Probabilities\\TPs for model- 10000 simulations 0.01 added\\tp_ward.RData")
tp_ward <<- tp_ward

load("Model Transition Probabilities\\TPs for model- 10000 simulations 0.01 added\\tp_ward_RRT.RData")
tp_ward_RRT <<- tp_ward_RRT

##Mean death rates in each health state
#hist(tp_norm[,"mort",], breaks=100, xlim=c(0,0.2))
#mean(tp_norm[,"mort",])

#NOTE: The AKI Registry does not give data on mortality or RRT post discharge.     
#See further down, after no AKI cohort values, where this is calculated based on the Lone et al. paper. 



#================== No AKI cohort transitions ==================#

#====== ICU mortality ======#

# Load AKI Registry ICU survival Kaplan Meier data 
# 10000 simulations of KM data- sample with replacement (nb: code for this in "Model Transition Probabilities\\AKI-Registry_survival analysis.R")
load("Model Transition Probabilities\\tp_ICUmort.RData")   
tp_ICUmort <<- tp_ICUmort

#RR for ICU and hospital mortality derived by model calibration.  
# Alter RR by trial and error to get relative mortality probabilities from model (AKI vs. no AKI) to match literature values. 
# Base case: OR from Susantiphong 2013 (OR=3.93), and mean of ratio observed in 3 papers (Abosaif 2005, De-Mendonca 2000 & Ostermann 2007) (mean RR= 3.92).  
# Assume that at day 10 mortality in AKI cohort should be ~3.93 times that of no AKI cohort. Get RR=0.30.Assume sd=rr/3.            
RR_mort <<- if (Model.pop=="Normal.pop") lnorm.draw(Nsim,0.30,0.30/3) else (if (Model.pop=="CVsurg.pop") lnorm.draw(Nsim,0.12, 0.12/3) else (NA))      
#hist(RR_mort, breaks=100)

#apply RR to ICU mortality
for (i in 1:Nsim){
 tp_ICUmort[i,] <- tp_ICUmort[i,3]*RR_mort[i]    
}
tp_ICUmort <<-tp_ICUmort
#summary(tp_ICUmort[,1])



#====== ICU discharge ======# 

# Apply RR to probability of discharge observed in AKI Registry, to get a 'No AKI' cohort probability of discharge. 
# Load AKI Registry ICU length of stay Kaplan Meier data (10000 simulations of KM data- sample with replacement)
load("Model Transition Probabilities\\tp_ICUdisch_alive.RData")
tp_ICUdisch_alive <- tp_ICUdisch_alive

# Apply RR between ICU LOS for AKI vs. no AKI cohorts observed across 3 UK studies: De-Mendonca 2000, Prowle 2014, Osterman 2007.
# Sensitivity analysis: Bastin et al. 2013
# RR=0.542 sd 0.151 i.e. (using simple mean)
RR_ICUdisch <<- if (Model.pop=="Normal.pop") lnorm.draw(Nsim,0.542,0.151) else (if (Model.pop=="CVsurg.pop") lnorm.draw(Nsim,0.216,0.216/3) else (NA))      
#hist(RR_ICUdisch, breaks=100)

#Apply RR to transitions from the KM data (note dividing here)
for (i in 1:Nsim){
tp_ICUdisch_alive[i,] <-tp_ICUdisch_alive[i,]/RR_ICUdisch[i]
}

#Adjust for mortality (have increased discharge rates, so may push discharged+mort probabilities to >1)
for (i in 1:Nsim){
  tp_ICUdisch_alive[i,] <-ifelse(tp_ICUdisch_alive[i,]+tp_ICUmort[i,]<=1, tp_ICUdisch_alive[i,], tp_ICUdisch_alive[i,]-((tp_ICUdisch_alive[i,]+tp_ICUmort[i,])-1))
}

#Probability discharged from ICU to ward vs. home                      
# Base case: from AKI Registry 29/60 patients discharged home on same day as discharged from ICU. 
p_disch <<- rbeta(Nsim, 29, 60-29)  

#Calculate ICU transition probabilities for No AKI cohort
#Nb: have to take away mortality proportion first as discharge data inlcuded discharge due to mortality. 
tp_ICUtoward <<- tp_ICUdisch_alive*(1-p_disch)
tp_ICUtodisch <<- tp_ICUdisch_alive*(p_disch)


#====== Hospital mortality ======#
#From AKI Registry, 3 out of 29 patients died in hospital, post ICU. Median LOS for those patients was 11 days.  
#Calculate daily mortality rate  #1-exp(--(log(1-(3/29))/11)*1) =0.009878096
tp_Hospmort <- beta.draw(Nsim, 0.009878096, 0.009878096/3)
#hist(tp_Hospmort, breaks=100)

#Apply RR to get No AKI cohort probability
tp_wardmort <- rep(NA, times=Nsim)
for (i in 1:Nsim){
tp_wardmort[i]<-tp_Hospmort[i]*RR_mort[i]
}
tp_wardmort<<-tp_wardmort
#hist(tp_wardmort,breaks=100)


#====== Hospital discharge ======#

#Ratio AKI/noAKI Hospital LOS
# Base case: Prowle et al. RR=1.447.
# Sensitivity analysis (post-cardiac population): Bastin et al. 2013. RR=1.914. 
# Median Hospital LOS from AKI Registry: 11 days 

#Calculate daily probability and define draw
m  <- (1-(exp(-(-(log(1-0.5)/(11)*1)))))*(if(Model.pop=="Normal.pop") 1.447 else(if(Model.pop=="CVsurg.pop") 1.914 else (NA)))       
tp_Hospdisch <- beta.draw(Nsim, m, m/3)
#hist(tp_Hospdisch, breaks=100)

#Get ward discharge (LOS data included deaths, so need to adjust for this)
for (i in 1:Nsim){
tp_wardmort[i]<- min(tp_wardmort[i], tp_Hospdisch[i])
}

tp_wardtodisch <<- tp_Hospdisch- tp_wardmort
tp_wardtodisch[tp_wardtodisch<0]<-0
tp_wardtodisch[tp_wardtodisch>1]<-1
#summary(tp_wardtodisch)
#summary(tp_Hospdisch)
#summary(tp_wardmort)


#====== Post discharge mortality ======#

#Base case: from Lone et al 2016. Reported 1-year mortality for 5,259 ICU cohort (post hospital discharge)= 10.9% (95% CI 10.0 to 11.7) 
# Convert yearly probability and 95% CIs to daily probabilities: 
# 1-exp(-(-log(1-0.109)/365)*1) #= 0.0003161441
# ((1-exp(-(-log(1-0.117)/365)*1)) - (1-exp(-(-log(1-0.10)/365)*1)))/3.92 #= 1.33237e-05
tp_dischmort <<- beta.draw(Nsim,0.0003161441,  1.33237e-05)
#hist(tp_dischmort, breaks=100)

#Discharged + RRT mortality
#Base case: derive from FUP model value for ESRD+dialysis mortality from 18th Annual Renal Registry Report 
temp <- rdirichlet(Nsim, c(3265,258,641))
tp_dialtomort   <- temp[,3]
# Convert yearly probability to daily probability
tp_dischmort_RRT <<- 1-exp(-(-log(1-tp_dialtomort)/365)*1)


#=============================== UTILITIES ========================================# 

#ICU utility 
# Base case: Kind et al, 1999. AKI utility assumed equal to that of an unconcious person. Mean= -0.402, sd= 0.201
# Mean =0 and =0.20 used in sensitivity analyses (assuming same sd as base case). 
u_ICU <<- rnorm(Nsim, -0.402, 0.201)    
#summary(u_ICU)
#hist(u_ICU)

#Non-ICU (general) hospital ward (dialysis independent)
# Base case: Hernandez et al, 2013. Baseline value post-ICU discharge. Mean= 0.44, sd=0.31. 
# Cannot fit beta distribution here, as sd is too large and function collapses. Also, we want to allow negative utilities. 
u_ward  <<- 1-gamma.draw(Nsim, 1-0.44, 0.31)   
#hist(u_ward)                           
                                                                           
#Discharged utility (up to 90 days) (dialysis independent)    
# Base case: Hernandez et al, 2013. Value for 6-month post-ICU discharge. Mean=0.62, sd=0.32. 
u_disch  <<- 1-gamma.draw(Nsim, 1-0.38, 0.32)  
#hist(u_disch)

#Decrement utility for dialysis dependence (any time point post ICU discharge)
# Base case: decrement from successful transplant taken from Wyld et al. 2012. Mean=0.11, sd=0.02. 
u_dialysis_dec  <- rnorm(Nsim, 0.11, 0.02)     
#hist(u_dialysis_dec)

#Non-ICU (general) hosital ward + on dialysis
u_ward_dialysis <<- u_ward - u_dialysis_dec 
#hist(u_ward_dialysis)

#Discharged utility (up to 90 days) + on dialysis
u_disch_dialysis <<- u_disch - u_dialysis_dec 
#hist(u_disch_dialysis)


#=================================== COSTS ========================================#

#NB: All costs in 2014/15 prices. 

#ICU normal kidney function   
# Base case: NHS Reference costs 2014/15. 
# Weighted mean of adult general critical care= £1,306.16, LQ=1019.53, UQ=1,508.72, Total N= 1012404. 
# Sensitivity analysis: weighted mean of post cardiac surgery critical care= £1,274.92, LQ=910.31, UQ=1,594.79. 
c_ICU  <<- if(Model.pop=="Normal.pop") lnorm.draw(Nsim, 1306.16, (1508.72- 1119.53)/(2*0.67)) else (if(Model.pop=="CVsurg.pop") lnorm.draw(Nsim, 1274.92, (1594.79- 910.31)/(2*0.67)) else (NA))         
#hist(c_ICU, breaks=100)

#Non-ICU hospital ward (dialysis independent)
# Base case: PSSRU 2015  
# Half cost of non-elective inpatient short stay (assume average short stay= 2 days). Mean= 304, LQ= 206, UQ=355. sd =  111.194 
c_ward  <<- lnorm.draw(Nsim, 304, 111.194)              
#hist(c_ward, breaks=100)

#Excess cost of AKI in hospital 
# Base case: NHS Reference costs 2014/15. 
# Weighted average of costs for AKI with Interventions across all cc bands. Mean= 265, LQ= 205, UQ=309, Total N=39883. 
c_AKI  <<- lnorm.draw(Nsim, 265, (308.55-205.24)/(2*0.67))            
#hist(c_AKI,breaks=100)

#Excess cost of dialysis in hospital 
# Base case: NHS Reference costs 2014/15. Weighted average of costs for Haemodialysis for AKI. Mean= 690.51, LQ= 176.18, UQ=1721.90. 
c_dialysis  <<- lnorm.draw(Nsim, 690.51, (1721.90-176.18)/(2*0.67))   
#hist(c_dialysis, breaks=100)

#ICU with AKI (S1-S3) 
c_ICU_AKI <<- c_ICU + c_AKI

#ICU on dialysis 
c_ICU_dialysis <<- c_ICU + c_AKI + c_dialysis

#Ward on dialysis
c_ward_dialysis <<- c_ward + c_dialysis

#Discharged from hospital (dialysis independent)
# Base case: Lone et al 2016. Reported hospital 1-year post discharge costs. (see fup costs c_fup_yr1)
c_disch  <<-lnorm.draw(Nsim, 6230.46, 190.2679)/365 
#hist(c_disch, breaks=100)  

#Discharged from hospital on dialysis  
# Base case: Kent et al, 2015. Assumed daily cost =annual cost of chronic dialysis/365 (see FUP parameters for working)
c_RRT <<- lnorm.draw(Nsim, 25035.34,87.49235)/365  
c_disch_dialysis <<- c_disch + c_RRT
#hist(c_disch_dialysis, breaks=100)
#hist(c_RRT)





#========================================================================================================================================#
#==================================================== INTERVENTION PARAMETERS ===========================================================#
#========================================================================================================================================#

#Impact of tests summary: 
              # Additional cost of testing applied to all patients
              # Additional treatment (i.e. consultation) cost added to all 'positive' tests
              # TP proportion in AKI cohort get reduced risks of AKI & mortality 


#========================= Test Costs ==========================#

##=== Nephrocheck ===#

#Astute Medical 140 Meter test costs (including cost of platform)
c_Nephro <<- rep(71.27, times=Nsim)
#hist(c_Nephro, breaks=100)


##=== NGAL ===#

#Two tests for NGAL: urine and plasma.
#Two current platforms for NGAL: BioPorto (can use plasma or urine) and Abbott (urine only). Assume BioPorto in base case. 

#BioPorto test costs
# Cost assumed not to be dependant on urine vs plasma test. 
# Per test direct cost from communications with BioPorto company, using Siemend ADVIA 1800 platform: 100 tests kit= £1,400 = £14/ test; 300 tests kit cost= £3,600 = £12/ test. 
# Plus data on laboratory overheads, staff costs and shipping costs: £17.55 (100 tests kit) or £14.98 (300 tests kit). 
#Base case: assume NHS would maximise efficiency and use cheapest cost. 
c_NGAL  <<- rep(14.98,  times=Nsim)

#Sensitivity analysis 
#Abbott test costs
#Kit cost £2000 for 100 kits, £20 per test. 
#Plust Cal, QC, shipping, staff and overhead costs = £28.14. 
#c_NGAL  <- rep(28.14, times=Nsim)


##=== Cystatin C ===#

# Communicated per-test cost from company, using Siemens ADVIA XPT platform                   
# Kit cost= £250.63 for 400 kits = £0.63/test. 
# Plust Cal, QC, staff & overhead costs = £4.00
c_CystC <<- rep(4.26, times=Nsim)  



#================== Early AKI Intervention Costs ===============#

#Base case: NHS Reference cost 2014/15. 
#Cost for critical care intensivist consultation. WF01B. Non-Admitted Face to Face Attendance, First. Critical Care Medicine. 
# Unit cost=£250.19, Lower quartile: 29.71; upper: 403.52 
c_early.treat   <<- lnorm.draw(Nsim,250.19,(403.52 -29.71)/1.34)             
#hist(c_early.treat, breaks=100)



#============= IMPACT OF EARLY INTERVENTION ON RISKS ================#        

#Base case
#From review of early nephrology involvement: Flores Gama et al, 2013. 
# Reported OR for developing AKI of 0.71 (25th & 75th percentiles: 0.53-0.95). 25.7% vs. 31.9%. 
# Use assumed controlled risk of 0.319 (the value reported for baseline arm) to convert OR to RR (see http://handbook.cochrane.org/chapter_12/12_5_4_4_computing_risk_ratio_from_an_odds_ratio_or.htm) 
#0.71/(1-0.319*(1-0.71)) #Convert OR to RR. Mean=0.782377767  Upper%=0.9653981, Lower%= 0.623478067. sd=0.2534618
RR_early.treat <<- lnorm.draw(Nsim,0.782377767,0.2534618)   #NB: a lot of values >1. 
#hist(RR_early.treat, breaks=100)

#Sensitivity analysis: reset values >1 to 1. Limitation: will distort mean value. 
# for (i in 1:Nsim){
# RR_early.treat[i] <- ifelse(RR_early.treat[i]>1, 1, RR_early.treat[i])
# }
# RR_early.treat<<- RR_early.treat

#Sensitivity analysis: truncate the variance to pull tail below 1. Limitation: Will make results look more certain and influence any value of information analysis. 
# RR_early.treat <<- beta.draw(Nsim,0.782377767,0.06) 



#============================== Test Accuracy ================================#       

#Base case: use outputs from meta-analysis (mean sensitivity and specificity, and variance-covariance matrix)


#========== Nephrocheck ==========#  

#===Base case===# 
mu <- c(0.89897078, 0.491559); mu <- array(mu, c(2,1))
Vcov  <- c(4.812714e-04, 4.238991e-05, 4.238991e-05, 3.116314e-04); Vcov <- array(Vcov, c(2,2))
temp <- mvrnorm(Nsim, mu, Vcov)

#=== Sensitivity analysis - post cardiac surgary ===#
#One paper identified from review: Meersch et al 2014
#Reported accuracy at 4 hours post surgery: sensitivity = 0.80; specificity= 0.83.
#No variance data- assume same variance-covariance matrix is the same as from the meta-analysis.
mu <- c(0.8, 0.83); mu <- array(mu, c(2,1))
Vcov  <- c(4.812714e-04, 4.238991e-05, 4.238991e-05, 3.116314e-04); Vcov <- array(Vcov, c(2,2))
temp_cardio <- mvrnorm(Nsim, mu, Vcov)

Sens_Nephro <<- if(Model.pop=="Normal.pop") temp[,1] else (if (Model.pop=="CVsurg.pop") temp_cardio[,1] else (NA))
Spec_Nephro <<- if(Model.pop=="Normal.pop") temp[,2] else (if (Model.pop=="CVsurg.pop") temp_cardio[,2] else (NA))

# hist(Sens_Nephro, breaks=100)
# hist(Spec_Nephro, breaks=100)
#length(Sens_Nephro[Sens_Nephro>1]) #0
#length(Spec_Nephro[Spec_Nephro>1]) #0


#========= NGAL-plasma =========#

#===Base case===#
mu <- c(0.72, 0.81); mu <- array(mu, c(2,1))
Vcov  <- c(0.001388834,-0.000305304 , -0.000305304, 0.000727345); Vcov <- array(Vcov, c(2,2))
temp <- mvrnorm(Nsim, mu, Vcov)

#=== Sensitivity analysis - post cardiac surgary ===#
mu <- c(0.61, 0.77); mu <- array(mu, c(2,1))
Vcov  <- c(0.00433586592537737,0.000673262, 0.000673262, 0.000221851); Vcov <- array(Vcov, c(2,2))
temp_cardio <- mvrnorm(Nsim, mu, Vcov)

Sens_NGALp <<- if(Model.pop=="Normal.pop") temp[,1] else (if (Model.pop=="CVsurg.pop") temp_cardio[,1] else (NA))
Spec_NGALp <<- if(Model.pop=="Normal.pop") temp[,2] else (if (Model.pop=="CVsurg.pop") temp_cardio[,2] else (NA))

# hist(Sens_NGALp, breaks=100)
# hist(Spec_NGALp, breaks=100)
#length(Sens_NGALp[Sens_NGALp>1]) #0
#length(Spec_NGALp[Spec_NGALp>1]) #0


#========= NGAL-urine =========#

#===Base case===#
mu <- c(0.7, 0.79); mu <- array(mu, c(2,1))
Vcov  <- c(0.002771809,0.00136847, 0.00136847, 0.001424614); Vcov <- array(Vcov, c(2,2))
temp <- mvrnorm(Nsim, mu, Vcov)

#=== Sensitivity analysis - post cardiac surgary ===#
mu <- c(0.66, 0.62); mu <- array(mu, c(2,1))
Vcov  <- c(0.003412758,0.000648401 ,0.000648401 , 0.009868229); Vcov <- array(Vcov, c(2,2))
temp_cardio <- mvrnorm(Nsim, mu, Vcov)

Sens_NGALu <<- if(Model.pop=="Normal.pop") temp[,1] else (if (Model.pop=="CVsurg.pop") temp_cardio[,1] else (NA))
Spec_NGALu <<- if(Model.pop=="Normal.pop") temp[,2] else (if (Model.pop=="CVsurg.pop") temp_cardio[,2] else (NA))

# hist(Sens_NGALu, breaks=100)
# hist(Spec_NGALu, breaks=100)
#length(Sens_NGALu[Sens_NGALu>1]) #0
#length(Spec_NGALu[Spec_NGALu>1]) #0

#========= NGAL-serum =========#

#===Base case===#
#Only one study, Chen 2012
#Use reported sensitivity and specificity, and assume same variance co-variance as for cardiac
mu <- c(0.92, 0.69); mu <- array(mu, c(2,1))
Vcov  <- c(0.017359294,-0.010533545 ,-0.010533545 , 0.007687521); Vcov <- array(Vcov, c(2,2))
temp <- mvrnorm(Nsim, mu, Vcov)

#=== Sensitivity analysis - post cardiac surgary ===#
mu <- c(0.84, 0.87); mu <- array(mu, c(2,1))
Vcov  <- c(0.017359294,-0.010533545 ,-0.010533545 , 0.007687521); Vcov <- array(Vcov, c(2,2))
temp_cardio <- mvrnorm(Nsim, mu, Vcov)

Sens_NGALs <<- if(Model.pop=="Normal.pop") temp[,1] else (if (Model.pop=="CVsurg.pop") temp_cardio[,1] else (NA))
Spec_NGALs <<- if(Model.pop=="Normal.pop") temp[,2] else (if (Model.pop=="CVsurg.pop") temp_cardio[,2] else (NA))

for (i in 1:Nsim){
  Sens_NGALs[i] <- ifelse(Sens_NGALs[i]>1, 1, Sens_NGALs[i])
  Spec_NGALs[i] <- ifelse(Spec_NGALs[i]>1, 1, Spec_NGALs[i])
}
Sens_NGALs <<- Sens_NGALs 
Spec_NGALs <<- Spec_NGALs

#hist(Sens_NGALs, breaks=100)   #lots>100
#hist(Spec_NGALs, breaks=100)   #a few >100
#length(Sens_NGALs[Sens_NGALs>1])   # Base case= 2635;   cardio=1135
#length(Sens_NGALs[Spec_NGALs>1])   # Base case= 2;      cardio=710


# #========= Cystatin C =========#

#test =31, Cystatin C plasma
#32 = urine
#33=serum

#===Base case plasma===# 
mu <- c(0.72, 0.74); mu <- array(mu, c(2,1))
Vcov  <- c(0.003488183,0.001272777 ,0.001272777 , 0.001716697); Vcov <- array(Vcov, c(2,2))
temp31 <- mvrnorm(Nsim, mu, Vcov)

#===Base case urine===# 
mu <- c(0.68, 0.76); mu <- array(mu, c(2,1))
Vcov  <- c(0.013701728,0.001430536 ,0.001430536, 0.003842759); Vcov <- array(Vcov, c(2,2))
temp32 <- mvrnorm(Nsim, mu, Vcov)

#===Base case serum===# 
mu <- c(0.76, 0.88); mu <- array(mu, c(2,1))
Vcov  <- c(0.006605824,0.00109115,0.00109115, 0.000620969); Vcov <- array(Vcov, c(2,2))
temp33 <- mvrnorm(Nsim, mu, Vcov)

#=== Sensitivity analysis - post cardiac surgary ===#
#Nb: no pooled analysis available for this for plasma test 

#only one paper for CystC plasma (Tziakas 2015). Use reported accuracy and assume var-covar matrix same as base case
mu <- c(0.61,0.56); mu <- array(mu, c(2,1))
Vcov  <- c(0.003488183,0.001272777 ,0.001272777 , 0.001716697); Vcov <- array(Vcov, c(2,2))
temp_c31 <- mvrnorm(Nsim, mu, Vcov)

mu <- c(0.52,0.72); mu <- array(mu, c(2,1))
Vcov  <- c(0.019440913,0.020658555 ,0.020658555 , 0.024791684); Vcov <- array(Vcov, c(2,2))
temp_c32 <- mvrnorm(Nsim, mu, Vcov)

mu <- c(0.73,0.72); mu <- array(mu, c(2,1))
Vcov  <- c(0.001390272,-4.45616E-05 ,-4.45616E-05 , 0.001754969); Vcov <- array(Vcov, c(2,2))
temp_c33 <- mvrnorm(Nsim, mu, Vcov)


Sens_CystCp <<- if(Model.pop=="Normal.pop") temp31[,1] else (if (Model.pop=="CVsurg.pop") temp_c31[,1] else NA) 
Spec_CystCp <<- if(Model.pop=="Normal.pop") temp31[,2] else (if (Model.pop=="CVsurg.pop") temp_c31[,2] else NA) 

Sens_CystCu <<- if(Model.pop=="Normal.pop") temp32[,1] else (if (Model.pop=="CVsurg.pop") temp_c32[,1] else NA) 
Spec_CystCu <<- if(Model.pop=="Normal.pop") temp32[,2] else (if (Model.pop=="CVsurg.pop") temp_c32[,2] else NA) 

Sens_CystCs <<- if(Model.pop=="Normal.pop") temp33[,1] else (if (Model.pop=="CVsurg.pop") temp_c33[,1] else NA) 
Spec_CystCs <<- if(Model.pop=="Normal.pop") temp33[,2] else (if (Model.pop=="CVsurg.pop") temp_c33[,2] else NA) 

# hist(Sens_CystCu, breaks=100)
# hist(Spec_CystCu, breaks=100)
# hist(Sens_CystCs, breaks=100)
# hist(Spec_CystCs, breaks=100)

#length(Sens_CystCu[Sens_CystCu>1]) # Base case=42;  cardio =1
#length(Spec_CystCu[Spec_CystCu>1]) # Base case= 7;  cardio =355
#length(Sens_CystCs[Sens_CystCs>1]) # Base Case=15;  cardio =0
#length(Spec_CystCs[Spec_CystCs>1]) # Base case=0;   cardio =0

#length(Sens_CystCp[Sens_CystCp>1]) # Base Case=0;   cardio =0
#length(Spec_CystCp[Spec_CystCp>1]) # Base case=0;   cardio =0

for (i in 1:Nsim){
  Sens_CystCu[i] <- ifelse(Sens_CystCu[i]>1, 1, Sens_CystCu[i])
  Spec_CystCu[i] <- ifelse(Spec_CystCu[i]>1, 1, Spec_CystCu[i])
  Sens_CystCs[i] <- ifelse(Sens_CystCs[i]>1, 1, Sens_CystCs[i])
  Spec_CystCs[i] <- ifelse(Spec_CystCs[i]>1, 1, Spec_CystCs[i])
}
Sens_CystCu <<- Sens_CystCu 
Spec_CystCu <<- Spec_CystCu
Sens_CystCs <<- Sens_CystCs 
Spec_CystCs <<- Spec_CystCs

#======================= Starting proportions for testing arms =======================# 

#Proportion who arrive with pre-existing AKI                        
#Base case: from communication with LTHT.  #59% of patients who experience AKI in ICU develop AKI before ICU admission 
#Of those, 65% Stage 1, 20% Stage 2, 15% Stage 3 (no data on RRT: assume 5% RRT, 10% stage 3) 

#For testing arms we need to split up the baseline model into patients with pre-existing AKI vs. no AKI on arrival. 
#Baseline proportions for AKI cohort come from 60 patients who have/develop AKI. 0.59*60= 35 already have AKI. 
#0.65*0.59*60= 23.01 stage 1, 0.2*0.59*60= 7.08 S2, 0.1*0.59*60= 3.54 S3, 0.05*0.59*60= 1.77 RRT.
#For nephrocheck, leave S0 and S1 (will both be tested); just take out S2+ (i.e. assume Anyone with S2, S3 or RRT would not be tested, but early cases [S1] would be)

#ASSUME TESTING DONE ON AKI S0-S1 PATIENTS
#Distributions on day 1 for AKI cohort, incorportating proportions not tested: 
day_0_AKI  <- c(18,25,8-7.08,7-3.54,2-1.77) #Distribution on day 0 with pre-existing AKI S2, S3 and RRT cases removed
day_0_preAKI  <- c(7.08,3.54,1.77)
p_day0_AKI  <<- data.frame(rdirichlet(Nsim, day_0_AKI))        #Start in AKI cohort (either have or destined to get AKI)
p_day0_preAKI <<-data.frame(rdirichlet(Nsim, day_0_preAKI))    #Start in AKI cohort, with pre-existing AKI on admission
#rowSums(p_day0_AKI)

#Probability arrive to ICU already with AKI S2+
p_preAKI <<- rbeta(Nsim, 12.39, 60-12.39)*p_AKI     #i.e. 7.08+3.54+1.77 = 12.39

#Starting distributions for AKI cohort, not tested. 
start.S2.preAKI         <<- p_preAKI*p_day0_preAKI[,1]    
start.S3.preAKI         <<- p_preAKI*p_day0_preAKI[,2]    
start.RRT.preAKI        <<- p_preAKI*p_day0_preAKI[,3]


#=== NEPHROCHECK ===#
#Starting distributions for AKI cohort, tested
start.S0.Nephro_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,1]*Sens_Nephro
start.S1.Nephro_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,2]*Sens_Nephro
start.S0.Nephro_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,1]*(1-Sens_Nephro)
start.S1.Nephro_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,2]*(1-Sens_Nephro)
start.S2.Nephro_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,3]*Sens_Nephro
start.S3.Nephro_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,4]*Sens_Nephro
start.RRT.Nephro_TP            <<- (p_AKI-p_preAKI)*p_day0_AKI[,5]*Sens_Nephro
start.S2.Nephro_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,3]*(1-Sens_Nephro)
start.S3.Nephro_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,4]*(1-Sens_Nephro)
start.RRT.Nephro_FN            <<- (p_AKI-p_preAKI)*p_day0_AKI[,5]*(1-Sens_Nephro)

#Starting distributions for no AKI cohort, tested
start.noAKI.Nephro_TN        <<- (1-p_AKI)*Spec_Nephro
start.noAKI.Nephro_FP        <<- (1-p_AKI)*(1-Spec_Nephro)
start.noAKI.Nephro           <<- start.noAKI.Nephro_TN+start.noAKI.Nephro_FP

#=== Sensitivity analysis recoding FP to TP S1 ===#  #DO IN MASTER CODE, HERE FOR REFERENCE ONLY
# start.S1.Nephro_TP <<-start.S1.Nephro_TP+0.1*start.noAKI.Nephro_FP
# #Re-calculte baseline probabilities to increas p_AKI
# p_AKI  <<- p_AKI +0.1*start.noAKI.Nephro_FP     
# day_0  <- c(18,25,8,7,2)  #raw data
# p_day0  <<- data.frame(rdirichlet(Nsim, day_0))
# start.noAKI        <<- 1-p_AKI                 
# start.S0           <<- p_AKI*p_day0[,1]    
# start.S1           <<- p_AKI*p_day0[,2]    
# start.S2           <<- p_AKI*p_day0[,3]    
# start.S3           <<- p_AKI*p_day0[,4]    
# start.RRT          <<- p_AKI*p_day0[,5]  
# #Set FP's to 0 in test arm
#start.noAKI.Nephro_TN        <<- (1-p_AKI)*Spec_Nephro   #NEED THIS AS WELL?
# start.noAKI.Nephro_FP <<- (1-0.1)*start.noAKI.Nephro_FP
# start.noAKI.Nephro    <<- start.noAKI.Nephro_TN+start.noAKI.Nephro_FP

#mean(p_AKI)


#=== NGAL plasma===#
#Starting distributions for AKI cohort, tested
start.S0.NGALp_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,1]*Sens_NGALp
start.S1.NGALp_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,2]*Sens_NGALp
start.S0.NGALp_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,1]*(1-Sens_NGALp)
start.S1.NGALp_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,2]*(1-Sens_NGALp)
start.S2.NGALp_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,3]*Sens_NGALp
start.S3.NGALp_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,4]*Sens_NGALp
start.RRT.NGALp_TP            <<- (p_AKI-p_preAKI)*p_day0_AKI[,5]*Sens_NGALp
start.S2.NGALp_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,3]*(1-Sens_NGALp)
start.S3.NGALp_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,4]*(1-Sens_NGALp)
start.RRT.NGALp_FN            <<- (p_AKI-p_preAKI)*p_day0_AKI[,5]*(1-Sens_NGALp)

#Starting distributions for no AKI cohort, tested
start.noAKI.NGALp_TN        <<- (1-p_AKI)*Spec_NGALp
start.noAKI.NGALp_FP        <<- (1-p_AKI)*(1-Spec_NGALp)
start.noAKI.NGALp           <<- start.noAKI.NGALp_TN+start.noAKI.NGALp_FP

#=== Sensitivity analysis recoding FP to TP S1 ===#  #DO IN MASTER CODE, HERE FOR REFERENCE ONLY
# start.S1.NGALp_TP <<-start.S1.NGALp_TP+0.1*start.noAKI.NGALp_FP
# #Re-calculte baseline probabilities to increas p_AKI
# p_AKI  <<- p_AKI +0.1*start.noAKI.NGALp_FP     
# day_0  <- c(18,25,8,7,2)  #raw data
# p_day0  <<- data.frame(rdirichlet(Nsim, day_0))
# start.noAKI        <<- 1-p_AKI                 
# start.S0           <<- p_AKI*p_day0[,1]    
# start.S1           <<- p_AKI*p_day0[,2]    
# start.S2           <<- p_AKI*p_day0[,3]    
# start.S3           <<- p_AKI*p_day0[,4]    
# start.RRT          <<- p_AKI*p_day0[,5]  
# #Set FP's to 0 in test arm
# start.noAKI.NGALp_FP <<- (1-0.1)*start.noAKI.NGALp_FP
# start.noAKI.NGALp    <<- start.noAKI.NGALp_TN+start.noAKI.NGALp_FP

#mean(p_AKI)


#=== NGAL urine===#
#Starting distributions for AKI cohort, tested
start.S0.NGALu_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,1]*Sens_NGALu
start.S1.NGALu_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,2]*Sens_NGALu
start.S0.NGALu_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,1]*(1-Sens_NGALu)
start.S1.NGALu_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,2]*(1-Sens_NGALu)
start.S2.NGALu_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,3]*Sens_NGALu
start.S3.NGALu_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,4]*Sens_NGALu
start.RRT.NGALu_TP            <<- (p_AKI-p_preAKI)*p_day0_AKI[,5]*Sens_NGALu
start.S2.NGALu_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,3]*(1-Sens_NGALu)
start.S3.NGALu_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,4]*(1-Sens_NGALu)
start.RRT.NGALu_FN            <<- (p_AKI-p_preAKI)*p_day0_AKI[,5]*(1-Sens_NGALu)

#Starting distributions for no AKI cohort, tested
start.noAKI.NGALu_TN        <<- (1-p_AKI)*Spec_NGALu
start.noAKI.NGALu_FP        <<- (1-p_AKI)*(1-Spec_NGALu)
start.noAKI.NGALu           <<- start.noAKI.NGALu_TN+start.noAKI.NGALu_FP

#=== NGAL serum===#
#Starting distributions for AKI cohort, tested
start.S0.NGALs_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,1]*Sens_NGALs
start.S1.NGALs_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,2]*Sens_NGALs
start.S0.NGALs_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,1]*(1-Sens_NGALs)
start.S1.NGALs_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,2]*(1-Sens_NGALs)
start.S2.NGALs_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,3]*Sens_NGALs
start.S3.NGALs_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,4]*Sens_NGALs
start.RRT.NGALs_TP            <<- (p_AKI-p_preAKI)*p_day0_AKI[,5]*Sens_NGALs
start.S2.NGALs_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,3]*(1-Sens_NGALs)
start.S3.NGALs_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,4]*(1-Sens_NGALs)
start.RRT.NGALs_FN            <<- (p_AKI-p_preAKI)*p_day0_AKI[,5]*(1-Sens_NGALs)

#Starting distributions for no AKI cohort, tested
start.noAKI.NGALs_TN        <<- (1-p_AKI)*Spec_NGALs
start.noAKI.NGALs_FP        <<- (1-p_AKI)*(1-Spec_NGALs)
start.noAKI.NGALs           <<- start.noAKI.NGALs_TN+start.noAKI.NGALs_FP

#=== Cystatin C plasma===#
#Starting distributions for AKI cohort, tested
start.S0.CystCp_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,1]*Sens_CystCp    
start.S1.CystCp_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,2]*Sens_CystCp 
start.S0.CystCp_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,1]*(1-Sens_CystCp)    
start.S1.CystCp_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,2]*(1-Sens_CystCp)  
start.S2.CystCp_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,3]*Sens_CystCp 
start.S3.CystCp_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,4]*Sens_CystCp   
start.RRT.CystCp_TP            <<- (p_AKI-p_preAKI)*p_day0_AKI[,5]*Sens_CystCp 
start.S2.CystCp_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,3]*(1-Sens_CystCp)
start.S3.CystCp_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,4]*(1-Sens_CystCp)  
start.RRT.CystCp_FN            <<- (p_AKI-p_preAKI)*p_day0_AKI[,5]*(1-Sens_CystCp) 

#Starting distributions for no AKI cohort, tested 
start.noAKI.CystCp_TN        <<- (1-p_AKI)*Spec_CystCp  
start.noAKI.CystCp_FP        <<- (1-p_AKI)*(1-Spec_CystCp)
start.noAKI.CystCp           <<- start.noAKI.CystCp_TN+start.noAKI.CystCp_FP

#=== Cystatin C urine===#
#Starting distributions for AKI cohort, tested
start.S0.CystCu_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,1]*Sens_CystCu    
start.S1.CystCu_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,2]*Sens_CystCu 
start.S0.CystCu_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,1]*(1-Sens_CystCu)    
start.S1.CystCu_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,2]*(1-Sens_CystCu)  
start.S2.CystCu_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,3]*Sens_CystCu 
start.S3.CystCu_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,4]*Sens_CystCu   
start.RRT.CystCu_TP            <<- (p_AKI-p_preAKI)*p_day0_AKI[,5]*Sens_CystCu 
start.S2.CystCu_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,3]*(1-Sens_CystCu)
start.S3.CystCu_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,4]*(1-Sens_CystCu)  
start.RRT.CystCu_FN            <<- (p_AKI-p_preAKI)*p_day0_AKI[,5]*(1-Sens_CystCu) 

#Starting distributions for no AKI cohort, tested 
start.noAKI.CystCu_TN        <<- (1-p_AKI)*Spec_CystCu  
start.noAKI.CystCu_FP        <<- (1-p_AKI)*(1-Spec_CystCu)
start.noAKI.CystCu           <<- start.noAKI.CystCu_TN+start.noAKI.CystCu_FP

#=== Cystatin C serum===#
#Starting distributions for AKI cohort, tested
start.S0.CystCs_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,1]*Sens_CystCs    
start.S1.CystCs_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,2]*Sens_CystCs 
start.S0.CystCs_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,1]*(1-Sens_CystCs)    
start.S1.CystCs_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,2]*(1-Sens_CystCs)  
start.S2.CystCs_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,3]*Sens_CystCs 
start.S3.CystCs_TP             <<- (p_AKI-p_preAKI)*p_day0_AKI[,4]*Sens_CystCs   
start.RRT.CystCs_TP            <<- (p_AKI-p_preAKI)*p_day0_AKI[,5]*Sens_CystCs 
start.S2.CystCs_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,3]*(1-Sens_CystCs)
start.S3.CystCs_FN             <<- (p_AKI-p_preAKI)*p_day0_AKI[,4]*(1-Sens_CystCs)  
start.RRT.CystCs_FN            <<- (p_AKI-p_preAKI)*p_day0_AKI[,5]*(1-Sens_CystCs) 

#Starting distributions for no AKI cohort, tested 
start.noAKI.CystCs_TN        <<- (1-p_AKI)*Spec_CystCs  
start.noAKI.CystCs_FP        <<- (1-p_AKI)*(1-Spec_CystCs)
start.noAKI.CystCs           <<- start.noAKI.CystCs_TN+start.noAKI.CystCs_FP

#=== Sensitivity analysis recoding FP to TP S1 ===#  #DO IN MASTER CODE, HERE FOR REFERENCE ONLY
# start.S1.CystCs_TP <<-start.S1.CystCs_TP+0.1*start.noAKI.CystCs_FP
# #Re-calculte baseline probabilities to increas p_AKI
# p_AKI  <<- p_AKI +0.1*start.noAKI.CystCs_FP     
# day_0  <- c(18,25,8,7,2)  #raw data
# p_day0  <<- data.frame(rdirichlet(Nsim, day_0))
# start.noAKI        <<- 1-p_AKI                 
# start.S0           <<- p_AKI*p_day0[,1]    
# start.S1           <<- p_AKI*p_day0[,2]    
# start.S2           <<- p_AKI*p_day0[,3]    
# start.S3           <<- p_AKI*p_day0[,4]    
# start.RRT          <<- p_AKI*p_day0[,5]  
# #Set FP's to 0 in test arm
# start.noAKI.CystCs_FP <<- (1-0.1)*start.noAKI.CystCs_FP
# start.noAKI.CystCs    <<- start.noAKI.CystCs_TN+start.noAKI.CystCs_FP

}