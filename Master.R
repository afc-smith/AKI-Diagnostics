#====================================== AKI DIAGNOSTICS ECONOMIC MODEL 2016 ======================================#

# This model and associated files were developed by Alison F Smith from the University of Leeds with supervision from Dr David Meads (Leeds) and Dr Peter Hall (University of Edinburgh)
# This model was developed as part of the NIHR funded Health Technology Assessment (HTA) 'AKI Diagnostics' project (http://www.nets.nihr.ac.uk/projects/hta/1311613)
# This HTA was supported by the National Institute for Health Research (NIHR) Diagnostic Evidence Co-operative Leeds. 
# The views expressed are those of the author(s) and not necessarily those of the NHS, the NIHR or the Department of Health.

#Nb: all lines beginning with a hash (#) are comments 
#To comment/uncomment multiple lines, select the lines and press Ctrl+Shift+C


#=========================================================================================================================#
#====================================== HOSPITAL PHASE MODEL ANALYSIS ====================================================# 
#=========================================================================================================================#

# Set working directory (location where model files are stored)
# NEW USERS: Replace "..." in the below code with the path location where model files are saved e.g. "C:\\User\\My Documents\\AKI Diagnostics")
# Alternatively can use R's project function to group all files in one project folder
setwd("...")

# upload required library packages (nb: if not already installed, these packages will need to be installed first)
library(MASS)     
library(gtools)  
library(plyr)   

 
#======= Load Global Model Data =======# 
 
# remove any stored data
rm(list=ls(all=TRUE))

# set global variables
seed       <- 10             # Set seed for random number generator. Keeping the same seed throughout analyses ensures any differences between results are not due to stochastic variation. 
disc       <- 0.035          # Discount rate for future costs and benefits (as per NICE 2013 methods guide)
Nsim       <- 10000          # Number of Monte Carlo simulations. Set to 10000 in final version.    
threshold  <- 20000          # NICE willingness to pay per QALY threshold (NICE 2013 methods guide)
  
# select which model population you want to run (put hash signs in front of all except the Model that you want to run)
Model.pop      <- "Normal.pop"   # Select to run model for ICU all-comers ('normal') population 
#Model.pop     <- "CVsurg.pop"   # Select to run model for subgroup of patients who are admitted to ICU post cardiovascular surgery



#======= Load Parameters =======#

# Load parameters code
source("Parameters_HOSP.R")

# Draw model parameters Nsim times, to use in the model     
parameter_draw_HOSP()



#=================================== Standard care ('baseline') arm of HOSPITAL model =========================#

# load hospital period model for baseline arm
source("Model_HOSP_baseline.R") 

# Run model
  sim.base <- array(NA,c(Nsim,16))            # Define array table in which costs and QALYs from model will be stored (16 rows relate to the costs, QALYs, and health state distibutions)
      for (i in 1:Nsim) {                     # Run model Nsim times 
        sim.base[i,]     <- model(i)     
      }
  costs.base.HOSP  <- sim.base[,2]             # Store the cost results           
  QALYs.base.HOSP  <- sim.base[,1]             # Store the QALY results    
  trace_hosp.D90   <<- sim.base[,3:16]         # Store the day 90 health state distributions (needed for FUP model)

#Check trace validity (rowSums should all add to 100)      
trace_hosp   <- apply(trace*100, c(1,2), mean)     #*100 to report full % rather than proportions (avoid e+ scientific notation)    
#rowSums(trace_hosp)

#Save results for lifetime analysis
CEresults.HOSP_base <- sim.base[,1:2]

#Plot trace over time (use as a validity check)
# plot(trace_hosp[,1], ylim=c(0, 100), xlim=c(0,90)) #No AKI cohort, In ICU
# lines(trace_hosp[,14])                            #Change value in second vector to plot alternative health states 


#==== Code for model calibration: to set RR_mort ===# 
#This code is kept for reference only - not needed to run the model
#The following ratio was used to identify which value of RR_mort produce ratio of mortalities close to ~4 (see HTA Report for full exaplanation)
#trace_hosp[10,14]/(sum(trace_hosp[10,5:14]))/(trace_hosp[10,4]/(sum(trace_hosp[10,1:4])))   
# Get ratio of mortalities ~=4. Last 10,000 run=4.000226




#================================================ TESTING ARMS OF HOSPITAL MODEL ============================================#

# Note: Model intervention parameter 'test':
# test = 1 for Nephrocheck test
# test = 21 for NGAL plasma
# test = 22 for NGAL urine
# test = 23 for NGAL serum
# test = 31 for Cycstatin C plasma
# test = 32 for Cycstatin C urine
# test = 33 for Cycstatin C serum

#For the testing arms of the model, the cohort is split into:
   #True positives (TP) receiving early treatment
   #TPs not receiving early treatment
   #True negatives (TN) + False Positives (FP) [i.e. 'no AKI' cohort]
   #Patients arriving to ICU with existing AKI and therefore not tested


#===================================== NEPHROCHECK ========================================#

#Set test value
test  <- 1

#=== Run TP arm ===#
source("Model_HOSP_test_TP_treat.R")

sim.test <- array(NA,c(Nsim,16))
  for (i in 1:Nsim) {
    sim.test[i,]     <- model(i,test)
  }
costs.Nephro.HOSP_TPtemp <- sim.test[,2]
QALYs.Nephro.HOSP_TPtemp <- sim.test[,1]
trace.D90.Nephro_TPtemp  <- sim.test[,3:16]

trace_temp_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp_full)  ##Nb: these won't add to 100 now, as the cohort is split across TP/FP etc.


source("Model_HOSP_test_TP_notreat.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.Nephro.HOSP_TP <- sim.test[,2] + costs.Nephro.HOSP_TPtemp
QALYs.Nephro.HOSP_TP <- sim.test[,1] + QALYs.Nephro.HOSP_TPtemp
trace.D90.Nephro_TP  <- sim.test[,3:16] + trace.D90.Nephro_TPtemp

trace_temp1_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp1_full)  


#=== Run FN arm ===#
source("Model_HOSP_test_FN.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.Nephro.HOSP_FN <- sim.test[,2]
QALYs.Nephro.HOSP_FN <- sim.test[,1]
trace.D90.Nephro_FN  <- sim.test[,3:16]

trace_temp2_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp2_full)  


#=== Run TN+FP (i.e. 'No AKI') arm ===#
source("Model_HOSP_test_noAKI.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.Nephro.HOSP_noAKI <- sim.test[,2]
QALYs.Nephro.HOSP_noAKI <- sim.test[,1]
trace.D90.Nephro_noAKI  <- sim.test[,3:16]

trace_temp3_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp3_full)  


#=== Run pre-existing AKI arm (i.e.'notest') ===#
source("Model_HOSP_notest.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.Nephro.HOSP_notest <- sim.test[,2]
QALYs.Nephro.HOSP_notest <- sim.test[,1]
trace.D90.Nephro_notest  <- sim.test[,3:16]

trace_temp4_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp4_full)  



#=== Combine outcomes ===#
costs.Nephro.HOSP      <- costs.Nephro.HOSP_TP + costs.Nephro.HOSP_FN + costs.Nephro.HOSP_noAKI + costs.Nephro.HOSP_notest
QALYs.Nephro.HOSP      <- QALYs.Nephro.HOSP_TP + QALYs.Nephro.HOSP_FN + QALYs.Nephro.HOSP_noAKI + QALYs.Nephro.HOSP_notest
CEresults.HOSP_Nephro  <<- data.frame(QALYs.Nephro.HOSP, costs.Nephro.HOSP)
trace.D90.Nephro       <<- trace.D90.Nephro_TP + trace.D90.Nephro_FN + trace.D90.Nephro_noAKI + trace.D90.Nephro_notest
trace_hosp_Nephro      <<- trace_temp_full +trace_temp1_full +trace_temp2_full +trace_temp3_full +trace_temp4_full


#Check full trace & explore differences
#rowSums(trace_hosp_Nephro)  #Should all be 100.
#rowSums(trace.D90.Nephro) #should all be 1.

# plot(trace_hosp_Nephro[,14]+trace_hosp_Nephro[,4], ylim=c(0, 100), xlim=c(0,90))    #Compare Nephro and base arms
# lines(trace_hosp[,14]+trace_hosp[,4])
#
# # #End distribution (as % of no AKI or AKI cohorts)
#  trace_hosp_Nephro[90,4]/sum(trace_hosp_Nephro[90,1:4])




#===================================== NGAL plasma ========================================#

#Set test value
test  <- 21

#=== Run TP arm ===#
source("Model_HOSP_test_TP_treat.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALp.HOSP_TPtemp <- sim.test[,2]
QALYs.NGALp.HOSP_TPtemp <- sim.test[,1]
trace.D90.NGALp_TPtemp  <- sim.test[,3:16]

trace_temp_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp_full)  #23.59101


source("Model_HOSP_test_TP_notreat.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALp.HOSP_TP <- sim.test[,2] + costs.NGALp.HOSP_TPtemp
QALYs.NGALp.HOSP_TP <- sim.test[,1] + QALYs.NGALp.HOSP_TPtemp
trace.D90.NGALp_TP  <- sim.test[,3:16] + trace.D90.NGALp_TPtemp

trace_temp1_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp1_full)  #2.454827


#=== Run FN arm ===#
source("Model_HOSP_test_FN.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALp.HOSP_FN <- sim.test[,2]
QALYs.NGALp.HOSP_FN <- sim.test[,1]
trace.D90.NGALp_FN  <- sim.test[,3:16]

trace_temp2_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp2_full)  #3.5687


#=== Run TN+FP (i.e. 'No AKI') arm ===#
source("Model_HOSP_test_noAKI.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALp.HOSP_noAKI <- sim.test[,2]
QALYs.NGALp.HOSP_noAKI <- sim.test[,1]
trace.D90.NGALp_noAKI  <- sim.test[,3:16]

trace_temp3_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp3_full)  #63.81555


#=== Run pre-existing AKI arm (i.e.'notest') ===#
source("Model_HOSP_notest.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALp.HOSP_notest <- sim.test[,2]
QALYs.NGALp.HOSP_notest <- sim.test[,1]
trace.D90.NGALp_notest  <- sim.test[,3:16]

trace_temp4_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp4_full)  #6.569913


#=== Combine outcomes ===#
costs.NGALp.HOSP      <- costs.NGALp.HOSP_TP + costs.NGALp.HOSP_FN + costs.NGALp.HOSP_noAKI + costs.NGALp.HOSP_notest
QALYs.NGALp.HOSP      <- QALYs.NGALp.HOSP_TP + QALYs.NGALp.HOSP_FN + QALYs.NGALp.HOSP_noAKI + QALYs.NGALp.HOSP_notest
CEresults.HOSP_NGALp  <<- data.frame(QALYs.NGALp.HOSP, costs.NGALp.HOSP)
trace.D90.NGALp       <<- trace.D90.NGALp_TP + trace.D90.NGALp_FN + trace.D90.NGALp_noAKI + trace.D90.NGALp_notest
trace_hosp_NGALp      <<- trace_temp_full +trace_temp1_full +trace_temp2_full +trace_temp3_full +trace_temp4_full

#rowSums(trace_hosp_NGALp)  #Should all be 100.  NOT ADDING TO 100!!! 
#rowSums(trace.D90.NGALp) #should all be 1. 

#plot(trace_hosp_NGALp[,14]+trace_hosp_NGALp[,4], ylim=c(0, 100), xlim=c(0,90))    #Compare Nephro and base arms
#lines(trace_hosp[,14]+trace_hosp[,4]) 



#===================================== NGAL urine ========================================#

#Set test value
test  <- 22

#=== Run TP arm ===#
source("Model_HOSP_test_TP_treat.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALu.HOSP_TPtemp <- sim.test[,2]
QALYs.NGALu.HOSP_TPtemp <- sim.test[,1]
trace.D90.NGALu_TPtemp  <- sim.test[,3:16]

trace_temp_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp_full)  #23.59101


source("Model_HOSP_test_TP_notreat.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALu.HOSP_TP <- sim.test[,2] + costs.NGALu.HOSP_TPtemp
QALYs.NGALu.HOSP_TP <- sim.test[,1] + QALYs.NGALu.HOSP_TPtemp
trace.D90.NGALu_TP  <- sim.test[,3:16] + trace.D90.NGALu_TPtemp

trace_temp1_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp1_full)  #2.454827


#=== Run FN arm ===#
source("Model_HOSP_test_FN.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALu.HOSP_FN <- sim.test[,2]
QALYs.NGALu.HOSP_FN <- sim.test[,1]
trace.D90.NGALu_FN  <- sim.test[,3:16]

trace_temp2_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp2_full)  #3.5687


#=== Run TN+FP (i.e. 'No AKI') arm ===#
source("Model_HOSP_test_noAKI.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALu.HOSP_noAKI <- sim.test[,2]
QALYs.NGALu.HOSP_noAKI <- sim.test[,1]
trace.D90.NGALu_noAKI  <- sim.test[,3:16]

trace_temp3_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp3_full)  #63.81555


#=== Run pre-existing AKI arm (i.e.'notest') ===#
source("Model_HOSP_notest.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALu.HOSP_notest <- sim.test[,2]
QALYs.NGALu.HOSP_notest <- sim.test[,1]
trace.D90.NGALu_notest  <- sim.test[,3:16]

trace_temp4_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp4_full)  #6.569913


#=== Combine outcomes ===#
costs.NGALu.HOSP      <- costs.NGALu.HOSP_TP + costs.NGALu.HOSP_FN + costs.NGALu.HOSP_noAKI + costs.NGALu.HOSP_notest
QALYs.NGALu.HOSP      <- QALYs.NGALu.HOSP_TP + QALYs.NGALu.HOSP_FN + QALYs.NGALu.HOSP_noAKI + QALYs.NGALu.HOSP_notest
CEresults.HOSP_NGALu  <<- data.frame(QALYs.NGALu.HOSP, costs.NGALu.HOSP)
trace.D90.NGALu       <<- trace.D90.NGALu_TP + trace.D90.NGALu_FN + trace.D90.NGALu_noAKI + trace.D90.NGALu_notest
trace_hosp_NGALu      <<- trace_temp_full +trace_temp1_full +trace_temp2_full +trace_temp3_full +trace_temp4_full

# rowSums(trace_hosp_NGALu)  #Should all be 100.
# rowSums(trace.D90.NGALu) #should all be 1.
#
# #plot(trace_hosp_NGALu[,14]+trace_hosp_NGALu[,4], ylim=c(0, 100), xlim=c(0,90))    #Compare Nephro and base arms
# #lines(trace_hosp[,14]+trace_hosp[,4])


#===================================== NGAL serum ========================================#

#Set test value
test  <- 23

#=== Run TP arm ===#
source("Model_HOSP_test_TP_treat.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALs.HOSP_TPtemp <- sim.test[,2]
QALYs.NGALs.HOSP_TPtemp <- sim.test[,1]
trace.D90.NGALs_TPtemp  <- sim.test[,3:16]

trace_temp_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp_full)  #23.59101


source("Model_HOSP_test_TP_notreat.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALs.HOSP_TP <- sim.test[,2] + costs.NGALs.HOSP_TPtemp
QALYs.NGALs.HOSP_TP <- sim.test[,1] + QALYs.NGALs.HOSP_TPtemp
trace.D90.NGALs_TP  <- sim.test[,3:16] + trace.D90.NGALs_TPtemp

trace_temp1_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp1_full)  #2.454827


#=== Run FN arm ===#
source("Model_HOSP_test_FN.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALs.HOSP_FN <- sim.test[,2]
QALYs.NGALs.HOSP_FN <- sim.test[,1]
trace.D90.NGALs_FN  <- sim.test[,3:16]

trace_temp2_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp2_full)  #3.5687


#=== Run TN+FP (i.e. 'No AKI') arm ===#
source("Model_HOSP_test_noAKI.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALs.HOSP_noAKI <- sim.test[,2]
QALYs.NGALs.HOSP_noAKI <- sim.test[,1]
trace.D90.NGALs_noAKI  <- sim.test[,3:16]

trace_temp3_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp3_full)  #63.81555


#=== Run pre-existing AKI arm (i.e.'notest') ===#
source("Model_HOSP_notest.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALs.HOSP_notest <- sim.test[,2]
QALYs.NGALs.HOSP_notest <- sim.test[,1]
trace.D90.NGALs_notest  <- sim.test[,3:16]

trace_temp4_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp4_full)  #6.569913


#=== Combine outcomes ===#
costs.NGALs.HOSP      <- costs.NGALs.HOSP_TP + costs.NGALs.HOSP_FN + costs.NGALs.HOSP_noAKI + costs.NGALs.HOSP_notest
QALYs.NGALs.HOSP      <- QALYs.NGALs.HOSP_TP + QALYs.NGALs.HOSP_FN + QALYs.NGALs.HOSP_noAKI + QALYs.NGALs.HOSP_notest
CEresults.HOSP_NGALs  <<- data.frame(QALYs.NGALs.HOSP, costs.NGALs.HOSP)
trace.D90.NGALs       <<- trace.D90.NGALs_TP + trace.D90.NGALs_FN + trace.D90.NGALs_noAKI + trace.D90.NGALs_notest
trace_hosp_NGALs      <<- trace_temp_full +trace_temp1_full +trace_temp2_full +trace_temp3_full +trace_temp4_full

# rowSums(trace_hosp_NGALs)  #Should all be 100.  NOT ADDING TO 100!!! 
# rowSums(trace.D90.NGALs) #should all be 1. 
# 
# #plot(trace_hosp_NGALs[,14]+trace_hosp_NGALs[,4], ylim=c(0, 100), xlim=c(0,90))    #Compare Nephro and base arms
# #lines(trace_hosp[,14]+trace_hosp[,4]) 


#======================================= Cystatin C plasma ==========================================#     

test  <- 31

#=== Run TP arm ===#
source("Model_HOSP_test_TP_treat.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.CystCp.HOSP_TPtemp <- sim.test[,2]
QALYs.CystCp.HOSP_TPtemp <- sim.test[,1]
trace.D90.CystCp_TPtemp  <- sim.test[,3:16]

trace_temp_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp_full)  #23.59101


source("Model_HOSP_test_TP_notreat.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.CystCp.HOSP_TP <- sim.test[,2] + costs.CystCp.HOSP_TPtemp
QALYs.CystCp.HOSP_TP <- sim.test[,1] + QALYs.CystCp.HOSP_TPtemp
trace.D90.CystCp_TP  <- sim.test[,3:16] + trace.D90.CystCp_TPtemp

trace_temp1_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp1_full)  #2.454827


#=== Run FN arm ===#
source("Model_HOSP_test_FN.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.CystCp.HOSP_FN <- sim.test[,2]
QALYs.CystCp.HOSP_FN <- sim.test[,1]
trace.D90.CystCp_FN  <- sim.test[,3:16]

trace_temp2_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp2_full)  #3.5687


#=== Run TN+FP (i.e. 'No AKI') arm ===#
source("Model_HOSP_test_noAKI.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.CystCp.HOSP_noAKI <- sim.test[,2]
QALYs.CystCp.HOSP_noAKI <- sim.test[,1]
trace.D90.CystCp_noAKI  <- sim.test[,3:16]

trace_temp3_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp3_full)  #63.81555


#=== Run pre-existing AKI arm (i.e.'notest') ===#
source("Model_HOSP_notest.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.CystCp.HOSP_notest <- sim.test[,2]
QALYs.CystCp.HOSP_notest <- sim.test[,1]
trace.D90.CystCp_notest  <- sim.test[,3:16]

trace_temp4_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp4_full)  #6.569913


#=== Combine outcomes ===#
costs.CystCp.HOSP      <- costs.CystCp.HOSP_TP + costs.CystCp.HOSP_FN + costs.CystCp.HOSP_noAKI + costs.CystCp.HOSP_notest
QALYs.CystCp.HOSP      <- QALYs.CystCp.HOSP_TP + QALYs.CystCp.HOSP_FN + QALYs.CystCp.HOSP_noAKI + QALYs.CystCp.HOSP_notest
CEresults.HOSP_CystCp  <<- data.frame(QALYs.CystCp.HOSP, costs.CystCp.HOSP)
trace.D90.CystCp       <<- trace.D90.CystCp_TP + trace.D90.CystCp_FN + trace.D90.CystCp_noAKI + trace.D90.CystCp_notest
trace_hosp_CystCp      <<- trace_temp_full +trace_temp1_full +trace_temp2_full +trace_temp3_full +trace_temp4_full


#======================================= Cystatin C urine ==========================================#     

test  <- 32

#=== Run TP arm ===#
source("Model_HOSP_test_TP_treat.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.CystCu.HOSP_TPtemp <- sim.test[,2]
QALYs.CystCu.HOSP_TPtemp <- sim.test[,1]
trace.D90.CystCu_TPtemp  <- sim.test[,3:16]

trace_temp_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp_full)  #23.59101


source("Model_HOSP_test_TP_notreat.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.CystCu.HOSP_TP <- sim.test[,2] + costs.CystCu.HOSP_TPtemp
QALYs.CystCu.HOSP_TP <- sim.test[,1] + QALYs.CystCu.HOSP_TPtemp
trace.D90.CystCu_TP  <- sim.test[,3:16] + trace.D90.CystCu_TPtemp

trace_temp1_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp1_full)  #2.454827


#=== Run FN arm ===#
source("Model_HOSP_test_FN.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.CystCu.HOSP_FN <- sim.test[,2]
QALYs.CystCu.HOSP_FN <- sim.test[,1]
trace.D90.CystCu_FN  <- sim.test[,3:16]

trace_temp2_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp2_full)  #3.5687


#=== Run TN+FP (i.e. 'No AKI') arm ===#
source("Model_HOSP_test_noAKI.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.CystCu.HOSP_noAKI <- sim.test[,2]
QALYs.CystCu.HOSP_noAKI <- sim.test[,1]
trace.D90.CystCu_noAKI  <- sim.test[,3:16]

trace_temp3_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp3_full)  #63.81555


#=== Run pre-existing AKI arm (i.e.'notest') ===#
source("Model_HOSP_notest.R")

sim.test <- array(NA,c(Nsim,16))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.CystCu.HOSP_notest <- sim.test[,2]
QALYs.CystCu.HOSP_notest <- sim.test[,1]
trace.D90.CystCu_notest  <- sim.test[,3:16]

trace_temp4_full   <- apply(trace*100, c(1,2), mean)
#rowSums(trace_temp4_full)  #6.569913


#=== Combine outcomes ===#
costs.CystCu.HOSP      <- costs.CystCu.HOSP_TP + costs.CystCu.HOSP_FN + costs.CystCu.HOSP_noAKI + costs.CystCu.HOSP_notest
QALYs.CystCu.HOSP      <- QALYs.CystCu.HOSP_TP + QALYs.CystCu.HOSP_FN + QALYs.CystCu.HOSP_noAKI + QALYs.CystCu.HOSP_notest
CEresults.HOSP_CystCu  <<- data.frame(QALYs.CystCu.HOSP, costs.CystCu.HOSP)
trace.D90.CystCu       <<- trace.D90.CystCu_TP + trace.D90.CystCu_FN + trace.D90.CystCu_noAKI + trace.D90.CystCu_notest
trace_hosp_CystCu      <<- trace_temp_full +trace_temp1_full +trace_temp2_full +trace_temp3_full +trace_temp4_full

 
#======================================= Cystatin C serum ==========================================#     

#Set test value
test  <- 33                                                           

#=== Run TP arm ===# 
source("Model_HOSP_test_TP_treat.R")                                             

sim.test <- array(NA,c(Nsim,16))               
for (i in 1:Nsim) {                     
  sim.test[i,]     <- model(i,test)   
}
costs.CystCs.HOSP_TPtemp <- sim.test[,2]               
QALYs.CystCs.HOSP_TPtemp <- sim.test[,1]
trace.D90.CystCs_TPtemp  <- sim.test[,3:16]

trace_temp_full   <- apply(trace*100, c(1,2), mean)   
#rowSums(trace_temp_full)  #23.59101


source("Model_HOSP_test_TP_notreat.R")                                             

sim.test <- array(NA,c(Nsim,16))               
for (i in 1:Nsim) {                     
  sim.test[i,]     <- model(i,test)   
}
costs.CystCs.HOSP_TP <- sim.test[,2] + costs.CystCs.HOSP_TPtemp          
QALYs.CystCs.HOSP_TP <- sim.test[,1] + QALYs.CystCs.HOSP_TPtemp
trace.D90.CystCs_TP  <- sim.test[,3:16] + trace.D90.CystCs_TPtemp

trace_temp1_full   <- apply(trace*100, c(1,2), mean)   
#rowSums(trace_temp1_full)  #2.454827


#=== Run FN arm ===#
source("Model_HOSP_test_FN.R")                                          

sim.test <- array(NA,c(Nsim,16))               
for (i in 1:Nsim) {                     
  sim.test[i,]     <- model(i,test)   
}
costs.CystCs.HOSP_FN <- sim.test[,2]               
QALYs.CystCs.HOSP_FN <- sim.test[,1]
trace.D90.CystCs_FN  <- sim.test[,3:16]

trace_temp2_full   <- apply(trace*100, c(1,2), mean)   
#rowSums(trace_temp2_full)  #3.5687


#=== Run TN+FP (i.e. 'No AKI') arm ===#
source("Model_HOSP_test_noAKI.R")                                          

sim.test <- array(NA,c(Nsim,16))               
for (i in 1:Nsim) {                     
  sim.test[i,]     <- model(i,test)   
}
costs.CystCs.HOSP_noAKI <- sim.test[,2]               
QALYs.CystCs.HOSP_noAKI <- sim.test[,1]
trace.D90.CystCs_noAKI  <- sim.test[,3:16]

trace_temp3_full   <- apply(trace*100, c(1,2), mean)   
#rowSums(trace_temp3_full)  #63.81555


#=== Run pre-existing AKI arm (i.e.'notest') ===#
source("Model_HOSP_notest.R")                                          

sim.test <- array(NA,c(Nsim,16))               
for (i in 1:Nsim) {                     
  sim.test[i,]     <- model(i,test)   
}
costs.CystCs.HOSP_notest <- sim.test[,2]               
QALYs.CystCs.HOSP_notest <- sim.test[,1]
trace.D90.CystCs_notest  <- sim.test[,3:16]

trace_temp4_full   <- apply(trace*100, c(1,2), mean)   
#rowSums(trace_temp4_full)  #6.569913


#=== Combine outcomes ===#
costs.CystCs.HOSP      <- costs.CystCs.HOSP_TP + costs.CystCs.HOSP_FN + costs.CystCs.HOSP_noAKI + costs.CystCs.HOSP_notest          
QALYs.CystCs.HOSP      <- QALYs.CystCs.HOSP_TP + QALYs.CystCs.HOSP_FN + QALYs.CystCs.HOSP_noAKI + QALYs.CystCs.HOSP_notest 
CEresults.HOSP_CystCs  <<- data.frame(QALYs.CystCs.HOSP, costs.CystCs.HOSP)
trace.D90.CystCs       <<- trace.D90.CystCs_TP + trace.D90.CystCs_FN + trace.D90.CystCs_noAKI + trace.D90.CystCs_notest 
trace_hosp_CystCs      <<- trace_temp_full +trace_temp1_full +trace_temp2_full +trace_temp3_full +trace_temp4_full 




#=========================================================================================================================#
#============================================= FOLLOWUP PHASE MODEL ANALYSIS =============================================# 
#=========================================================================================================================#

#======= Load Model & Parameters =======#

source("Parameters_FUP.R")      

#======== Run standard care ('baseline') arm of model =======#

#Set test (have to do before running parameter draw for FUP model)
test <- 0  

#Draw model parameters Nsim times     
parameter_draw_FUP()

#Load follow up model
source("Model_FUP.R")

  sim.base <- array(NA,c(Nsim,2))
      for (i in 1:Nsim) {
        sim.base[i,]     <- model(i,test)       
      }
  costs.base.FUP <- sim.base[,2]
  QALYs.base.FUP <- sim.base[,1]       
  trace_FUP <<- apply(trace*100, c(1,2), mean)

#Save outputs
CEresults.FUP_base  <<- sim.base

#Check Trace
#rowSums(trace_FUP)
# plot(trace_FUP[,1], ylim=c(0, 100), xlim=c(0,40))                    
# lines(trace_FUP[,2])                                              

# #Final states as proportions of AKI/no AKI cohorts
# trace_FUP[T,2]/(trace_FUP[T,2]+trace_FUP[T,1])  



#=============== Run Nephrocheck arms of model ===============#

#Set test value
test <- 1  

#Load code for test-specific parameters
source("Parameters_FUP_tests.R")

#Draw model parameters Nsim times     
parameter_draw_FUP()

#Load follow up model
source("Model_FUP.R")

  sim.test <- array(NA,c(Nsim,2))
      for (i in 1:Nsim) {
        sim.test[i,]     <- model(i,test)   
      }
  costs.Nephro.FUP  <- sim.test[,2] 
  QALYs.Nephro.FUP  <- sim.test[,1] 
  trace_Nephro_FUP <<- apply(trace*100, c(1,2), mean)

#Save outputs
CEresults.FUP_Nephro <<- sim.test
#rowSums(trace_Nephro_FUP)


#=============== Run NGAL arms of model ===============#

#=== PLASMA TEST ===#
#Set test value
test <- 21  

#Load code for test-specific parameters
source("Parameters_FUP_tests.R")

#Draw model parameters Nsim times     
parameter_draw_FUP()               #Error in trace_hosp[i, 1] : incorrect number of dimensions

#Load follow up model
source("Model_FUP.R")

sim.test <- array(NA,c(Nsim,2))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)   
}
costs.NGALp.FUP  <- sim.test[,2] 
QALYs.NGALp.FUP  <- sim.test[,1] 
trace_NGALp_FUP <<- apply(trace*100, c(1,2), mean)

#Save outputs
CEresults.FUP_NGALp <<- sim.test
#write.csv(CEresults.FUP_Nephro ,"Model Outputs\\FUP Model outputs\\Outputs_Nephro.csv",row.names=FALSE)
#rowSums(trace_NGALp_FUP)


#=== URINE TEST ===#

#Set test value
test <- 22

#Load code for test-specific parameters
source("Parameters_FUP_tests.R")

#Draw model parameters Nsim times
parameter_draw_FUP()

#Load follow up model
source("Model_FUP.R")

sim.test <- array(NA,c(Nsim,2))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALu.FUP  <- sim.test[,2]
QALYs.NGALu.FUP  <- sim.test[,1]
trace_NGALu_FUP <<- apply(trace*100, c(1,2), mean)

#Save outputs
CEresults.FUP_NGALu <<- sim.test
#rowSums(trace_NGALu_FUP)

#=== SERUM TEST ===#
#Set test value
test <- 23

#Load code for test-specific parameters
source("Parameters_FUP_tests.R")

#Draw model parameters Nsim times
parameter_draw_FUP()

#Load follow up model
source("Model_FUP.R")

sim.test <- array(NA,c(Nsim,2))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.NGALs.FUP  <- sim.test[,2]
QALYs.NGALs.FUP  <- sim.test[,1]
trace_NGALs_FUP <<- apply(trace*100, c(1,2), mean)

#Save outputs
CEresults.FUP_NGALs <<- sim.test
#rowSums(trace_NGALs_FUP)


#=============== Run Cystatin C arms of model ===============#

#==PLASMA TEST===#

#Set test value
test <- 31

#Load code for test-specific parameters
source("Parameters_FUP_tests.R")

#Draw model parameters Nsim times
parameter_draw_FUP()

#Load follow up model
source("Model_FUP.R")

sim.test <- array(NA,c(Nsim,2))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.CystCp.FUP  <- sim.test[,2]
QALYs.CystCp.FUP  <- sim.test[,1]
trace_CystCp_FUP <<- apply(trace*100, c(1,2), mean)

#Save outputs
CEresults.FUP_CystCp <<- sim.test


#==URINE TEST===#

#Set test value
test <- 32

#Load code for test-specific parameters
source("Parameters_FUP_tests.R")

#Draw model parameters Nsim times
parameter_draw_FUP()

#Load follow up model
source("Model_FUP.R")

sim.test <- array(NA,c(Nsim,2))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.CystCu.FUP  <- sim.test[,2]
QALYs.CystCu.FUP  <- sim.test[,1]
trace_CystCu_FUP <<- apply(trace*100, c(1,2), mean)

#Save outputs
CEresults.FUP_CystCu <<- sim.test


#==SERUM TEST===#

#Set test value
test <- 33

#Load code for test-specific parameters
source("Parameters_FUP_tests.R")

#Draw model parameters Nsim times
parameter_draw_FUP()

#Load follow up model
source("Model_FUP.R")

sim.test <- array(NA,c(Nsim,2))
for (i in 1:Nsim) {
  sim.test[i,]     <- model(i,test)
}
costs.CystCs.FUP  <- sim.test[,2]
QALYs.CystCs.FUP  <- sim.test[,1]
trace_CystCs_FUP <<- apply(trace*100, c(1,2), mean)

#Save outputs
CEresults.FUP_CystCs <<- sim.test



#================== Get list of model parameters =====================#

#Run code to store all model parameters (needed for EVPPI & EVSI analyses)
source("theta_code.R")
theta.code()   #6567 parameters
#write.csv(theta_HOSP, "theta_HOSP_CystCs.csv")   #change name of saved file as required
#write.csv(theta_FUP, "theta_FUP_CystCs.csv")
#write.csv(theta_HOSP_medium, "theta_HOSP_medium.csv")
#write.csv(theta_HOSP_large, "theta_HOSP_large.csv")

#theta_HOSP <- read.csv("theta_HOSP.csv")  #Use 'read.csv' code to read in saved files, if running at a later date rather than directly after running the model
#theta_FUP <- read.csv("theta_HOSP.csv")


#=========================================================================================================================#
#============================================== TWO-WAY LIFETIME ANALYSIS ================================================#
#=========================================================================================================================#

#============ NEPHROCHECK VS. STANDARD CARE ===========#

#=== Save hops & FUP results ===#
Nephro_split_results  <- data.frame(CEresults.HOSP_base, CEresults.HOSP_Nephro, CEresults.FUP_base, CEresults.FUP_Nephro)
colnames(Nephro_split_results)  <- c("Hosp_Q_base", "Hosp_C_base", "Hosp_Q_Nephro", "Hosp_C_Nephro", "FUP_Q_base", "FUP_C_base", "FUP_Q_Nephro", "FUP_C_Nephro")

#=== Combine HOSP and FUP model results ===# 
CEresults.life_base    <- CEresults.HOSP_base + CEresults.FUP_base  
CEresults.life_Nephro  <- CEresults.HOSP_Nephro +CEresults.FUP_Nephro

#=== Two- way Incremental Net Health Benefit analysis ===#
threshold  <- 20000
NB.base    <- CEresults.life_base[,1] - CEresults.life_base[,2]/threshold
NB.Nephro  <- CEresults.life_Nephro[,1] - CEresults.life_Nephro[,2]/threshold
INB_Nephro  <- NB.Nephro-NB.base
mean(INB_Nephro)

#=== Compile and save results ===# 
Nephro_results <- data.frame(CEresults.life_base, CEresults.life_Nephro, NB.base, NB.Nephro, INB_Nephro)
Nephro_results[,8]<-  Nephro_results[,3]- Nephro_results[,1]
Nephro_results[,9]<-   Nephro_results[,4]- Nephro_results[,2]
Nephro_results[,10]  <- ifelse(Nephro_results[,8]>0,1,0)
Nephro_results[,11]  <- ifelse(Nephro_results[,9]<0,1,0)
colnames(Nephro_results) <- c("QALYs_base","Costs_base","QALYs_Nephro","Costs_Nephro","NB_base","NB_Nephro","INB", "Incr.QALY", "Incr.Cost", "p_eff", "p_save")
Nephro_results <- data.frame(Nephro_results, Nephro_split_results)
#write.csv(Nephro_results, "Nephro_basecase.csv")

# #==== 95% CIs ===#
# Conf_ints <- array(NA, c(length(Nephro_results[1,]),2))
# for (i in 1:19){
#   Conf_ints[i,1] <- quantile(Nephro_results[,i],c(0.025))
#   Conf_ints[i,2] <- quantile(Nephro_results[,i], c(0.975))
# }
# write.csv(Conf_ints, "Conf_ints_cost_saving_100.csv")

#=== Traces ====#
# write.csv(trace_hosp, "trace_hosp.csv")
# write.csv(trace_FUP, "trace_FUP.csv")
# write.csv(trace_hosp_Nephro, "trace_hosp_Nephro.csv")
# write.csv(trace_Nephro_FUP, "trace_Nephro_FUP.csv")


#============ NGAL VS. STANDARD CARE ===========#

#=== Save hops & FUP results ===#
 NGALp_split_results  <- data.frame(CEresults.HOSP_base, CEresults.HOSP_NGALp, CEresults.FUP_base, CEresults.FUP_NGALp)
 colnames(NGALp_split_results)  <- c("Hosp_Q_base", "Hosp_C_base", "Hosp_Q_NGALp", "Hosp_C_NGALp", "FUP_Q_base", "FUP_C_base", "FUP_Q_NGALp", "FUP_C_NGALp")
# 
 NGALu_split_results  <- data.frame(CEresults.HOSP_base, CEresults.HOSP_NGALu, CEresults.FUP_base, CEresults.FUP_NGALu)
 colnames(NGALu_split_results)  <- c("Hosp_Q_base", "Hosp_C_base", "Hosp_Q_NGALu", "Hosp_C_NGALu", "FUP_Q_base", "FUP_C_base", "FUP_Q_NGALu", "FUP_C_NGALu")

 NGALs_split_results  <- data.frame(CEresults.HOSP_base, CEresults.HOSP_NGALs, CEresults.FUP_base, CEresults.FUP_NGALs)
 colnames(NGALs_split_results)  <- c("Hosp_Q_base", "Hosp_C_base", "Hosp_Q_NGALs", "Hosp_C_NGALs", "FUP_Q_base", "FUP_C_base", "FUP_Q_NGALs", "FUP_C_NGALs")


# #=== Combine HOSP and FUP model results ===# 
 CEresults.life_base    <- CEresults.HOSP_base + CEresults.FUP_base 
 CEresults.life_NGALp   <- CEresults.HOSP_NGALp +CEresults.FUP_NGALp
 CEresults.life_NGALu   <- CEresults.HOSP_NGALu +CEresults.FUP_NGALu
 CEresults.life_NGALs   <- CEresults.HOSP_NGALs +CEresults.FUP_NGALs

# #=== Two- way Incremental Net Health Benefit analysis ===#
 threshold  <- 20000
 NB.base    <- CEresults.life_base[,1] - CEresults.life_base[,2]/threshold
 NB.NGALp   <- CEresults.life_NGALp[,1] - CEresults.life_NGALp[,2]/threshold
 INB_NGALp  <- NB.NGALp-NB.base
 NB.NGALu   <- CEresults.life_NGALu[,1] - CEresults.life_NGALu[,2]/threshold
 INB_NGALu  <- NB.NGALu-NB.base
 NB.NGALs   <- CEresults.life_NGALs[,1] - CEresults.life_NGALs[,2]/threshold
 INB_NGALs  <- NB.NGALs-NB.base

#=== Compile and save results ===# 
NGALp_results       <- data.frame(CEresults.life_base, CEresults.life_NGALp, NB.base, NB.NGALp, INB_NGALp)
NGALp_results[,8]   <-  NGALp_results[,3]- NGALp_results[,1]
NGALp_results[,9]   <-  NGALp_results[,4]- NGALp_results[,2]
NGALp_results[,10]  <- ifelse(NGALp_results[,8]>0,1,0)
NGALp_results[,11]  <- ifelse(NGALp_results[,9]<0,1,0)
colnames(NGALp_results) <- c("QALYs_base","Costs_base","QALYs_NGALp","Costs_NGALp","NB_base","NB_NGALp","INB", "Incr.QALY", "Incr.Cost", "p_eff", "p_save")
NGALp_results       <- data.frame(NGALp_results, NGALp_split_results) 
#write.csv(NGALp_results, "NGALp_results_basecase.csv")
 
NGALu_results       <- data.frame(CEresults.life_base, CEresults.life_NGALu, NB.base, NB.NGALu, INB_NGALu)
NGALu_results[,8]   <-  NGALu_results[,3]- NGALu_results[,1]
NGALu_results[,9]   <-   NGALu_results[,4]- NGALu_results[,2]
NGALu_results[,10]  <- ifelse(NGALu_results[,8]>0,1,0)
NGALu_results[,11]  <- ifelse(NGALu_results[,9]<0,1,0)
colnames(NGALu_results) <- c("QALYs_base","Costs_base","QALYs_NGALu","Costs_NGALu","NB_base","NB_NGALu","INB", "Incr.QALY", "Incr.Cost", "p_eff", "p_save")
NGALu_results       <- data.frame(NGALu_results, NGALu_split_results)
#write.csv(NGALu_results, "NGALu_results_basecase.csv")

NGALs_results       <- data.frame(CEresults.life_base, CEresults.life_NGALs, NB.base, NB.NGALs, INB_NGALs)
NGALs_results[,8]   <-  NGALs_results[,3]- NGALs_results[,1]
NGALs_results[,9]   <-   NGALs_results[,4]- NGALs_results[,2]
NGALs_results[,10]  <- ifelse(NGALs_results[,8]>0,1,0)
NGALs_results[,11]  <- ifelse(NGALs_results[,9]<0,1,0)
colnames(NGALs_results) <- c("QALYs_base","Costs_base","QALYs_NGALs","Costs_NGALs","NB_base","NB_NGALs","INB", "Incr.QALY", "Incr.Cost", "p_eff", "p_save")
NGALs_results       <- data.frame(NGALs_results, NGALs_split_results)
# write.csv(NGALs_results, "NGALs_results_basecase.csv")

#==== 95% CIs ===#
# Conf_ints <- array(NA, c(length(NGALp_results[1,]),2))
# for (i in 1:11){
#   Conf_ints[i,1] <- quantile(NGALp_results[,i],c(0.025))
#   Conf_ints[i,2] <- quantile(NGALp_results[,i], c(0.975))
# }
#write.csv(Conf_ints, "Conf_ints_NGALp.csv")

# Conf_ints <- array(NA, c(length(NGALu_results[1,]),2))
# for (i in 1:11){
#   Conf_ints[i,1] <- quantile(NGALu_results[,i],c(0.025))
#   Conf_ints[i,2] <- quantile(NGALu_results[,i], c(0.975))
# }
# write.csv(Conf_ints, "Conf_ints_NGALu.csv")

# Conf_ints <- array(NA, c(length(NGALs_results[1,]),2))
# for (i in 1:11){
#   Conf_ints[i,1] <- quantile(NGALs_results[,i],c(0.025))
#   Conf_ints[i,2] <- quantile(NGALs_results[,i], c(0.975))
# }
# write.csv(Conf_ints, "Conf_ints_NGALs.csv")

#=== Traces ====#
# write.csv(trace_hosp, "trace_hosp.csv")
# write.csv(trace_FUP, "trace_FUP.csv")

# write.csv(trace_hosp_NGALp, "trace_hosp_NGALp.csv")
# write.csv(trace_NGALp_FUP, "trace_NGALp_FUP.csv")
 
#  write.csv(trace_hosp_NGALu, "trace_hosp_NGALu.csv")
#  write.csv(trace_NGALu_FUP, "trace_NGALu_FUP.csv")

#  write.csv(trace_hosp_NGALs, "trace_hosp_NGALs.csv")
#  write.csv(trace_NGALs_FUP, "trace_NGALs_FUP.csv")

#============ CystatinC VS. STANDARD CARE ===========#

CystCp_split_results  <- data.frame(CEresults.HOSP_base, CEresults.HOSP_CystCp, CEresults.FUP_base, CEresults.FUP_CystCp)
colnames(CystCp_split_results)  <- c("Hosp_Q_base", "Hosp_C_base", "Hosp_Q_CystCp", "Hosp_C_CystCp", "FUP_Q_base", "FUP_C_base", "FUP_Q_CystCp", "FUP_C_CystCp")

CystCu_split_results  <- data.frame(CEresults.HOSP_base, CEresults.HOSP_CystCu, CEresults.FUP_base, CEresults.FUP_CystCu)
colnames(CystCu_split_results)  <- c("Hosp_Q_base", "Hosp_C_base", "Hosp_Q_CystCu", "Hosp_C_CystCu", "FUP_Q_base", "FUP_C_base", "FUP_Q_CystCu", "FUP_C_CystCu")

CystCs_split_results  <- data.frame(CEresults.HOSP_base, CEresults.HOSP_CystCs, CEresults.FUP_base, CEresults.FUP_CystCs)
colnames(CystCs_split_results)  <- c("Hosp_Q_base", "Hosp_C_base", "Hosp_Q_CystCs", "Hosp_C_CystCs", "FUP_Q_base", "FUP_C_base", "FUP_Q_CystCs", "FUP_C_CystCs")
 
#=== Combine HOSP and FUP model results ===# 
CEresults.life_base    <- CEresults.HOSP_base + CEresults.FUP_base  
CEresults.life_CystCp  <- CEresults.HOSP_CystCp +CEresults.FUP_CystCp
CEresults.life_CystCu  <- CEresults.HOSP_CystCu +CEresults.FUP_CystCu
CEresults.life_CystCs  <- CEresults.HOSP_CystCs +CEresults.FUP_CystCs
 
#=== Two- way Incremental Net Health Benefit analysis ===#
threshold    <- 20000
NB.base      <- CEresults.life_base[,1] - CEresults.life_base[,2]/threshold
NB.CystCp    <- CEresults.life_CystCp[,1] - CEresults.life_CystCp[,2]/threshold
INB_CystCp   <- NB.CystCp-NB.base
 NB.CystCu   <- CEresults.life_CystCu[,1] - CEresults.life_CystCu[,2]/threshold
INB_CystCu   <- NB.CystCu-NB.base
NB.CystCs    <- CEresults.life_CystCs[,1] - CEresults.life_CystCs[,2]/threshold
INB_CystCs   <- NB.CystCs-NB.base

#=== Compile and save results ===# 
CystCp_results       <- data.frame(CEresults.life_base, CEresults.life_CystCp, NB.base, NB.CystCp, INB_CystCp)
CystCp_results[,8]   <-  CystCp_results[,3]- CystCp_results[,1]
CystCp_results[,9]   <-   CystCp_results[,4]- CystCp_results[,2]
CystCp_results[,10]  <- ifelse(CystCp_results[,8]>0,1,0)
CystCp_results[,11]  <- ifelse(CystCp_results[,9]<0,1,0)
colnames(CystCp_results) <- c("QALYs_base","Costs_base","QALYs_CystCp","Costs_CystCp","NB_base","NB_CystCp","INB", "Incr.QALY", "Incr.Cost", "p_eff", "p_save")
CystCp_results       <- data.frame(CystCp_results, CystCp_split_results)
# write.csv(CystCp_results, "CystCp_results_cardio.csv")

CystCu_results      <- data.frame(CEresults.life_base, CEresults.life_CystCu, NB.base, NB.CystCu, INB_CystCu)
CystCu_results[,8]  <-  CystCu_results[,3]- CystCu_results[,1]
CystCu_results[,9]  <-   CystCu_results[,4]- CystCu_results[,2]
CystCu_results[,10] <- ifelse(CystCu_results[,8]>0,1,0)
CystCu_results[,11] <- ifelse(CystCu_results[,9]<0,1,0)
CystCu_results      <- data.frame(CystCu_results, CystCu_split_results)
colnames(CystCu_results) <- c("QALYs_base","Costs_base","QALYs_CystCu","Costs_CystCu","NB_base","NB_CystCu","INB", "Incr.QALY", "Incr.Cost", "p_eff", "p_save")
#write.csv(CystCu_results, "CystCu_base_NEW.csv")

CystCs_results      <- data.frame(CEresults.life_base, CEresults.life_CystCs, NB.base, NB.CystCs, INB_CystCs)
CystCs_results[,8]  <-  CystCs_results[,3]- CystCs_results[,1]
CystCs_results[,9]  <-   CystCs_results[,4]- CystCs_results[,2]
CystCs_results[,10] <- ifelse(CystCs_results[,8]>0,1,0)
CystCs_results[,11] <- ifelse(CystCs_results[,9]<0,1,0)
colnames(CystCs_results) <- c("QALYs_base","Costs_base","QALYs_CystCs","Costs_CystCs","NB_base","NB_CystCs","INB", "Incr.QALY", "Incr.Cost", "p_eff", "p_save")
CystCs_results      <- data.frame(CystCs_results, CystCs_split_results) 
#write.csv(CystCs_results, "CystCs_pAKI_010.csv")


#==== 95% CIs ===#
# Conf_ints <- array(NA, c(length(CystCp_results[1,]),2))
# for (i in 1:11){
#   Conf_ints[i,1] <- quantile(CystCp_results[,i],c(0.025))
#   Conf_ints[i,2] <- quantile(CystCp_results[,i], c(0.975))
# }
# write.csv(Conf_ints, "Conf_CystCp_cardio.csv")
# 
# Conf_ints <- array(NA, c(length(CystCu_results[1,]),2))
# for (i in 1:11){
#   Conf_ints[i,1] <- quantile(CystCu_results[,i],c(0.025))
#   Conf_ints[i,2] <- quantile(CystCu_results[,i], c(0.975))
# }
# write.csv(Conf_ints, "Conf_CystCu_new.csv")
# 
# Conf_ints <- array(NA, c(length(CystCs_results[1,]),2))
# for (i in 1:11){
#   Conf_ints[i,1] <- quantile(CystCs_results[,i],c(0.025))
#   Conf_ints[i,2] <- quantile(CystCs_results[,i], c(0.975))
# }
# write.csv(Conf_ints, "Conf_CystCs_cardio.csv")

#=== Traces ====#
# write.csv(trace_hosp, "trace_hosp.csv")
# write.csv(trace_FUP, "trace_FUP.csv")
# write.csv(trace_hosp_NGALp, "trace_hosp_NGALp.csv")
# write.csv(trace_NGALp_FUP, "trace_NGALp_FUP.csv")
#  write.csv(trace_hosp_CystCu, "trace_hosp_CystCu_cardio.csv")
#  write.csv(trace_CystCu_FUP, "trace_CystCu_FUP_Cardio.csv")
#  write.csv(trace_hosp_CystCs, "trace_hosp_CystCs.csv")
#  write.csv(trace_CystCs_FUP, "trace_CystCs_FUP.csv")



#=========================================================================================================================#
#============================================ SCATTER PLOT & CEAF ================================================#
#=========================================================================================================================#

#=== Scatter Plot ===#

#Copy relevant code into plot code below
    #Nephro_results
    #NGALp_results
    #NGALu_results
    #NGALs_results
    #CystCp_results
    #CystCu_results
    #CystCs_results 

#===NEPHROCHECK===#
plot(Nephro_results[c(8,9)], xlim=c(-0.7,0.7), ylim=c(-6000, 5000),  xlab = "Incremental QALY",
     ylab = "Incremental Cost (£)", yaxs="i", col="deepskyblue3")
abline(a=0, b=20000, col = "orange", lwd = 2)
abline(a=0, b=0, col = "black", lwd = 0.5)
abline(v=0)
points(x=mean(Nephro_results$Incr.QALY), y=mean(Nephro_results$Incr.Cost),pch=24, bg="white", col="black", cex=1.2)
labels <- c("£20,000/ˆ†QALY threshold", "Mean value (0.016 QALY, £225)")
cols <- c("orange","black")
legend(0.05,-3600, labels, col = cols, cex = 0.8, lwd=2, lty=c(1,NA), pch=c(NA, 24))

ICER <- mean(Nephro_results$Incr.Cost)/mean(Nephro_results$Incr.QALY)
ICER

#mean(Nephro_results$Incr.QALY)
#mean(Nephro_results$p_eff)

##===NGAL===#
plot(NGALp_results[c(8,9)], xlim=c(-0.7,0.7), ylim=c(-6000, 5000),  xlab = "Incremental QALY",
     ylab = "Incremental Cost (£)", yaxs="i", col="deepskyblue3")
abline(a=0, b=20000, col = "orange", lwd = 2)
abline(a=0, b=0, col = "black", lwd = 0.5)
abline(v=0)
points(x=mean(NGALp_results$Incr.QALY), y=mean(NGALp_results$Incr.Cost),pch=24, bg="white", col="black", cex=1.2)
labels <- c("£20,000/QALY threshold", "Mean value (0.016 QALY, £225)")
cols <- c("orange","black")
legend(0.05,-3600, labels, col = cols, cex = 0.8, lwd=2, lty=c(1,NA), pch=c(NA, 24))


#===CYSTATIN C PLASMA===#
plot(CystCp_results[c(8,9)], xlim=c(-0.7,0.7), ylim=c(-6000, 5000),  xlab = "Incremental QALY",
     ylab = "Incremental Cost (£)", yaxs="i", col="deepskyblue3")
abline(a=0, b=20000, col = "orange", lwd = 2)
abline(a=0, b=0, col = "black", lwd = 0.5)
abline(v=0)
points(x=mean(CystCp_results$Incr.QALY), y=mean(CystCp_results$Incr.Cost),pch=24, bg="white", col="black", cex=1.2)
labels <- c("£20,000/QALY threshold", "Mean value (0.016 QALY, £225)")
cols <- c("orange","black")
legend(0.05,-3600, labels, col = cols, cex = 0.8, lwd=2, lty=c(1,NA), pch=c(NA, 24))

#===CYSTATIN C URINE===#
plot(CystCu_results[c(8,9)], xlim=c(-0.7,0.7), ylim=c(-6000, 5000),  xlab = "Incremental QALY",
     ylab = "Incremental Cost (£)", yaxs="i", col="deepskyblue3")
abline(a=0, b=20000, col = "orange", lwd = 2)
abline(a=0, b=0, col = "black", lwd = 0.5)
abline(v=0)
points(x=mean(CystCu_results$Incr.QALY), y=mean(CystCu_results$Incr.Cost),pch=24, bg="white", col="black", cex=1.2)
labels <- c("£20,000/QALY threshold", "Mean value (0.007 QALY, £73)")
cols <- c("orange","black")
legend(0.05,-3600, labels, col = cols, cex = 0.8, lwd=2, lty=c(1,NA), pch=c(NA, 24))

#===CYSTATIN C SERUM===#
plot(CystCs_results[c(8,9)], xlim=c(-0.7,0.7), ylim=c(-6000, 5000),  xlab = "Incremental QALY",
     ylab = "Incremental Cost (£)", yaxs="i", col="deepskyblue3")
abline(a=0, b=20000, col = "orange", lwd = 2)
abline(a=0, b=0, col = "black", lwd = 0.5)
abline(v=0)
points(x=mean(CystCs_results$Incr.QALY), y=mean(CystCs_results$Incr.Cost),pch=24, bg="white", col="black", cex=1.2)
labels <- c("£20,000/QALY threshold", "Mean value (0.011 QALY, £102)")
cols <- c("orange","black")
legend(0.05,-3600, labels, col = cols, cex = 0.8, lwd=2, lty=c(1,NA), pch=c(NA, 24))



#==================== CEAF  ======================#                   

# ALL
CystCs_results <<- read.csv("Cystatin C\\Base case\\Serum\\CystCs_base_new.csv")   #NEW USERS: update pathway code here as necessary to read in appropriate results files
CystCp_results <<- read.csv("Cystatin C\\Base case\\Plasma\\CystCp_base_new.csv")
CystCu_results <<- read.csv("Cystatin C\\Base case\\Urine\\CystCu_base_new.csv")
NGALs_results  <<- read.csv("NGAL\\Base case\\Serum\\NGALs_base_new.csv")
NGALp_results  <<- read.csv("NGAL\\Base case\\plasma\\NGALp_base_new.csv")
NGALu_results  <<- read.csv("NGAL\\Base case\\Urine\\NGALu_base_new.csv")
Nephro_results <<- read.csv("Nephrocheck\\Base case\\Nephro_base_new.csv")

CystCs_results <<- read.csv("Cystatin C\\Post cardiac subgroup\\CystCs_results_cardio.csv")
CystCp_results <<- read.csv("Cystatin C\\Post cardiac subgroup\\CystCp_results_cardio.csv")
CystCu_results <<- read.csv("Cystatin C\\Post cardiac subgroup\\CystCu_results_cardio.csv")
NGALs_results  <<- read.csv("NGAL\\Post cardiac subgroup\\NGALs_results_cardio.csv")
NGALp_results  <<- read.csv("NGAL\\Post cardiac subgroup\\NGALp_results_cardio.csv")
NGALu_results  <<- read.csv("NGAL\\Post cardiac subgroup\\NGALu_results_cardio.csv")
Nephro_results <<- read.csv("Nephrocheck\\Post cardiac subgroup\\Nephro_cardio_new.csv")


OUT <- list()
thresh <- seq(100,200000,100)

for (i in 1:length(thresh)) {

  NB.Nephro                 <- Nephro_results$QALYs_Nephro - Nephro_results$Costs_Nephro/thresh[i]  #Note NB in QALYS (i.e. net health benefit, NOT monetary)
  NB.NGALp                  <- NGALp_results$QALYs_NGALp - NGALp_results$Costs_NGALp/thresh[i]  #Note NB in QALYS (i.e. net health benefit, NOT monetary)
  NB.NGALu                  <- NGALu_results$QALYs_NGALu - NGALu_results$Costs_NGALu/thresh[i]  #Note NB in QALYS (i.e. net health benefit, NOT monetary)
  NB.NGALs                  <- NGALs_results$QALYs_NGALs - NGALs_results$Costs_NGALs/thresh[i]  #Note NB in QALYS (i.e. net health benefit, NOT monetary)
   NB.CystCp                <- CystCp_results$QALYs_CystCp - CystCp_results$Costs_CystCp/thresh[i]  #Note NB in QALYS (i.e. net health benefit, NOT monetary)
   NB.CystCu                <- CystCu_results$QALYs_CystCu - CystCu_results$Costs_CystCu/thresh[i]  #Note NB in QALYS (i.e. net health benefit, NOT monetary)
   NB.CystCs                <- CystCs_results$QALYs_CystCs - CystCs_results$Costs_CystCs/thresh[i]  #Note NB in QALYS (i.e. net health benefit, NOT monetary)
   NB.base                  <- Nephro_results$QALYs_base - Nephro_results$Costs_base/thresh[i]
  maxNB                    <- ifelse(NB.base>= NB.Nephro, NB.base, NB.Nephro)
  maxNB                    <- ifelse(maxNB>= NB.NGALp, maxNB, NB.NGALp)
  maxNB                    <- ifelse(maxNB>= NB.NGALu, maxNB, NB.NGALu)
  maxNB                    <- ifelse(maxNB>= NB.NGALs, maxNB, NB.NGALs)
  maxNB                    <- ifelse(maxNB>= NB.CystCp, maxNB, NB.CystCp)
  maxNB                    <- ifelse(maxNB>= NB.CystCu, maxNB, NB.CystCu)
  maxNB                    <- ifelse(maxNB>= NB.CystCs, maxNB, NB.CystCs)
  CE.base                  <- ifelse(NB.base == maxNB, 1, 0)
  CE.Nephro                  <- ifelse(NB.Nephro == maxNB, 1, 0)
  CE.NGALp                  <- ifelse(NB.NGALp == maxNB, 1, 0)
  CE.NGALu                  <- ifelse(NB.NGALu == maxNB, 1, 0)
  CE.NGALs                  <- ifelse(NB.NGALs == maxNB, 1, 0)
  CE.CystCp                <- ifelse(NB.CystCp == maxNB, 1, 0)
  CE.CystCs                  <- ifelse(NB.CystCs == maxNB, 1, 0)
  CE.CystCu                  <- ifelse(NB.CystCu == maxNB, 1, 0)
  prob.CE.base             <- mean(CE.base)
  prob.CE.Nephro             <- mean(CE.Nephro)
  prob.CE.NGALp             <- mean(CE.NGALp)
  prob.CE.NGALs             <- mean(CE.NGALs)
  prob.CE.NGALu             <- mean(CE.NGALu)
  prob.CE.CystCp             <- mean(CE.CystCp )
  prob.CE.CystCs             <- mean(CE.CystCs )
  prob.CE.CystCu             <- mean(CE.CystCu )
  mean.NB.base             <- mean(NB.base)
  mean.NB.Nephro             <- mean(NB.Nephro)
  mean.NB.NGALp             <- mean(NB.NGALp)
  mean.NB.NGALs             <- mean(NB.NGALs)
  mean.NB.NGALu             <- mean(NB.NGALu)
  mean.NB.CystCp             <- mean(NB.CystCp )
  mean.NB.CystCs              <- mean(NB.CystCs )
  mean.NB.CystCu              <- mean(NB.CystCu )
  
  mean.NB.max              <- max(mean.NB.base,mean.NB.Nephro,mean.NB.NGALp,mean.NB.NGALs,mean.NB.NGALu ,mean.NB.CystCp ,mean.NB.CystCs ,mean.NB.CystCu )
#  ceaf                     <- ifelse(mean.NB.max==mean.NB.base, prob.CE.base,prob.CE.test)
  ceaf                     <- ifelse(mean.NB.max==mean.NB.base, prob.CE.base,0)
  ceaf                     <- ifelse(mean.NB.max== mean.NB.Nephro,  mean.NB.Nephro,ceaf )
  ceaf                     <- ifelse(mean.NB.max==mean.NB.NGALp,prob.CE.NGALp,ceaf )
  ceaf                     <- ifelse(mean.NB.max==mean.NB.NGALs, prob.CE.NGALs,ceaf )
  ceaf                     <- ifelse(mean.NB.max==mean.NB.NGALu, prob.CE.NGALu,ceaf )
  ceaf                     <- ifelse(mean.NB.max== mean.NB.CystCp, prob.CE.CystCp,ceaf )
  ceaf                     <- ifelse(mean.NB.max== mean.NB.CystCs, prob.CE.CystCs,ceaf )
  ceaf                     <- ifelse(mean.NB.max== mean.NB.CystCu, prob.CE.CystCu,ceaf )
  ceaf.test                     <- ifelse(mean.NB.max==mean.NB.Nephro, 1,0)
  ceaf.test                     <- ifelse(mean.NB.max==mean.NB.NGALp,21,0)
  ceaf.test                     <- ifelse(mean.NB.max==mean.NB.NGALs, 23,0)
  ceaf.test                     <- ifelse(mean.NB.max==mean.NB.NGALu, 22,0)
  ceaf.test                     <- ifelse(mean.NB.max== mean.NB.CystCp, 31,0)
  ceaf.test                     <- ifelse(mean.NB.max== mean.NB.CystCs, 33,0)
  ceaf.test                     <- ifelse(mean.NB.max== mean.NB.CystCu, 32,0)

  OUT[[i]] <- list(NB.base=NB.base,NB.Nephro=NB.Nephro,maxNB=maxNB,
                   NB.NGALp=      NB.NGALp,           
                   NB.NGALu =      NB.NGALu  ,         
                   NB.NGALs  =     NB.NGALs,           
                   NB.CystCp  =       NB.CystCp  ,       
                   NB.CystCu  =           NB.CystCu  ,   
                   NB.CystCs  =             NB.CystCs ,  
                   
                   CE.base=CE.base,
                   CE.Nephro=CE.Nephro,
                   CE.NGALp=CE.NGALp,
                   CE.NGALs=CE.NGALs,
                   CE.NGALu=CE.NGALu,
                   CE.CystCp=CE.CystCp,
                   CE.CystCs=CE.CystCs,
                   CE.CystCu=CE.CystCu,
                   prob.CE.base=prob.CE.base, prob.CE.Nephro=prob.CE.Nephro,
                   prob.CE.NGALp=prob.CE.NGALp,
                   prob.CE.NGALs=prob.CE.NGALs,
                   prob.CE.NGALu=prob.CE.NGALu,
                   prob.CE.CystCp=prob.CE.CystCp,
                   prob.CE.CystCs=prob.CE.CystCs,
                   prob.CE.CystCu=prob.CE.CystCu,
                   mean.NB.base = mean.NB.base, mean.NB.Nephro= mean.NB.Nephro,
                   mean.NB.NGALp=mean.NB.NGALp,
                   mean.NB.NGALs=mean.NB.NGALs,
                   mean.NB.NGALu=mean.NB.NGALu,
                   mean.NB.CystCp=mean.NB.CystCp,
                   mean.NB.CystCs=mean.NB.CystCs,
                   mean.NB.CystCu=mean.NB.CystCu,
                   ceaf=ceaf, 
                   ceaf.test=ceaf.test)
}

CEAF <- data.frame(threshold=thresh,Probability.CE = NA,test=NA)
for (i in 1:length(thresh)){
  CEAF$Probability.CE[i]  <- OUT[[i]]$ceaf
  CEAF$CE.base[i]         <- OUT[[i]]$prob.CE.base
  CEAF$CE.Neprho[i]       <- OUT[[i]]$prob.CE.Nephro
  CEAF$CE.NGALp[i]       <- OUT[[i]]$prob.CE.NGALp
  CEAF$CE.NGALs[i]       <- OUT[[i]]$prob.CE.NGALs
  CEAF$CE.NGALu[i]       <- OUT[[i]]$prob.CE.NGALu
  CEAF$CE.CystCp[i]       <- OUT[[i]]$prob.CE.CystCp
  CEAF$CE.CystCs[i]       <- OUT[[i]]$prob.CE.CystCs
  CEAF$CE.CystCu[i]       <- OUT[[i]]$prob.CE.CystCu
  CEAF$Probability.CE[i] <- OUT[[i]]$ceaf
  
  }


plot(CEAF$threshold,CEAF$CE.base,type="l", ylim = c(0,1), xlim = c(0,50000), xaxp  = c(0, 50000, 5), col = "deepskyblue3", lwd=3, xlab = "Willingness-to-pay threshold (£ per QALY)", ylab = "Probability cost-effective")
lines(CEAF$threshold,CEAF$CE.Nephro,type="l", ylim = c(0,1), xlim = c(0,50000), col = "black",lwd=2)
lines(CEAF$threshold,CEAF$CE.NGALp,type="l", ylim = c(0,1), xlim = c(0,50000), col = "aquamarine",lwd=2)  ###
lines(CEAF$threshold,CEAF$CE.NGALs,type="l", ylim = c(0,1), xlim = c(0,50000), col = "orange",lwd=2)  #COST-EFF
lines(CEAF$threshold,CEAF$CE.NGALu,type="l", ylim = c(0,1), xlim = c(0,50000), col = "deeppink3",lwd=2)
lines(CEAF$threshold,CEAF$CE.CystCp,type="l", ylim = c(0,1), xlim = c(0,50000), col = "darkmagenta",lwd=2)
lines(CEAF$threshold,CEAF$CE.CystCs,type="l", ylim = c(0,1), xlim = c(0,50000), col = "cyan3",lwd=2)  #COST-EFF 
lines(CEAF$threshold,CEAF$CE.CystCu,type="l", ylim = c(0,1), xlim = c(0,50000), col = "grey",lwd=2)     #
lines(CEAF$threshold[seq(1,500,6)],CEAF$Probability.CE[seq(1,500,6)],type="p", ylim = c(0,1))
labels <- c("Standard Care","Nephrocheck", "NGAL plasma", "NGAL serum", "NGAL urine", "Cystatin-C plasma", "Cystatin-C serum", "Cystatin-C urine", "Cost-effectiveness frontier")
cols <- c("deepskyblue3", "black", "aquamarine", "orange", "deeppink3", "darkmagenta", "cyan3", "grey")
legend(30000,1,labels, col = cols, lty = c(1,1,1,1,1,1,1,1,NA), pch=c(NA, NA,NA,NA,NA,NA,NA,NA, 1), cex = 1, lwd=c(2,2,2,2,2,2,2,2,1))
#write.csv(CEAF, "CEAF_ALL_cardio.csv", row.names=FALSE)

##NB: to get CEAF data into excel columns, copy over csv data to excel, then go to data tab -> text to Columns-> follow instructions to get data converted in to table columns


#===== NGAL =====#
OUT <- list()
thresh <- seq(100,200000,100)

for (i in 1:length(thresh)) {
  
  NB.test                  <- NGALp_results$QALYs_NGALp - NGALp_results$Costs_NGALp/thresh[i]  #Note NB in QALYS (i.e. net health benefit, NOT monetary)
  NB.base                  <- NGALp_results$QALYs_base - NGALp_results$Costs_base/thresh[i]  
  maxNB                    <- ifelse(NB.base>= NB.test, NB.base, NB.test)   
  CE.base                  <- ifelse(NB.base == maxNB, 1, 0)
  CE.test                  <- ifelse(NB.test == maxNB, 1, 0)
  prob.CE.base             <- mean(CE.base) 
  prob.CE.test             <- mean(CE.test) 
  mean.NB.base             <- mean(NB.base)
  mean.NB.test             <- mean(NB.test)
  mean.NB.max              <- max(mean.NB.base,mean.NB.test)    
  ceaf                     <- ifelse(mean.NB.max==mean.NB.base, prob.CE.base,prob.CE.test)
  
  OUT[[i]] <- list(NB.base=NB.base,NB.test=NB.test,maxNB=maxNB,
                   CE.base=CE.base,CE.test=CE.test, prob.CE.base=prob.CE.base, prob.CE.test=prob.CE.test,
                   mean.NB.base = mean.NB.base, mean.NB.test= mean.NB.test,ceaf=ceaf)
}

CEAF <- data.frame(threshold=thresh)
for (i in 1:length(thresh)){
  CEAF$Probability.CE[i]  <- OUT[[i]]$ceaf
  CEAF$CE.base[i]         <- OUT[[i]]$prob.CE.base
  CEAF$CE.test[i]       <- OUT[[i]]$prob.CE.test
}


plot(CEAF$threshold,CEAF$CE.base,type="l", ylim = c(0,1), xlim = c(0,50000), xaxp  = c(0, 50000, 5), col = "deepskyblue3", lwd=2, xlab = "Willingness-to-pay threshold (£ per QALY)", ylab = "Probability cost-effective")
lines(CEAF$threshold,CEAF$CE.test,type="l", ylim = c(0,1), xlim = c(0,50000), col = "orange",lwd=2)
lines(CEAF$threshold[seq(1,500,6)],CEAF$Probability.CE[seq(1,500,6)],type="p", ylim = c(0,1))
labels <- c("Standard Care","NGAL (plasma)", "Cost-effectiveness frontier")
cols <- c("deepskyblue3","orange", "black")
legend(20000,0.30,labels, col = cols, lty = c(1,1,NA), pch=c(NA, NA, 1), cex = 1, lwd=c(2,2,1))
write.csv(CEAF, "CEAF_NGALp.csv", row.names=FALSE)




#CYSTATIN C
OUT <- list()
thresh <- seq(100,200000,100)

for (i in 1:length(thresh)) {

  NB.test                  <- CystCp_results$QALYs_CystCp - CystCp_results$Costs_CystCp/thresh[i]  #Note NB in QALYS (i.e. net health benefit, NOT monetary)
  NB.base                  <- CystCp_results$QALYs_base - CystCp_results$Costs_base/thresh[i]
  maxNB                    <- ifelse(NB.base>= NB.test, NB.base, NB.test)
  CE.base                  <- ifelse(NB.base == maxNB, 1, 0)
  CE.test                  <- ifelse(NB.test == maxNB, 1, 0)
  prob.CE.base             <- mean(CE.base)
  prob.CE.test             <- mean(CE.test)
  mean.NB.base             <- mean(NB.base)
  mean.NB.test             <- mean(NB.test)
  mean.NB.max              <- max(mean.NB.base,mean.NB.test)
  ceaf                     <- ifelse(mean.NB.max==mean.NB.base, prob.CE.base,prob.CE.test)

  OUT[[i]] <- list(NB.base=NB.base,NB.test=NB.test,maxNB=maxNB,
                   CE.base=CE.base,CE.test=CE.test, prob.CE.base=prob.CE.base, prob.CE.test=prob.CE.test,
                   mean.NB.base = mean.NB.base, mean.NB.test= mean.NB.test,ceaf=ceaf)
}

CEAF <- data.frame(threshold=thresh)
for (i in 1:length(thresh)){
  CEAF$Probability.CE[i]  <- OUT[[i]]$ceaf
  CEAF$CE.base[i]         <- OUT[[i]]$prob.CE.base
  CEAF$CE.test[i]       <- OUT[[i]]$prob.CE.test
}


plot(CEAF$threshold,CEAF$CE.base,type="l", ylim = c(0,1), xlim = c(0,50000), xaxp  = c(0, 50000, 5), col = "deepskyblue3", lwd=2, xlab = "Willingness-to-pay threshold (£ per QALY)", ylab = "Probability cost-effective")
lines(CEAF$threshold,CEAF$CE.test,type="l", ylim = c(0,1), xlim = c(0,50000), col = "orange",lwd=2)
lines(CEAF$threshold[seq(1,500,6)],CEAF$Probability.CE[seq(1,500,6)],type="p", ylim = c(0,1))
labels <- c("Standard Care","Cystatin C (urine)", "Cost-effectiveness frontier")
cols <- c("deepskyblue3","orange", "black")
legend(20000,0.30,labels, col = cols, lty = c(1,1,NA), pch=c(NA, NA, 1), cex = 1, lwd=c(2,2,1))
write.csv(CEAF, "CEAF_CystCu_NEW.csv", row.names=FALSE)


#================================================================================================================================#
#================================= VALUE OF INFORMATION: NON PARAMETRIC REGRESSION ANALYSIS =====================================#
#================================================================================================================================#

#Load packages 
#install.packages("grid")#, repos="https://www.math.ntnu.no/inla/R/stable")  #not available for R version 3.2.1
#install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable")  
#install.packages("earth")
library(expm)
library(devtools)
library(ldr)
library(splancs)
library(INLA)
library(grid)
library(ggplot2)
library(BCEA)
library(mgcv)
library(earth)


#===================== BCEA PACKAGE for PSA and EVPI outputs =====================#

#create effects (e) and cost (c) data frames
e <- data.frame(NGALs_results$QALYs_base,NGALs_results$QALYs_NGALs) 
c <- data.frame(NGALs_results$Costs_base,NGALs_results$Costs_NGALs) 
e <- data.matrix(e, rownames.force = NA)  #must be numeric data, so force into numeric 
c <- data.matrix(c, rownames.force = NA)
#dim(e)
#dim(CEresults.life_base)
#is.numeric(c)

#Run BCEA
m <- bcea(e,c, ref = 2, interventions = NULL)

#Outputs & plots
# summary(m)
# plot(m)
# ceplane.plot(m)  
# eib.plot(m, plot.cri=FALSE)
# contour(m)
# m$ib <- m$ib*20000
ib.plot(m, wtp=20000)
# ceac.plot(m)
# 
# mce <- multi.ce(m)
# ceaf.plot(mce)
# mce.plot(mce)
# 
# #ceef.plot(mce)
# 
# m$evi  <- m$evi/20000
#evi.plot(m) 

plot(m$k,m$evi, xlab="Willingness to pay (£/QALY)", ylab="Per patient EVPI (£)", pch=20, cex=0.5)     

m$evi[m$k==20000]


#=================================== EVPPI: single parameters ============================#

library(mgcv)

#Get data
setwd("N:\\Faculty-of-Medicine-and-Health\\LIHS\\LIHS-VOL1\\HealthEconomics\\PROJECTS\\ON GOING PROJECTS_FUndED\\AKI-Diagnostics\\Economic Model\\R Model\\Model Results")

CystCs_results <<- read.csv("Cystatin C\\Base case\\Serum\\CystCs_results_base.csv")
theta_HOSP <<- read.csv("Cystatin C\\theta_HOSP_CystCs.csv")
theta_FUP <<- read.csv("Cystatin C\\theta_FUP_CystCs.csv")

NGALs_results <<- read.csv("NGAL\\Base case\\Serum\\NGALs_results_basecase.csv")
theta_HOSP <<- read.csv("NGAL\\theta_HOSP_NGALs.csv")
theta_FUP <<- read.csv("NGAL\\theta_FUP_NGALs.csv")

Nephro_results <<- read.csv("Nephrocheck\\Base case\\Nephro_base_newdisc.csv")
theta_HOSP <<- read.csv("Nephrocheck\\VOI\\theta_HOSP.csv")
theta_FUP <<- read.csv("Nephrocheck\\VOI\\theta_FUP.csv")

theta_HOSP <- theta_HOSP[-1]
theta_FUP <- theta_FUP[-1]
theta <- data.frame(theta_HOSP, theta_FUP)
#is.numeric(theta_HOSP)
theta_HOSP <- data.matrix(theta_HOSP, rownames.force = NA)
theta_FUP <- data.matrix(theta_FUP, rownames.force = NA)
theta <- data.matrix(theta, rownames.force = NA)
theta <- theta[,-1]

#Get INB and matrix of theta's 
INB <- CystCs_results$INB*20000 #get in monetary terms
INB <- NGALs_results$INB*20000 
INB <- Nephro_results$INB*20000
#is.numeric(INB)


#=== Using a loop ===#
# initialise a vector to hold the results.
evppi_vector <- numeric(length(theta[1,]))  #Number of parameters (87 not including large matrices)

# Run loop
for (i in (1:length(evppi_vector))){  #[-c(7, 16)]) { # theta 7 and 16 are not defined
  parameter_of_interest <- theta[1:Nsim, i]
  model <- gam(INB[1:Nsim] ~ te(parameter_of_interest))
  fittedValues <- fitted(model)
  evppi_vector[i] <- mean(pmax(0, fittedValues)) - max(0, mean(fittedValues))
}
names(evppi_vector) <- colnames(theta)

#Sort into order
sort(evppi_vector, decreasing=TRUE)

#Save only the non-zero values
evppi_vector_top <- evppi_vector[evppi_vector!=0]
print(sort(evppi_vector_top,decreasing=TRUE),1) 

#Save EVPPI
#write.csv(evppi_vector_top, "evppi_perperson_INMB.csv", row.names=TRUE)

#population values   
#  I <- NA # incident population over 10 years
#  for (t in 1:10) {
#    I[t] <- 258956/(1.035^(t-1))
#  }
# evppi_vector_top <- evppi_vector_top*sum(I)   #in millions
# 
# write.csv(evppi_vector_top, "POP_evppi_258956_10yrs_INMB.csv")

#barchart
# x.points <- barplot(evppi_vector_top/1000000, ylab = "Population EVPPI (£millions)",
#                     xaxt = "n", ylim = c(0, 400), col="deepskyblue3")
# axis(side = 1, at = x.points, labels = names(evppi_vector_top),
#      cex.axis=0.7, tick = FALSE, hadj = 0.8, las = 3)


# SINGLE PARAMETER AT A TIME
##Run regression on INB
# model <- gam(INB ~ te(theta[,1]))   #change the parameter 
# 
# #Check residuals for any structure indicating potential bias
# plot(fitted(model)[1:Nsim], residuals(model)[1:Nsim])
# 
# #Calculate EVPPI
# fittedValues <- fitted(model)
# evppi_theta <- mean(pmax(0, fittedValues)) - max(0, mean(fittedValues))
# evppi_theta

#=================================== EVPPI: multiple parameters ============================#

#==== Using GAM model ====#
#use for up to groups of 5

model <- gam(INB[1:Nsim] ~ te(theta_FUP[1:Nsim,    25:28]))  
fittedValues <- fitted(model)
evppi_group <- mean(pmax(0, fittedValues)) - max(0, mean(fittedValues))


#========== GAM/earth package ==========#
#select groups of parameters using either c(1,2,3) or 1:3
#use 'earth' package to enable groupings of >5 parameters (gam does not cope well with groups of parameters)
model <- earth(theta_FUP[,   29:38], INB, degree = 10,
                           nk = 1000, linpreds = TRUE, thres = 1e-4, Use.beta.cache=FALSE)

fittedValues <- fitted(model)
evppi_group1 <- mean(pmax(0, fittedValues)) - max(0, mean(fittedValues))
 
#plot residuals to check for any structure indicating bias 
plot(fitted(model)[1:1000], residuals(model)[1:1000])

# AKI_and_stages <- evppi_group1        #hosp 1:10 - arrive with AKI/ get AKI/ starting stages of AKI  EARTH
# RR_AKI_vs_NoAKI  <- evppi_group       #hsop 11:12 mort and discharge RRs
# ward_mort_and_disch  <- evppi_group   #hosp 14:15
# Hospital_costs  <- evppi_group        #hosp 22:25
# Disch_costs  <- evppi_group           #hosp 29:30
# Hospital_utilities <- evppi_group     #hosp 17:18
# Disch_utilities  <- evppi_group       #hosp 19,21
# NGALs_accuracy  <- evppi_group        #hosp 34:35
# NGALs_all  <- evppi_group             #hosp 32:35 (treatment cost, treatment impact, test accuracy)
# Early_treat <- evppi_group            #hosp 32:33 (early treatment impact and cost)

#Followup_start <- evppi_group1         #fup 1:5  EARTH
#Fup_5yr_mort  <- evppi_group           #fup 6:7
#Fup_risk_CKD  <- evppi_group           #fup 8:9  (RR_CKD and tp_fuptoCKD) 
#Fup_postCKD_risks  <- evppi_group1     #fup 10:21  EARTH
#Fup_utility  <- evppi_group            #fup 22:24
#Fup_CKD_utilities  <- evppi_group      #fup 25:28
#Followup_costs  <- evppi_group1        #fup 29:38 EARTH
#Fup_CKD_costs  <- evppi_group1         #fup 39:43 EARTH

#Hosp_daily_mort_and_disch  <- evppi_group   #HOSP_medium ALL (270 variables)

EVPPI_groups <- t(data.frame(#HOSP_dialy_AKI,Hosp_daily_mort_and_disch, 
                  AKI_and_stages, RR_AKI_vs_NoAKI, ward_mort_and_disch, Hospital_costs,Disch_costs,
                  Hospital_utilities,Disch_utilities, NGALs_accuracy,NGALs_all,Early_treat,Followup_start,
                  Fup_5yr_mort,Fup_risk_CKD, Fup_postCKD_risks, Followup_costs,Fup_utility, 
                  Fup_CKD_utilities, Fup_CKD_costs))

#write.csv(EVPPI_groups, "EVPPI_groups_NGALs.csv")

x.points <- barplot(EVPPI_groups[,1], ylab = "partial EVPI (per patient QALY)",
                    xaxt = "n", ylim = c(0, max(EVPPI_groups)), col="deepskyblue3")
axis(side = 1, at = x.points, labels = rownames(EVPPI_groups),
     cex.axis=0.7, tick = FALSE, hadj = 0.8, las = 3)






