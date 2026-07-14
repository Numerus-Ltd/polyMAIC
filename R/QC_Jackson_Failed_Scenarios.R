# Numerus Ltd. 
# Purpose : Performs Jackson weighting for pseudo-random selection of 5 scenarios/simulations
#           to validate production calculations 
# Notes   : Focuses on 2 scenarios where it does work and 3 where it does
#           not (i.e. no estimable IPWs, imputed instead with 1)


# PROGRAM SET-UP #######################################################################

## Load libraries 
if(!require(dplyr))    {install.packages("dplyr")    ; library(dplyr)}
if(!require(sandwich)) {install.packages("sandwich") ; library(sandwich)}
if(!require(limSolve)) {install.packages("limSolve") ; library(limSolve)}
if(!require(haven))    {install.packages("haven")    ; library(haven)}
if(!require(Hmisc))    {install.packages("Hmisc")    ; library(Hmisc)}

## Start with clean program  
rm(list=ls())
options(stringsAsFactors = FALSE)
#turn off scientific notation for all variables (mainly regarding output as csv)
options(scipen= 999)
# set the seed for reproducibility 
set.seed(872634)

## load Jackson functions
source("path/QC_Jackson_MAIC.R")

# READ IN DATA #########################################################################

# IPD and targets from SAS files (all 2000 scenarios)
data.full <- read_sas("path/data.sas7bdat")

# diagnostics from production runs (to be used to identify scenarios that did not produce a solution, weights imputed with 1)
diag.full <- read.csv("path/FINAL_DIAGNOSTICSALL.csv")

# production weights
wgt.full <- read.csv("path/FINAL_WEIGHTSALL.csv")


# SUBSET DATA ##########################################################################

# Identify scenarios where Jackson could not produce a solution: just take simulations 15 (N=50), 247 (N=100),
# and 320 (N=250). NB: Ns per treatment arm
diag.jack.no <- diag.full[diag.full$METHOD2=="JACKMAIC" & diag.full$ONEFL=="Y",]

# Identify scenarios where Jackson did produce a solution: just take simulations 1439 (N=50) and 1 (N=250).
# NB: Ns per treatment arm
diag.jack.yes <- diag.full[diag.full$METHOD2=="JACKMAIC" & diag.full$ONEFL!="Y",]

# columns in data to keep
vars <- c("Scenario", "SUBJID", "TREAT", "Xdata1", "Xdata2", "Xdata3", "Xdata4", "Xdata5", "TREATT", 
          "TMT1", "TST1", "TMT2", "TST2", "TMT3", "TST3", "TMT4", "TST4", "TMT5", "TST5",
          "TMP1", "TSP1", "TMP2", "TSP2", "TMP3", "TSP3", "TMP4", "TSP4", "TMP5", "TSP5")

# extracting the IPD and targets for the above 5 chosen scenarios (save in separate data frames)
data.15   <- data.full[data.full$Scenario==15, vars] 
data.247  <- data.full[data.full$Scenario==247, vars] 
data.320  <- data.full[data.full$Scenario==320, vars] 
data.1439 <- data.full[data.full$Scenario==1439, vars] 
data.1    <- data.full[data.full$Scenario==1, vars] 

# check distribution of TREAT and TREATT
#table(data.15$TREAT); table(data.15$TREATT)

# subset wgt.full to selected scenarios (Jackson weights only)
wgt.subset      <- wgt.full[wgt.full$SIM %in% c(15,247,320,1439,1),c(1,2,6)]

# sort by subject within simulation
wgt.subset.sort <- wgt.subset[with(wgt.subset, order(SIM, SUBJID)),] 

# check number of subjects per scenario
table(wgt.subset.sort$SIM)

# clean up workspace
rm(wgt.full, wgt.subset, diag.full, data.full, vars)

# ESTIMATE WEIGHTS #####################################################################

### Scenario 1439 (separate weights per treatment). Only a single continuous covariate (Mean [SD]) ###

## Active treatment ##
IPD.1439.1      <- data.1439[data.1439$TREAT==1,1:8]  # IPD
means.1439.1    <- as.numeric(data.1439[1,c("TMT1")]) # Target means
SDs.1439.1      <- as.numeric(data.1439[1,c("TST1")]) # Target SDs
MAIC.1439.1     <- MAIC(the_data   = IPD.1439.1,      # estimate weights. NB: estimates both SIGN (weight.maic) and JACK (weight.opt)
                        id         = "SUBJID", 
                        match_mean = "Xdata1", 
                        agg_means  = means.1439.1, 
                        agg_sds    = SDs.1439.1, 
                        accuracy   = 0.000000000001)
the_data.1439.1 <- MAIC.1439.1$the_data               # extract the data frame from the list object

#hist(the_data.1439.1$weight.opt)
#dev.off()

## Control treatment ##
IPD.1439.0      <- data.1439[data.1439$TREAT==0,1:8]  # IPD
means.1439.0    <- as.numeric(data.1439[1,c("TMP1")]) # Target means
SDs.1439.0      <- as.numeric(data.1439[1,c("TSP1")]) # Target SDs
MAIC.1439.0     <- MAIC(the_data   = IPD.1439.0,      # estimate weights. NB: estimates both SIGN (weight.maic) and JACK (weight.opt)
                        id         = "SUBJID", 
                        match_mean = "Xdata1", 
                        agg_means  = means.1439.0, 
                        agg_sds    = SDs.1439.0, 
                        accuracy   = 0.000000000001)
the_data.1439.0 <- MAIC.1439.0$the_data               # extract the data frame from the list object

#hist(the_data.1439.0$weight.opt)
#dev.off()

## Check the matching (Xdata1, mean and standard deviation), separately by treatment ##
cat("\n ----- SIMULATION 1439 -----")

# Active treatment
cat("\n --- Active ---")
cat("\n Unweighted Summary Stats \n")
data.frame("Xdata1.mean"=mean(the_data.1439.1$Xdata1), "Xdata1.SD"=sd(the_data.1439.1$Xdata1))
cat("\n Weighted Summary Stats \n")
data.frame("ESS"=ESS(the_data.1439.1$weight.opt))
data.frame("Xdata1.mean"=wtd.mean(the_data.1439.1$Xdata1, the_data.1439.1$weight.opt), 
           "Xdata1.SD"  =sqrt(wtd.var(the_data.1439.1$Xdata1, the_data.1439.1$weight.opt, normwt=TRUE)))
cat("\n Target Summary Stats \n")
data.frame("Xdata1.mean"=means.1439.1, "Xdata1.SD"=SDs.1439.1)

# control treatment
cat("\n --- Control ---")
cat("\n Unweighted Summary Stats \n")
data.frame("Xdata1.mean"=mean(the_data.1439.0$Xdata1), "Xdata1.SD"=sd(the_data.1439.0$Xdata1))
cat("\n Weighted Summary Stats \n")
data.frame("ESS"=ESS(the_data.1439.0$weight.opt))
data.frame("Xdata1.mean"=wtd.mean(the_data.1439.0$Xdata1, the_data.1439.0$weight.opt), 
           "Xdata1.SD"  =sqrt(wtd.var(the_data.1439.0$Xdata1, the_data.1439.0$weight.opt, normwt=TRUE)))
cat("\n Target Summary Stats \n")
data.frame("Xdata1.mean"=means.1439.0, "Xdata1.SD"=SDs.1439.0)

# clean the workspace
rm(data.1439, IPD.1439.1, IPD.1439.0, MAIC.1439.1, MAIC.1439.0, means.1439.1, means.1439.0, SDs.1439.1, SDs.1439.0)


#####################################################################################


### Scenario 1 (separate weights per treatment). 4 continuous covariates (Mean [SD]) and one binary (Mean=prop) [xdata4] ###

## Active treatment ##
IPD.1.1      <- data.1[data.1$TREAT==1,1:8]                                     # IPD
means.1.1    <- as.numeric(data.1[1,c("TMT1", "TMT2", "TMT3", "TMT4", "TMT5")]) # Target means
SDs.1.1      <- as.numeric(data.1[1,c("TST1", "TST2", "TST3", "TST4", "TST5")]) # Target SDs
SDs.1.1[4]   <- NA                                                              # Xdata4 is binary, so set TST4 to NA (i.e. no match on SD)
MAIC.1.1     <- MAIC(the_data   = IPD.1.1,                                      # estimate weights. NB: estimates both SIGN (weight.maic) and JACK (weight.opt)
                     id         = "SUBJID", 
                     match_mean = c("Xdata1", "Xdata2", "Xdata3", "Xdata4", "Xdata5"), 
                     agg_means  = means.1.1, 
                     agg_sds    = SDs.1.1, 
                     accuracy   = 0.000000000001)
the_data.1.1 <- MAIC.1.1$the_data                                               # extract the data frame from the list object

#hist(the_data.1.1$weight.opt)
#dev.off()

## Control treatment ##
IPD.1.0      <- data.1[data.1$TREAT==0,1:8]  # IPD
means.1.0    <- as.numeric(data.1[1,c("TMP1", "TMP2", "TMP3", "TMP4", "TMP5")]) # Target means
SDs.1.0      <- as.numeric(data.1[1,c("TSP1", "TSP2", "TSP3", "TSP4", "TSP5")]) # Target SDs
SDs.1.0[4]   <- NA                                                              # Xdata4 is binary, so set TST4 to NA (i.e. no match on SD)
MAIC.1.0     <- MAIC(the_data   = IPD.1.0,                                      # estimate weights. NB: estimates both SIGN (weight.maic) and JACK (weight.opt)
                     id         = "SUBJID", 
                     match_mean = c("Xdata1", "Xdata2", "Xdata3", "Xdata4", "Xdata5"), 
                     agg_means  = means.1.0, 
                     agg_sds    = SDs.1.0, 
                     accuracy   = 0.000000000001)
the_data.1.0 <- MAIC.1.0$the_data                                               # extract the data frame from the list object

#hist(the_data.1.0$weight.opt)
#dev.off()

## Check the matching (Xdata1/2/3/5, mean and standard deviation + Xdata4, mean), separately by treatment ##
cat("\n ----- SIMULATION 1 -----")

# Active treatment
cat("\n --- Active ---")
cat("\n Unweighted Summary Stats \n")
data.frame("Xdata1.mean"=mean(the_data.1.1$Xdata1), "Xdata1.SD"=sd(the_data.1.1$Xdata1),
           "Xdata2.mean"=mean(the_data.1.1$Xdata2), "Xdata2.SD"=sd(the_data.1.1$Xdata2),
           "Xdata3.mean"=mean(the_data.1.1$Xdata3), "Xdata3.SD"=sd(the_data.1.1$Xdata3),
           "Xdata4.mean"=mean(the_data.1.1$Xdata4), "Xdata4.SD"=NA,                      # 07 JUL 2026: correction to assigned column name (Xdata4.mean vs Xdata3.mean. No effect on result)
           "Xdata5.mean"=mean(the_data.1.1$Xdata5), "Xdata5.SD"=sd(the_data.1.1$Xdata5)) # 07 JUL 2026: correction to assigned column name (Xdata5.mean vs Xdata4.mean. No effect on result)
cat("\n Weighted Summary Stats \n")
data.frame("ESS"=ESS(the_data.1.1$weight.opt))
data.frame("Xdata1.mean"=wtd.mean(the_data.1.1$Xdata1, the_data.1.1$weight.opt), 
           "Xdata1.SD"  =sqrt(wtd.var(the_data.1.1$Xdata1, the_data.1.1$weight.opt, normwt=TRUE)),
           "Xdata2.mean"=wtd.mean(the_data.1.1$Xdata2, the_data.1.1$weight.opt), 
           "Xdata2.SD"  =sqrt(wtd.var(the_data.1.1$Xdata2, the_data.1.1$weight.opt, normwt=TRUE)),
           "Xdata3.mean"=wtd.mean(the_data.1.1$Xdata3, the_data.1.1$weight.opt), 
           "Xdata3.SD"  =sqrt(wtd.var(the_data.1.1$Xdata3, the_data.1.1$weight.opt, normwt=TRUE)),
           "Xdata4.mean"=wtd.mean(the_data.1.1$Xdata4, the_data.1.1$weight.opt), 
           "Xdata4.SD"  =NA,
           "Xdata5.mean"=wtd.mean(the_data.1.1$Xdata5, the_data.1.1$weight.opt), 
           "Xdata5.SD"  =sqrt(wtd.var(the_data.1.1$Xdata5, the_data.1.1$weight.opt, normwt=TRUE)))
cat("\n Target Summary Stats \n")
data.frame("Xdata1.mean"=means.1.1[1], "Xdata1.SD"=SDs.1.1[1],
           "Xdata2.mean"=means.1.1[2], "Xdata2.SD"=SDs.1.1[2],
           "Xdata3.mean"=means.1.1[3], "Xdata3.SD"=SDs.1.1[3],
           "Xdata4.mean"=means.1.1[4], "Xdata4.SD"=SDs.1.1[4],
           "Xdata5.mean"=means.1.1[5], "Xdata5.SD"=SDs.1.1[5])


# control treatment
cat("\n --- Control ---")
cat("\n Unweighted Summary Stats \n")
data.frame("Xdata1.mean"=mean(the_data.1.0$Xdata1), "Xdata1.SD"=sd(the_data.1.0$Xdata1),
           "Xdata2.mean"=mean(the_data.1.0$Xdata2), "Xdata2.SD"=sd(the_data.1.0$Xdata2),
           "Xdata3.mean"=mean(the_data.1.0$Xdata3), "Xdata3.SD"=sd(the_data.1.0$Xdata3),
           "Xdata4.mean"=mean(the_data.1.0$Xdata4), "Xdata4.SD"=NA,                       # 07 JUL 2026: correction to assigned column name (Xdata4.mean vs Xdata3.mean. No effect on result)
           "Xdata5.mean"=mean(the_data.1.0$Xdata5), "Xdata5.SD"=sd(the_data.1.0$Xdata5))  # 07 JUL 2026: correction to assigned column name (Xdata5.mean vs Xdata4.mean. No effect on result)
cat("\n Weighted Summary Stats \n")
data.frame("ESS"=ESS(the_data.1.0$weight.opt))
data.frame("Xdata1.mean"=wtd.mean(the_data.1.0$Xdata1, the_data.1.0$weight.opt), 
           "Xdata1.SD"  =sqrt(wtd.var(the_data.1.0$Xdata1, the_data.1.0$weight.opt, normwt=TRUE)),
           "Xdata2.mean"=wtd.mean(the_data.1.0$Xdata2, the_data.1.0$weight.opt), 
           "Xdata2.SD"  =sqrt(wtd.var(the_data.1.0$Xdata2, the_data.1.0$weight.opt, normwt=TRUE)),
           "Xdata3.mean"=wtd.mean(the_data.1.0$Xdata3, the_data.1.0$weight.opt), 
           "Xdata3.SD"  =sqrt(wtd.var(the_data.1.0$Xdata3, the_data.1.0$weight.opt, normwt=TRUE)),
           "Xdata4.mean"=wtd.mean(the_data.1.0$Xdata4, the_data.1.0$weight.opt), 
           "Xdata4.SD"  =NA,
           "Xdata5.mean"=wtd.mean(the_data.1.0$Xdata5, the_data.1.0$weight.opt), 
           "Xdata5.SD"  =sqrt(wtd.var(the_data.1.0$Xdata5, the_data.1.0$weight.opt, normwt=TRUE)))
cat("\n Target Summary Stats \n")
data.frame("Xdata1.mean"=means.1.0[1], "Xdata1.SD"=SDs.1.0[1],
           "Xdata2.mean"=means.1.0[2], "Xdata2.SD"=SDs.1.0[2],
           "Xdata3.mean"=means.1.0[3], "Xdata3.SD"=SDs.1.0[3],
           "Xdata4.mean"=means.1.0[4], "Xdata4.SD"=SDs.1.0[4],
           "Xdata5.mean"=means.1.0[5], "Xdata5.SD"=SDs.1.0[5])


# clean the workspace
rm(data.1, IPD.1.1, IPD.1.0, MAIC.1.1, MAIC.1.0, means.1.1, means.1.0, SDs.1.1, SDs.1.0)


#####################################################################################


### Scenario 15 (separate weights per treatment). 4 continuous covariates (Mean [SD]) and one binary (Mean=prop) [xdata2] ###

## Active treatment ##
IPD.15.1      <- data.15[data.15$TREAT==1,1:8]                                    # IPD
means.15.1    <- as.numeric(data.15[1,c("TMT1", "TMT2", "TMT3", "TMT4", "TMT5")]) # Target means
SDs.15.1      <- as.numeric(data.15[1,c("TST1", "TST2", "TST3", "TST4", "TST5")]) # Target SDs
SDs.15.1[2]   <- NA                                                               # Xdata2 is binary, so set TST2 to NA (i.e. no match on SD)
MAIC.15.1     <- MAIC(the_data   = IPD.15.1,                                      # estimate weights. NB: estimates both SIGN (weight.maic) and JACK (weight.opt)
                      id         = "SUBJID", 
                      match_mean = c("Xdata1", "Xdata2", "Xdata3", "Xdata4", "Xdata5"), 
                      agg_means  = means.15.1, 
                      agg_sds    = SDs.15.1, 
                      accuracy   = 0.000000000001)
the_data.15.1 <- MAIC.15.1$the_data                                               # extract the data frame from the list object
# Warning message:
# In lsei(A = A, B = B, E = E, F = F, G = G, H = H) :
# LSEI error: inequalities contradictory

## Control treatment ##
IPD.15.0      <- data.15[data.15$TREAT==0,1:8]                                   # IPD
means.15.0    <- as.numeric(data.15[1,c("TMP1", "TMP2", "TMP3", "TMP4", "TMP5")]) # Target means
SDs.15.0      <- as.numeric(data.15[1,c("TSP1", "TSP2", "TSP3", "TSP4", "TSP5")]) # Target SDs
SDs.15.0[2]   <- NA                                                              # Xdata2 is binary, so set TST2 to NA (i.e. no match on SD)
MAIC.15.0     <- MAIC(the_data   = IPD.15.0,                                     # estimate weights. NB: estimates both SIGN (weight.maic) and JACK (weight.opt)
                     id         = "SUBJID", 
                     match_mean = c("Xdata1", "Xdata2", "Xdata3", "Xdata4", "Xdata5"), 
                     agg_means  = means.15.0, 
                     agg_sds    = SDs.15.0, 
                     accuracy   = 0.000000000001)
the_data.15.0 <- MAIC.15.0$the_data                                              # extract the data frame from the list object
# Warning message:
# In lsei(A = A, B = B, E = E, F = F, G = G, H = H) :
# LSEI error: inequalities contradictory

# Pool the treatments and then try estimating the weights
IPD.15      <- data.15[,1:8]                                                                    # IPD
means.15    <- print(((dim(IPD.15)[1]/2)*means.15.1 + (dim(IPD.15)[1]/2)*means.15.0)/(dim(IPD.15)[1])) # Target means (pooled)
means.15.1 ; means.15.0
SDs.15      <- print(sqrt((((dim(IPD.15)[1]/2)-1)*SDs.15.1^2 + ((dim(IPD.15)[1]/2)-1)*SDs.15.0^2)/(dim(IPD.15)[1] - 2))) # Target SDs (pooled)
SDs.15.1 ; SDs.15.0
MAIC.15     <- MAIC(the_data   = IPD.15,                                                        # estimate weights. NB: estimates both SIGN (weight.maic) and JACK (weight.opt)
                      id         = "SUBJID", 
                      match_mean = c("Xdata1", "Xdata2", "Xdata3", "Xdata4", "Xdata5"), 
                      agg_means  = means.15, 
                      agg_sds    = SDs.15, 
                      accuracy   = 0.000000000001)
the_data.15 <- MAIC.15$the_data                                                                 # extract the data frame from the list object

## Check the matching (Xdata1/3/4/5, mean and standard deviation + Xdata2, mean) ##
cat("\n ----- SIMULATION 15 -----")

# Pooled treatments
cat("\n --- Pooled ---")
cat("\n Unweighted Summary Stats \n")
data.frame("Xdata1.mean"=mean(the_data.15$Xdata1), "Xdata1.SD"=sd(the_data.15$Xdata1),
           "Xdata2.mean"=mean(the_data.15$Xdata2), "Xdata2.SD"=NA,
           "Xdata3.mean"=mean(the_data.15$Xdata3), "Xdata3.SD"=sd(the_data.15$Xdata3),
           "Xdata4.mean"=mean(the_data.15$Xdata4), "Xdata4.SD"=sd(the_data.15$Xdata4),
           "Xdata5.mean"=mean(the_data.15$Xdata5), "Xdata5.SD"=sd(the_data.15$Xdata5))
cat("\n Weighted Summary Stats \n")
data.frame("ESS"=ESS(the_data.15$weight.opt))
data.frame("Xdata1.mean"=wtd.mean(the_data.15$Xdata1, the_data.15$weight.opt), 
           "Xdata1.SD"  =sqrt(wtd.var(the_data.15$Xdata1, the_data.15$weight.opt, normwt=TRUE)),
           "Xdata2.mean"=wtd.mean(the_data.15$Xdata2, the_data.15$weight.opt), 
           "Xdata2.SD"  =NA,
           "Xdata3.mean"=wtd.mean(the_data.15$Xdata3, the_data.15$weight.opt), 
           "Xdata3.SD"  =sqrt(wtd.var(the_data.15$Xdata3, the_data.15$weight.opt, normwt=TRUE)),
           "Xdata4.mean"=wtd.mean(the_data.15$Xdata4, the_data.15$weight.opt), 
           "Xdata4.SD"  =sqrt(wtd.var(the_data.15$Xdata4, the_data.15$weight.opt, normwt=TRUE)),
           "Xdata5.mean"=wtd.mean(the_data.15$Xdata5, the_data.15$weight.opt), 
           "Xdata5.SD"  =sqrt(wtd.var(the_data.15$Xdata5, the_data.15$weight.opt, normwt=TRUE)))
cat("\n Target Summary Stats \n")
data.frame("Xdata1.mean"=means.15[1], "Xdata1.SD"=SDs.15[1],
           "Xdata2.mean"=means.15[2], "Xdata2.SD"=SDs.15[2],
           "Xdata3.mean"=means.15[3], "Xdata3.SD"=SDs.15[3],
           "Xdata4.mean"=means.15[4], "Xdata4.SD"=SDs.15[4],
           "Xdata5.mean"=means.15[5], "Xdata5.SD"=SDs.15[5])
cat("\n Errors \n")
sprintf("%.20f", data.frame("Xdata1.mean"=wtd.mean(the_data.15$Xdata1, the_data.15$weight.opt) - means.15[1], 
                            "Xdata1.SD"  =sqrt(wtd.var(the_data.15$Xdata1, the_data.15$weight.opt, normwt=TRUE)) - SDs.15[1],
                            "Xdata2.mean"=wtd.mean(the_data.15$Xdata2, the_data.15$weight.opt) - means.15[2], 
                            "Xdata2.SD"  =NA,
                            "Xdata3.mean"=wtd.mean(the_data.15$Xdata3, the_data.15$weight.opt) - means.15[3], 
                            "Xdata3.SD"  =sqrt(wtd.var(the_data.15$Xdata3, the_data.15$weight.opt, normwt=TRUE)) - SDs.15[3],
                            "Xdata4.mean"=wtd.mean(the_data.15$Xdata4, the_data.15$weight.opt) - means.15[4], 
                            "Xdata4.SD"  =sqrt(wtd.var(the_data.15$Xdata4, the_data.15$weight.opt, normwt=TRUE)) - SDs.15[4],
                            "Xdata5.mean"=wtd.mean(the_data.15$Xdata5, the_data.15$weight.opt) - means.15[5], 
                            "Xdata5.SD"  =sqrt(wtd.var(the_data.15$Xdata5, the_data.15$weight.opt, normwt=TRUE)) - SDs.15[5]))


#hist(the_data.15$weight.opt)
#dev.off()

# clean the workspace
rm(data.15, IPD.15.1, IPD.15.0, MAIC.15.1, MAIC.15.0, means.15.1, means.15.0, SDs.15.1, SDs.15.0,
   IPD.15, means.15, SDs.15)


#####################################################################################


### Scenario 247 (separate weights per treatment). 2 continuous covariates (Mean [SD]) and two binary (Mean=prop) [xdata1/2] ###

## Active treatment ##
IPD.247.1     <- data.247[data.247$TREAT==1,1:8]                                   # IPD
means.247.1   <- as.numeric(data.247[1,c("TMT1", "TMT2", "TMT3", "TMT5")])         # Target means
SDs.247.1     <- as.numeric(data.247[1,c("TST1", "TST2", "TST3", "TST5")])         # Target SDs
SDs.247.1[1]  <- SDs.247.1[2] <- NA                                                # Xdata1/2 are binary, so set TST1/2 to NA (i.e. no match on SD)
MAIC.247.1    <- MAIC(the_data   = IPD.247.1,                                      # estimate weights. NB: estimates both SIGN (weight.maic) and JACK (weight.opt)
                      id         = "SUBJID", 
                      match_mean = c("Xdata1", "Xdata2", "Xdata3", "Xdata5"), 
                      agg_means  = means.247.1, 
                      agg_sds    = SDs.247.1, 
                      accuracy   = 0.000000000001)
the_data.247.1 <- MAIC.247.1$the_data                                              # extract the data frame from the list object

## Control treatment ##
IPD.247.0      <- data.247[data.247$TREAT==0,1:8]                                   # IPD
means.247.0    <- as.numeric(data.247[1,c("TMP1", "TMP2", "TMP3", "TMP5")])         # Target means
SDs.247.0      <- as.numeric(data.247[1,c("TSP1", "TSP2", "TSP3", "TSP5")])         # Target SDs
SDs.247.0[1]   <- SDs.247.0[2] <- NA                                                # Xdata1/2 are binary, so set TST1/2 to NA (i.e. no match on SD)
MAIC.247.0     <- MAIC(the_data   = IPD.247.0,                                      # estimate weights. NB: estimates both SIGN (weight.maic) and JACK (weight.opt)
                      id         = "SUBJID", 
                      match_mean = c("Xdata1", "Xdata2", "Xdata3", "Xdata5"), 
                      agg_means  = means.247.0, 
                      agg_sds    = SDs.247.0, 
                      accuracy   = 0.000000000001)
the_data.247.0 <- MAIC.247.0$the_data                                               # extract the data frame from the list object
# Warning message:
# In lsei(A = A, B = B, E = E, F = F, G = G, H = H) :
# LSEI error: inequalities contradictory

# Pool the treatments and then try estimating the weights
IPD.247      <- data.247[,1:8]                                                                    # IPD
means.247    <- print(((dim(IPD.247)[1]/2)*means.247.1 + (dim(IPD.247)[1]/2)*means.247.0)/(dim(IPD.247)[1])) # Target means (pooled)
means.247.1 ; means.247.0
SDs.247      <- print(sqrt((((dim(IPD.247)[1]/2)-1)*SDs.247.1^2 + ((dim(IPD.247)[1]/2)-1)*SDs.247.0^2)/(dim(IPD.247)[1] - 2))) # Target SDs (pooled)
SDs.247.1 ; SDs.247.0
MAIC.247     <- MAIC(the_data   = IPD.247,                                                        # estimate weights. NB: estimates both SIGN (weight.maic) and JACK (weight.opt)
                    id         = "SUBJID", 
                    match_mean = c("Xdata1", "Xdata2", "Xdata3", "Xdata5"), 
                    agg_means  = means.247, 
                    agg_sds    = SDs.247, 
                    accuracy   = 0.000000000001)
the_data.247 <- MAIC.247$the_data                                                                 # extract the data frame from the list object

## Check the matching (Xdata3/5, mean and standard deviation + Xdata1/2, mean) ##
cat("\n ----- SIMULATION 247 -----")

# Pooled treatments
cat("\n --- Pooled ---")
cat("\n Unweighted Summary Stats \n")
data.frame("Xdata1.mean"=mean(the_data.247$Xdata1), "Xdata1.SD"=NA,
           "Xdata2.mean"=mean(the_data.247$Xdata2), "Xdata2.SD"=NA,
           "Xdata3.mean"=mean(the_data.247$Xdata3), "Xdata3.SD"=sd(the_data.247$Xdata3),
           "Xdata5.mean"=mean(the_data.247$Xdata5), "Xdata5.SD"=sd(the_data.247$Xdata5))
cat("\n Weighted Summary Stats \n")
data.frame("ESS"=ESS(the_data.247$weight.opt))
data.frame("Xdata1.mean"=wtd.mean(the_data.247$Xdata1, the_data.247$weight.opt), 
           "Xdata1.SD"  =NA,
           "Xdata2.mean"=wtd.mean(the_data.247$Xdata2, the_data.247$weight.opt), 
           "Xdata2.SD"  =NA,
           "Xdata3.mean"=wtd.mean(the_data.247$Xdata3, the_data.247$weight.opt), 
           "Xdata3.SD"  =sqrt(wtd.var(the_data.247$Xdata3, the_data.247$weight.opt, normwt=TRUE)),
           "Xdata5.mean"=wtd.mean(the_data.247$Xdata5, the_data.247$weight.opt), 
           "Xdata5.SD"  =sqrt(wtd.var(the_data.247$Xdata5, the_data.247$weight.opt, normwt=TRUE)))
cat("\n Target Summary Stats \n")
data.frame("Xdata1.mean"=means.247[1], "Xdata1.SD"=SDs.247[1],
           "Xdata2.mean"=means.247[2], "Xdata2.SD"=SDs.247[2],
           "Xdata3.mean"=means.247[3], "Xdata3.SD"=SDs.247[3],
           "Xdata5.mean"=means.247[4], "Xdata5.SD"=SDs.247[4])
cat("\n Errors \n")
sprintf("%.20f", data.frame("Xdata1.mean"=wtd.mean(the_data.247$Xdata1, the_data.247$weight.opt) - means.247[1], 
                            "Xdata1.SD"  =NA,
                            "Xdata2.mean"=wtd.mean(the_data.247$Xdata2, the_data.247$weight.opt) - means.247[2], 
                            "Xdata2.SD"  =NA,
                            "Xdata3.mean"=wtd.mean(the_data.247$Xdata3, the_data.247$weight.opt) - means.247[3], 
                            "Xdata3.SD"  =sqrt(wtd.var(the_data.247$Xdata3, the_data.247$weight.opt, normwt=TRUE)) - SDs.247[3],
                            "Xdata5.mean"=wtd.mean(the_data.247$Xdata5, the_data.247$weight.opt) - means.247[4], 
                            "Xdata5.SD"  =sqrt(wtd.var(the_data.247$Xdata5, the_data.247$weight.opt, normwt=TRUE)) - SDs.247[4]))

#hist(the_data.247$weight.opt)
#dev.off()

# clean the workspace
rm(data.247, IPD.247.1, IPD.247.0, MAIC.247.1, MAIC.247.0, means.247.1, means.247.0, SDs.247.1, SDs.247.0,
   IPD.247, means.247, SDs.247)


#####################################################################################


### Scenario 320 (separate weights per treatment). 3 continuous covariates (Mean [SD]) and two binary (Mean=prop) [xdata2/3] ###

## Active treatment ##
IPD.320.1     <- data.320[data.320$TREAT==1,1:8]                                    # IPD
means.320.1   <- as.numeric(data.320[1,c("TMT1", "TMT2", "TMT3", "TMT4", "TMT5")])  # Target means
SDs.320.1     <- as.numeric(data.320[1,c("TST1", "TST2", "TST3", "TST4", "TST5")])  # Target SDs
SDs.320.1[2]  <- SDs.320.1[3] <- NA                                                 # Xdata2/3 are binary, so set TST2/3 to NA (i.e. no match on SD)
MAIC.320.1    <- MAIC(the_data   = IPD.320.1,                                       # estimate weights. NB: estimates both SIGN (weight.maic) and JACK (weight.opt)
                      id         = "SUBJID", 
                      match_mean = c("Xdata1", "Xdata2", "Xdata3", "Xdata4", "Xdata5"), 
                      agg_means  = means.320.1, 
                      agg_sds    = SDs.320.1, 
                      accuracy   = 0.000000000001)
the_data.320.1 <- MAIC.320.1$the_data                                               # extract the data frame from the list object

## Control treatment ##
IPD.320.0     <- data.320[data.320$TREAT==0,1:8]                                    # IPD
means.320.0   <- as.numeric(data.320[1,c("TMP1", "TMP2", "TMP3", "TMP4", "TMP5")])  # Target means
SDs.320.0     <- as.numeric(data.320[1,c("TSP1", "TSP2", "TSP3", "TSP4", "TSP5")])  # Target SDs
SDs.320.0[2]  <- SDs.320.0[3] <- NA                                                 # Xdata2/3 are binary, so set TST2/3 to NA (i.e. no match on SD)
MAIC.320.0    <- MAIC(the_data   = IPD.320.0,                                       # estimate weights. NB: estimates both SIGN (weight.maic) and JACK (weight.opt)
                      id         = "SUBJID", 
                      match_mean = c("Xdata1", "Xdata2", "Xdata3", "Xdata4", "Xdata5"), 
                      agg_means  = means.320.0, 
                      agg_sds    = SDs.320.0, 
                      accuracy   = 0.000000000001)
the_data.320.0 <- MAIC.320.0$the_data                                               # extract the data frame from the list object
# Warning message:
# In lsei(A = A, B = B, E = E, F = F, G = G, H = H) :
# LSEI error: inequalities contradictory

# Pool the treatments and then try estimating the weights
IPD.320      <- data.320[,1:8]                                                                    # IPD
means.320    <- print(((dim(IPD.320)[1]/2)*means.320.1 + (dim(IPD.320)[1]/2)*means.320.0)/(dim(IPD.320)[1])) # Target means (pooled)
means.320.1 ; means.320.0
SDs.320      <- print(sqrt((((dim(IPD.320)[1]/2)-1)*SDs.320.1^2 + ((dim(IPD.320)[1]/2)-1)*SDs.320.0^2)/(dim(IPD.320)[1] - 2))) # Target SDs (pooled)
SDs.320.1 ; SDs.320.0
MAIC.320     <- MAIC(the_data   = IPD.320,                                                        # estimate weights. NB: estimates both SIGN (weight.maic) and JACK (weight.opt)
                    id         = "SUBJID", 
                    match_mean = c("Xdata1", "Xdata2", "Xdata3", "Xdata4", "Xdata5"), 
                    agg_means  = means.320, 
                    agg_sds    = SDs.320, 
                    accuracy   = 0.000000000001)
the_data.320 <- MAIC.320$the_data   

## Check the matching (Xdata1/4/5, mean and standard deviation + Xdata2/3, mean) ##
cat("\n ----- SIMULATION 320 -----")

# Pooled treatments
cat("\n --- Pooled ---")
cat("\n Unweighted Summary Stats \n")
data.frame("Xdata1.mean"=mean(the_data.320$Xdata1), "Xdata1.SD"=sd(the_data.320$Xdata1),
           "Xdata2.mean"=mean(the_data.320$Xdata2), "Xdata2.SD"=NA,
           "Xdata3.mean"=mean(the_data.320$Xdata3), "Xdata3.SD"=NA,
           "Xdata4.mean"=mean(the_data.320$Xdata4), "Xdata4.SD"=sd(the_data.320$Xdata4),
           "Xdata5.mean"=mean(the_data.320$Xdata5), "Xdata5.SD"=sd(the_data.320$Xdata5))
cat("\n Weighted Summary Stats \n")
data.frame("ESS"=ESS(the_data.320$weight.opt))
data.frame("Xdata1.mean"=wtd.mean(the_data.320$Xdata1, the_data.320$weight.opt), 
           "Xdata1.SD"  =sqrt(wtd.var(the_data.320$Xdata1, the_data.320$weight.opt, normwt=TRUE)),
           "Xdata2.mean"=wtd.mean(the_data.320$Xdata2, the_data.320$weight.opt), 
           "Xdata2.SD"  =NA,
           "Xdata3.mean"=wtd.mean(the_data.320$Xdata3, the_data.320$weight.opt), 
           "Xdata3.SD"  =NA,
           "Xdata4.mean"=wtd.mean(the_data.320$Xdata4, the_data.320$weight.opt), 
           "Xdata4.SD"  =sqrt(wtd.var(the_data.320$Xdata4, the_data.320$weight.opt, normwt=TRUE)),           
           "Xdata5.mean"=wtd.mean(the_data.320$Xdata5, the_data.320$weight.opt), 
           "Xdata5.SD"  =sqrt(wtd.var(the_data.320$Xdata5, the_data.320$weight.opt, normwt=TRUE)))
cat("\n Target Summary Stats \n")
data.frame("Xdata1.mean"=means.320[1], "Xdata1.SD"=SDs.320[1],
           "Xdata2.mean"=means.320[2], "Xdata2.SD"=SDs.320[2],
           "Xdata3.mean"=means.320[3], "Xdata3.SD"=SDs.320[3],
           "Xdata4.mean"=means.320[4], "Xdata4.SD"=SDs.320[4],           
           "Xdata5.mean"=means.320[5], "Xdata5.SD"=SDs.320[5])
cat("\n Errors \n")
sprintf("%.20f", data.frame("Xdata1.mean"=wtd.mean(the_data.320$Xdata1, the_data.320$weight.opt) - means.320[1], 
                            "Xdata1.SD"  =sqrt(wtd.var(the_data.320$Xdata1, the_data.320$weight.opt, normwt=TRUE)) - SDs.320[1],
                            "Xdata2.mean"=wtd.mean(the_data.320$Xdata2, the_data.320$weight.opt) - means.320[2], 
                            "Xdata2.SD"  =NA,
                            "Xdata3.mean"=wtd.mean(the_data.320$Xdata3, the_data.320$weight.opt) - means.320[3], 
                            "Xdata3.SD"  =NA,
                            "Xdata4.mean"=wtd.mean(the_data.320$Xdata4, the_data.320$weight.opt) - means.320[4], 
                            "Xdata4.SD"  =sqrt(wtd.var(the_data.320$Xdata4, the_data.320$weight.opt, normwt=TRUE)) - SDs.320[4],
                            "Xdata5.mean"=wtd.mean(the_data.320$Xdata5, the_data.320$weight.opt) - means.320[5], 
                            "Xdata5.SD"  =sqrt(wtd.var(the_data.320$Xdata5, the_data.320$weight.opt, normwt=TRUE)) - SDs.320[5]))


#hist(the_data.320$weight.opt)
#dev.off()

# clean the workspace
rm(data.320, IPD.320.1, IPD.320.0, MAIC.320.1, MAIC.320.0, means.320.1, means.320.0, SDs.320.1, SDs.320.0,
   IPD.320, means.320, SDs.320)




#####################################################################################
 

# EXTRACT WEIGHTS ###################################################################

# stack weights across scenarios, align variable names with production, keep only relevant variables,
# and then sort by SIM and SUBJID
WEIGHTSALL <- dplyr:: bind_rows(the_data.1439.1, the_data.1439.0, # Scenario 1439 (Active, Control)
                                the_data.1.1   , the_data.1.0   , # Scenario 1 (Active, Control)
                                the_data.15,                      # Scenario 15 (Pooled)
                                the_data.247,                     # Scenario 247 (Pooled)
                                the_data.320) %>%                 # Scenario 320 (Pooled)
                      mutate(SIM=Scenario) %>%
                      select(SIM, SUBJID, TREAT, weight.opt) %>%
                      mutate(POOLFL=ifelse(SIM %in% c(1,1439),"","Y")) %>%
                      arrange(SIM, SUBJID)
#table(WEIGHTSALL$SIM, WEIGHTSALL$POOLFL)

# calculate the sum of the weights and number of subjects (for use in re-scaling them to have a mean of 1)
SUM.WEIGHTS1 <- WEIGHTSALL %>% 
                filter(POOLFL!="Y") %>%
                group_by(SIM, TREAT) %>%
                summarise(SUMWGTS=sum(weight.opt), N=n()) # treatment-specific weighting scenarios
SUM.WEIGHTS2 <- WEIGHTSALL %>% 
                filter(POOLFL=="Y") %>%
                group_by(SIM) %>%
                summarise(SUMWGTS=sum(weight.opt), N=n()) # pooled weighting scenarios


# merge them back on to the unscaled IPWs and calculate the re-scaled IPWs
V_FINAL_WEIGHTSALL1 <- full_join(x = WEIGHTSALL[WEIGHTSALL$POOLFL!="Y",],
                                 y = SUM.WEIGHTS1,
                                 by = c("SIM" = "SIM", "TREAT" = "TREAT")
                                 ) %>%
                       mutate(JACK_WEIGHTS=N*(weight.opt/SUMWGTS)) %>%
                       select(SIM, SUBJID, TREAT, JACK_WEIGHTS, POOLFL) # treatment-specific weighting scenarios

V_FINAL_WEIGHTSALL2 <- full_join(x = WEIGHTSALL[WEIGHTSALL$POOLFL=="Y",],
                                 y = SUM.WEIGHTS2,
                                 by = c("SIM" = "SIM")
                                 ) %>%
                                 mutate(JACK_WEIGHTS=N*(weight.opt/SUMWGTS)) %>%
                                 select(SIM, SUBJID, TREAT, JACK_WEIGHTS, POOLFL) # treatment-specific weighting scenarios


# set the scenarios back together and sort
V_FINAL_WEIGHTSALL <- rbind(V_FINAL_WEIGHTSALL1, V_FINAL_WEIGHTSALL2) %>%
                      arrange(SIM, SUBJID)

# COLLATE DIAGNOSTIC MEASURES #######################################################

V_FINAL_DIAGNOSTICSALL <- V_FINAL_WEIGHTSALL %>% 
                          group_by(SIM, TREAT) %>%              # perform summaries by simulation and treatment
                          summarise(N       =n(),               # Number of patients
                                    ESS     =ESS(JACK_WEIGHTS), # ESS
                                    ESSP    =ESS/N,             # ESS as proportion of number of patients
                                    MAXW    =max(JACK_WEIGHTS), # maximum weight
                                    VARIANCE=var(JACK_WEIGHTS)  # variance of weights
                                    ) %>%
                          mutate(POOLFL=ifelse(SIM %in% c(1,1439),"","Y")) %>%
                          arrange(SIM, TREAT)


## END OF PROGRAM ##
