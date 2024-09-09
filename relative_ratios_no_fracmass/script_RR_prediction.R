#script for the prediction of the relative ratios of lipids

#libraries
library(flexmix)
 
#data preparation
database <-  read.csv(file="./Lipid_Maps_and_SWISS_lipids.csv", header=TRUE)
reduced_database$Halogens <- database$Br + database$Cl
reduced_database$decimal <- reduced_database$mass1 - floor(reduced_database$mass1)
reduced_database <- reduced_database[order(reduced_database$mass1),]
reduced_database$m0 <- reduced_database$mass1/1000

reduced_noS <- subset(reduced_database, S==0)
reduced_noS_noX <- subset(reduced_noS, Halogens==0)
reduced_S <- subset(reduced_database, S!=0)
reduced_S_noX <- subset(reduced_S, Halogens==0)


#loading fitted models
fModel21.noS <- readRDS("resultRR21.noS.rds")
fModel32.noS <- readRDS("resultRR32.noS.rds")
fModel43.noS <- readRDS("resultRR43.noS.rds")
fModel54.noS <- readRDS("resultRR54.noS.rds")

fModel21.S <- readRDS("resultRR21.S.rds")
fModel32.S <- readRDS("resultRR32.S.rds")
fModel43.S <- readRDS("resultRR43.S.rds")
fModel54.S <- readRDS("resultRR54.S.rds")

#predicting
toPred <- data.frame("m0"=reduced_noS$m0[1100:1105], "m21"=reduced_noS$m21[1100:1105], "m32"=reduced_noS$m32[1100:1105])

resRR21.noS <- predRR21.agg(x=toPred, fModel21.noS, conf.int = TRUE, conf.level=0.95, S=FALSE)
resRR21.S <- predRR21.agg(x=toPred, fModel21.S, conf.int = TRUE, conf.level=0.95, S=TRUE)

#Halogens are excluded in RR32 predictor training
resRR32.noS <- predRR32.agg(x=toPred, fModel32.noS, conf.int = TRUE, conf.level=0.95, S=FALSE)
resRR32.S <- predRR32.agg(x=toPred, fModel32.S, conf.int = TRUE, conf.level=0.95, S=TRUE)

resRR43.noS <- predRR43.agg(x=toPred, fModel43.noS, conf.int = TRUE, conf.level=0.95, S=FALSE)
resRR43.S <- predRR43.agg(x=toPred, fModel43.S, conf.int = TRUE, conf.level=0.95, S=TRUE)

#Halogens are excluded in RR54 predictor training
resRR54.noS <- predRR54.agg(x=toPred, fModel54.noS, conf.int = TRUE, conf.level=0.95, S=FALSE)
resRR54.S <- predRR54.agg(x=toPred, fModel54.S, conf.int = TRUE, conf.level=0.95, S=TRUE)

#check results
resRR21.noS$truth <- reduced_noS$RR21[1100:1105]
which(resRR21.noS$predLower <= resRR21.noS$truth & resRR21.noS$predUpper >= resRR21.noS$truth)

resRR21.S$truth <- reduced_noS$RR21[1100:1105]
which(resRR21.S$predLower <= resRR21.S$truth & resRR21.S$predUpper >= resRR21.S$truth)

#check results wrt to element composition (C, H, N, O, S, P, X)
