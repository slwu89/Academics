######################################
######PH252D FINAL PROJET SCRIPT######
######################################


######################
###DATA & LIBRARIES###
######################

#analysis
library(parallel)
library(SuperLearner)
library(doSNOW)
library(snow)
library(doParallel)

#visualization
library(ggplot2)
library(gridExtra)
library(reshape2)
library(plyr)
library(tmle)

#data
setwd("C:/Users/WuS/Dropbox/Academics/Spring 2016/PH252D/final_project/")
setwd("C:/Users/ASRock Z87 Pro4/Dropbox/Academics/Spring 2016/PH252D/final_project/")
data <- readRDS(file="data.no.outliers.RDS")
data <- data[,c(5,4,1,2,3,6,7,8,9,10,11,12,13)]