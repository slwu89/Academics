######################################
######PH252D FINAL PROJET SCRIPT######
######################################


######################
###DATA & LIBRARIES###
######################

#analysis (please load these packages to run 2 part bootstrap procedure)
library(parallel)
library(SuperLearner)
library(doSNOW)
library(snow)
library(doParallel)
library(tmle)

#visualization (please wait until after running part 1 and 2 of bootstrap process to load these packages,
#exporting them to the parallel cluster will slow down bootstrapping loop considerably)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(plyr)

#data (set working directory to where 'data.no.outliers.RDS' is)
setwd("C:/Users/WuS/Dropbox/Academics/Spring 2016/PH252D/final_project/")
setwd("C:/Users/ASRock Z87 Pro4/Dropbox/Academics/Spring 2016/PH252D/final_project/")

#load data and resort columns for TMLE package
data <- readRDS(file="data.no.outliers.RDS")
data_index <- c("Y","A","sex","age","v02","age_sq","age_cos","age_sin","age_log","v02_sq","v02_cos","v02_sin","v02_log")
data <- data[,match(data_index,names(data))]

#SuperLearner library
slLib <- c("SL.polymars","SL.glmnet","SL.gam","SL.glm","SL.rpartPrune")


#########################################
###2 PART BOOTSTRAP PROCEDURE FOR TMLE###
#########################################

###Part 1: Estimate Q0 and g0###
boot_partA <- function(data,slLib,verbose=TRUE){
  
  ###step 1 of TMLE algorthm: estimate E0(Y|A,W) by Qn(A,W)###
  #prepare data
  n <- nrow(data)
  data$Y <- as.numeric(data$Y) #SL does not like integers
  data$A <- as.numeric(data$A) #SL does not like integers
  dat_a1 <- data ; dat_a1$A <- 1
  dat_a0 <- data ; dat_a0$A <- 0
  dat_countFact <- rbind(data,dat_a1,dat_a0) ; dat_countFact <- dat_countFact[,-1]
  
  #create initial estimate Qn_0 (untargeted estimator)
  Qn_0 <- SuperLearner(Y=data$Y,X=data[,!names(data) %in% c("Y")],newX=dat_countFact[,!names(dat_countFact) %in% c("Y")],SL.library=slLib,family="gaussian",verbose=verbose,cvControl=list(V=10,shuffle=TRUE))
  
  #extract predictions from Qn_0 (untargeted estimator)
  Qn_0_predA <- Qn_0$SL.predict[1:n] #predicted Y given observed covariates
  Qn_0_pred1 <- Qn_0$SL.predict[(n+1):(2*n)] #counterfactual Y1
  Qn_0_pred0 <- Qn_0$SL.predict[(2*n+1):(3*n)] #counterfactual Y0
  
  #evaluate simple substitution estimator
  sub_est <- mean(Qn_0_pred1 - Qn_0_pred0)
  
  ####step 2 of TMLE algorithm: estimate P0(A|W) by Gn(A|W)###
  #estimation
  gn_hat <- SuperLearner(Y=data$A,X=data[,!names(data) %in% c("Y","A")],SL.library=slLib[grep("SL.glm.[1-9]",slLib,invert=TRUE)],family="binomial",verbose=verbose,cvControl=list(V=10,shuffle=TRUE))
  gn_hat1 <- gn_hat$SL.predict #propensity scores
  gn_hat0 <- (1 - gn_hat1) #propensity scores
  
  #get predicted probability of A given W for each subject
  gn_hatAW <- rep(NA,n)
  gn_hatAW[data$A==1] <- gn_hat1[data$A==1]
  gn_hatAW[data$A==0] <- gn_hat0[data$A==0]
  
  #evaluate IPTW and stabilized IPTW estimator
  iptw_est <- mean(as.numeric(data$A==1) * data$Y / gn_hatAW) - mean(as.numeric(data$A==0) * data$Y / gn_hatAW)
  iptw_estStable <- mean(as.numeric(data$A==1) * data$Y / gn_hatAW)/mean(as.numeric(data$A==1) / gn_hatAW) - mean(as.numeric(data$A==0) * data$Y / gn_hatAW)/mean(as.numeric(data$A==1) / gn_hatAW)
  
  #return all output as list
  results <- list(sub_est=sub_est,iptw_est=iptw_est,iptw_estStable=iptw_estStable,gAW=gn_hatAW,gn_hat1=gn_hat1,gn_hat0=gn_hat0,Qn_0_pred1=Qn_0_pred1,Qn_0_pred0=Qn_0_pred0,Qn_0_predA=Qn_0_predA)
  return(results)
}

#Bootstrap B samples for Q0 and g0 estimates
pkg_export <- c(sessionInfo()$basePkgs,names(sessionInfo()$otherPkgs))
B <- 5e3

cl <- makeCluster(spec=detectCores())
registerDoParallel(cl)

system.time(bootA_out <- foreach(i=1:B,.packages=pkg_export,.export=c("boot_partA","slLib","data"),.verbose=TRUE) %dopar% {
  data_b <- data[sample(nrow(data),replace=TRUE),]
  result <- boot_partA(data=data_b,slLib=slLib)
  result$boot_b <- data_b #return the bootstrap data set
  result
})

stopCluster(cl)
rm(cl)

saveRDS(object=bootA_out,file="bootA_out.rds")


####Part 2: Estimate TMLE###
bootA_out <- readRDS(file="bootA_out.rds")
boot_partB_index <- 1:length(bootA_out)

bootB_out <- foreach(i=boot_partB_index,.verbose=TRUE) %do% {
  if(sum(unlist(sapply(bootA_out[[i]],is.na)))>0){
    return(NULL)
  }
  tmle(Y=bootA_out[[i]]$boot_b$Y,A=bootA_out[[i]]$boot_b$A,W=bootA_out[[i]]$boot_b[,3:13],Q=cbind(bootA_out[[i]]$Qn_0_pred0,bootA_out[[i]]$Qn_0_pred1),g1W=bootA_out[[i]]$gn_hat1,family="gaussian")
}

#remove bad samples
bootB_out <- bootB_out[!sapply(bootB_out,is.null)]

saveRDS(object=bootB_out,file="bootB_out.rds")


###################################
###VISUALIZE AND EXTRACT RESULTS###
###################################

bootB_out <- readRDS(file="bootB_out.rds")

#extract data
boot_iptw <- sapply(bootA_out,function(x) {x$iptw_est}) #IPTW estimator
boot_iptwS <- sapply(bootA_out,function(x) {x$iptw_estStable}) #Stabilized IPTW estimator
boot_sub <- sapply(bootA_out,function(x) {x$sub_est}) #Simple substitution estimator
boot_gAW <- lapply(bootA_out,function(x) {x$gAW}) #Estimation of treatment mechanism (gAW)
boot_wt <- lapply(bootA_out,function(x) {1/(x$gAW)}) #Weights
boot_tmle <- sapply(bootB_out,function(x) {x$estimates$ATE$psi}) #TMLE estimator
boot_tmleVar <- sapply(bootB_out,function(x) {x$estimates$ATE$var.psi}) #TMLE variance
boot_tmleCI <- sapply(bootB_out,function(x) {x$estimates$ATE$CI}) #TMLE influence curve CI
boot_tmleCI <- as.data.frame(t(boot_tmleCI))
boot_epsilon <- sapply(bootB_out,function(x) {x$epsilon}) #epsilon
boot_epsilon <- as.data.frame(t(boot_epsilon))
boot_epsilon <- melt(boot_epsilon)
boot_epsilon$sample <- c(1:(nrow(boot_epsilon)/2),1:(nrow(boot_epsilon)/2))
boot_epsilon$value <- boot_epsilon$value + 1e-6 #push away from 0 for log transform

#calculate clever covariate from bootA_out
boot_g1w <- lapply(bootA_out,function(x) {x$gn_hat1})
boot_g0w <- lapply(bootA_out,function(x) {x$gn_hat0})
boot_hAW <- foreach(i=1:length(boot_partB_index),.verbose=TRUE) %do% {
  ans <- as.numeric(data$A)/as.vector(boot_g1w[[i]]) - as.numeric(data$A)/as.vector(boot_g0w[[i]])
  return(ans)
}
boot_hAW <- lapply(boot_hAW,function(x) {
  bad_ind <- which(is.na(x) | is.nan(x) | is.infinite(x))
  return(x[-bad_ind])
})

#IPTW estimator
plot_iptw <- ggplot() +
  geom_histogram(data=as.data.frame(boot_iptw),aes(boot_iptw),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(boot_iptw,na.rm=TRUE),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(boot_iptw,na.rm=TRUE),colour="tomato",lty=3,size=1.05) +
  labs(x="IPTW Estimator (5000 Bootstrap Samples)") +
  theme_bw() +
  theme(axis.title=element_text(size=13),axis.title.y=element_blank())

#Stabilized IPTW estimator
plot_iptwS <- ggplot() +
  geom_histogram(data=as.data.frame(boot_iptwS),aes(boot_iptwS),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(boot_iptwS,na.rm=TRUE),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(boot_iptwS,na.rm=TRUE),colour="tomato",lty=3,size=1.05) +
  labs(x="Stabilized IPTW Estimator (5000 Bootstrap Samples)") +
  theme_bw() +
  theme(axis.title=element_text(size=13),axis.title.y=element_blank())

#Simple substitution estimator
plot_sub <- ggplot() +
  geom_histogram(data=as.data.frame(boot_sub),aes(boot_sub),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(boot_sub,na.rm=TRUE),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(boot_sub,na.rm=TRUE),colour="tomato",lty=3,size=1.05) +
  labs(x="Substitution Estimator (5000 Bootstrap Samples)") +
  theme_bw() +
  theme(axis.title=element_text(size=13),axis.title.y=element_blank())

#Treatment mechanism estimator (gAW)
plot_gAW <- ggplot() +
  geom_histogram(data=melt(boot_gAW),aes(value,colour=as.factor(L1),group=L1,fill=as.factor(L1)),position="identity",alpha=0.5) +
  labs(x="Predicted Probabilities (5000 Bootstrap Samples)") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw() +
  theme(axis.title=element_text(size=13),axis.title.y=element_blank())

#IPTW weights
plot_wt <- ggplot() +
  geom_histogram(data=melt(boot_wt),aes(value,colour=as.factor(L1),group=L1,fill=as.factor(L1)),position="identity",alpha=0.5) +
  labs(x="Weights (5000 Bootstrap Samples)") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw() +
  theme(axis.title=element_text(size=13),axis.title.y=element_blank())

#clever covariate
plot_hAW <- ggplot() +
  geom_histogram(data=melt(boot_hAW),aes(value,colour=as.factor(L1),group=L1,fill=as.factor(L1)),position="identity",alpha=0.5) +
  labs(x="Clever Covariate (5000 Bootstrap Samples)") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw() +
  theme(axis.title=element_text(size=13),axis.title.y=element_blank())

grid.arrange(plot_iptw,plot_iptwS,plot_sub,ncol=3)
grid.arrange(plot_gAW,plot_wt,ncol=2)

#TMLE estimator
plot_tmle <- ggplot() +
  geom_histogram(data=data.frame(boot_tmle),aes(boot_tmle),fill="steelblue",colour="black",bins=30) +
  geom_vline(xintercept=mean(boot_tmle,na.rm=TRUE),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(boot_tmle,na.rm=TRUE),colour="tomato",lty=3,size=1.05) +
  labs(x="TMLE Estimator (5000 Bootstrap Samples)") +
  theme_bw() +
  theme(axis.title=element_text(size=13),axis.title.y=element_blank())

#TMLE variance
plot_tmleVar <- ggplot() +
  geom_boxplot(data=data.frame(boot_tmleVar),aes(x=1,y=boot_tmleVar),fill="steelblue",width=0.70) +
  labs(x="Variance of TMLE Estimator (log scale)") +
  theme_bw() +
  theme(axis.title=element_text(size=13),axis.title.y=element_blank()) +
  coord_cartesian(xlim = c(0,2)) +
  scale_x_continuous(breaks=NULL,minor_breaks=NULL,labels=NULL) +
  scale_y_log10()

#TMLE CI distribution
plot_tmleCI <- ggplot() +
  geom_boxplot(data=boot_tmleCI,aes(x=1,y=V1),fill="steelblue",width=0.75) +
  geom_boxplot(data=boot_tmleCI,aes(x=2,y=V2),fill="tomato",width=0.75) +
  labs(y="Distribution of 95% CI Lower and Upper Bounds (log scale)") +
  theme_bw() +
  theme(axis.title=element_text(size=13),axis.title.y=element_blank()) +
  scale_x_continuous(breaks=NULL,minor_breaks=NULL,labels=NULL) +
  scale_y_log10() +
  coord_flip()

#TMLE epsilon convergence
plot_epsilon <- ggplot(data=boot_epsilon) +
  geom_point(aes(x=variable,y=value,colour=as.factor(sample)),alpha=0.5) +
  geom_line(aes(x=variable,y=value,group=sample,colour=as.factor(sample)),alpha=0.5) +
  guides(colour=FALSE) +
  labs(y="Epsilon") +
  theme_bw() +
  theme(axis.title=element_text(size=13),axis.text=element_text(size=13),axis.title.x=element_blank()) +
  scale_y_log10()

grid.arrange(plot_tmle,plot_tmleVar,plot_tmleCI,plot_epsilon,ncol=2)

#EDA Plots
data_raw <- read.csv("dat.csv")

#LOESS regression Y on A
eda_1 <- ggplot(data=data_raw) +
  geom_point(aes(x=daily_mod.vigPA.min.,y=WMHvolume.voxels.)) +
  geom_smooth(method="loess",aes(x=daily_mod.vigPA.min.,y=WMHvolume.voxels.),colour="steelblue") +
  labs(x="Daily Moderate/Vigorous Physical Activity",y="White Matter Leison Volume (log scale)") +
  theme_bw() +
  theme(axis.title=element_text(size=13),axis.text=element_text(size=13)) +
  scale_y_log10()
  
#Boxplots
eda_2 <- ggplot(data=melt(data[,c("A","Y")],id.vars="A"),aes(as.factor(A),value,group=as.factor(A),fill=as.factor(A))) +
  geom_boxplot() +
  scale_fill_manual(values=c("tomato","steelblue")) +
  guides(fill=FALSE) +
  labs(x="A",y="") +
  theme_bw() +
  theme(axis.title=element_text(size=13),axis.text=element_text(size=13),axis.title.y=element_blank()) +
  scale_y_log10()
  
grid.arrange(eda_1,eda_2,ncol=2)

#summary statistics of all estimators
#IPTW
mean(boot_iptw,na.rm=T) + (1.96 * sd(boot_iptw,na.rm=T))
mean(boot_iptw,na.rm=T) - (1.96 * sd(boot_iptw,na.rm=T))
quantile(boot_iptw,probs=c(.025,.975),na.rm=T)

#Stabilized IPTW
mean(boot_iptwS,na.rm=T) + (1.96 * sd(boot_iptwS,na.rm=T))
mean(boot_iptwS,na.rm=T) - (1.96 * sd(boot_iptwS,na.rm=T))
quantile(boot_iptwS,probs=c(.025,.975),na.rm=T)

#Simple Substitution Estimator
mean(boot_sub,na.rm=T) + (1.96 * sd(boot_sub,na.rm=T))
mean(boot_sub,na.rm=T) - (1.96 * sd(boot_sub,na.rm=T))
quantile(boot_sub,probs=c(.025,.975),na.rm=T)

#TMLE estimator
quantile(boot_tmle,probs=c(.025,.975),na.rm=T)