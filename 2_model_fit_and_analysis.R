##' This script 
##' - loads the pre-processed behavioral data 
##' - loads the DDM fits of both the training and the testing phase
##' - Fits the prior belief parameter (Vs) to the feedback received in the training
##' - Generates model predictions of RT/accuracy and confidence in the testing phase
##' - Computes stat tests on the predictions and fitted parameters

curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)
source("1_preprocessing.R")
library(reshape)
library(effects)
library(lmerTest)
library(scales)
library(DEoptim)
# FIX (reproducibility report, Sept 2026): see 1_preprocessing.R - 'prob' is
# archived on CRAN and unused by this pipeline; load it only if available.
# if (requireNamespace("prob", quietly = TRUE)) library(prob)
library(car)
library(MALDIquant)
library(Rcpp)
sourceCpp("DDM_with_confidence_slow.cpp")
sourceCpp("DDM_with_confidence_slow_fullconfRT.cpp")
source("quantile_fit_DDM.R")
source("build_hm.R")


stat_test <- F
# Global parameters --------------------------------------------------------------
## Heat map resolution
dt <- .001; ev_bound <- .5; ev_window <- dt*10; upperRT <- 5
ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
timesteps <- upperRT/dt

## List of heat maps
v_s_min <- dt; v_s_max <- .5; step <- dt
v_s_all <- seq(v_s_min,v_s_max,step)

nsim <- 1 # Number of times the dataset is simulated in model predictions
nrepeat <- 24 # Number of times each v_s is estimated
# EXP 1 -------------------------------------------------------------------
## Data load =====
subs_1 <- sort(unique(Training_exp1$sub)); Nsub_1 <- length(subs_1) 
cond_1 <- sort(unique(Training_exp1$fbcond)); Ncond_1 <- length(cond_1)
trialdifflevel <- sort(unique(Training_exp1$trialdifflevel));Ndiff <- length(trialdifflevel)



# Load fitted DDM parameters in the training phase 
bound_train <- matrix(NA,Nsub_1,Ncond_1)
ter_train <- matrix(NA,Nsub_1,Ncond_1) 
# 3 difficulty levels so 3 drift rates to fit
v_train <- matrix(NA,Nsub_1,Ncond_1)
v2_train <- matrix(NA,Nsub_1,Ncond_1) 
v3_train <- matrix(NA,Nsub_1,Ncond_1)
resid_train <- matrix(NA,Nsub_1,Ncond_1)
for (i in 1:Nsub_1) {
  tempAll <- subset(Training_exp1,sub==subs_1[i])
  for (cond in 1:Ncond_1) {
    print(paste('Running participant',i,'from',Nsub_1,"condition",cond))
    tempDat <- subset(tempAll,fbcond==cond_1[cond])
    file_name <- paste0('Fits/Exp1/Train/trainfit',cond_1[cond],subs_1[i],'.Rdata')
    if(file.exists(file_name)){
      load(file_name)
    }
    else{ #if not, fit the model
      optimal_params <- DEoptim(quantile_optim_DDM, # function to optimize
                                # bound,ter,z,ntrials,sigma,dt,t2time,vratio,drift(s)
                                lower = c( 0, 0, 0, 5000, .1, dt, 0,1,0,0,0), 
                                upper = c(.2, 2, 0, 5000, .1, dt, 0,1,.5,.5,.5), 
                                observations = tempDat,
                                control=c(itermax=1000,steptol=100,reltol=dt,NP=50), 
                                returnFit = 1)
      results <- summary(optimal_params)
      #save individual results
      save(results, file=file_name)
    }
    bound_train[i,cond] <- results$optim$bestmem[1]
    ter_train[i,cond] <- results$optim$bestmem[2]
    v_train[i,cond] <- results$optim$bestmem[9]
    v2_train[i,cond] <- results$optim$bestmem[10]
    v3_train[i,cond] <- results$optim$bestmem[11]
    resid_train[i,cond] <- results$optim$bestval
  }  
}

# Load fitted DDM parameters in the testing phase 
bound <- matrix(NA,Nsub_1,Ncond_1)
ter <- matrix(NA,Nsub_1,Ncond_1)
v <- matrix(NA,Nsub_1,Ncond_1)
v2 <- matrix(NA,Nsub_1,Ncond_1)
v3 <- matrix(NA,Nsub_1,Ncond_1) 
resid <- matrix(NA,Nsub_1,Ncond_1)
for(i in 1:Nsub_1){
  for(cond in 1:Ncond_1){
    print(paste('Running participant',i,'from',Nsub_1,"condition",cond))
    file_name <- paste0('Fits/Exp1/Test/testfit',cond_1[cond],subs_1[i],'.Rdata')
    if(file.exists(file_name)){
      load(file_name)
    }
    else{ #if not, fit the model
      optimal_params <- DEoptim(quantile_optim_DDM, # function to optimize
                                # bound,ter,z,ntrials,sigma,dt,t2time,vratio,drift(s)
                                lower = c( 0, 0, 0, 5000, .1, dt, 0,1,0,0,0), 
                                upper = c(.2, 2, 0, 5000, .1, dt, 0,1,.5,.5,.5),
                                observations = tempDat,
                                control=c(itermax=1000,steptol=100,reltol=dt,NP=50), 
                                returnFit = 1)
      results <- summary(optimal_params)
      #save individual results
      save(results, file=file_name)
    }
    bound[i,cond] <- results$optim$bestmem[1]
    ter[i,cond] <- results$optim$bestmem[2]
    v[i,cond] <- results$optim$bestmem[9]
    v2[i,cond] <-   results$optim$bestmem[10]
    v3[i,cond] <-   results$optim$bestmem[11]
    resid[i,cond] <- results$optim$bestval
  }
}

param_ddm_test_exp1 <- data.frame(drift = c(v,v2,v3),bound=rep(bound,Ndiff),ter=rep(ter,Ndiff),
                      sub=rep(subs_1,Ndiff*Ncond_1),
                      condition=rep(cond_1,each=Nsub_1,length.out=Nsub_1*Ncond_1*Ndiff),
                      difflevel=rep(trialdifflevel,each=Nsub_1*Ncond_1),exp=1,resid=rep(resid,Ndiff))
## Fit subjective drift ====
if (!(file.exists("Data/Aggregated/cost_vs_exp1.csv"))) {
  s <- 1; cond <- 1
  while (s <= Nsub_1) {
    while (cond <= Ncond_1) {
      print(paste("Running participant",s,"of",Nsub_1,"condition",cond))
      cost_conf <- matrix(NA,nrow=nrepeat,ncol=length(v_s_all))
      tempDat <- subset(Training_exp1,sub==subs_1[s]&fbcond==cond_1[cond])
      tempDat_test <- subset(Data_exp1,sub==subs_1[s]&fbcond==cond_1[cond])
      ntrial <- dim(tempDat_test)[1]
      ntrial_train <- dim(tempDat)[1]
      temp_par <- c(bound_train[s,cond],ter_train[s,cond],0,nrepeat,
                    .1,dt,1,v_train[s,cond],v2_train[s,cond],v3_train[s,cond])
      
      # First generate trials from estimated DDM parameters in the training phase
      temp <- quantile_optim_DDM_fullconfRT(temp_par,observations=tempDat_test,returnFit=0)
      
      # match to the heatmap
      temp$closest_evdnc2 <- match.closest(temp$evidence2,ev_mapping)
      temp$rt2[temp$rt2>5] <- 5 #heatmap doesn't go higher (5 seconds)
      temp$rt2 <- temp$rt2*timesteps/5 #scale with the heatmap
      
      # We compute the cost function several times to take noise into account
      temp$Nrep <- rep(1:nrepeat,each=ntrial/Ndiff,length.out=ntrial*nrepeat)        
      
      bar <- txtProgressBar(0,length(v_s_all),style=3,char="#")
      for (d in 1:length(v_s_all)) {
        hm_up <- build_hm(v_s_all[d])
        hm_low <- 1-hm_up
        hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
        
        temp[temp$resp==1,]$cj <- hmvec_up[(temp[temp$resp==1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==1,]$rt2)]
        temp[temp$resp==-1,]$cj <- hmvec_low[(temp[temp$resp==-1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==-1,]$rt2)]
        
        
        for (i in 1:nrepeat) {
          pred_sample <- subset(temp,Nrep==i)
          #' Randomly sample 40 trials for each difficulty level to get the same number of 
          #' trials as in the training
          pred_sample <- do.call(rbind,
                                 lapply(split(pred_sample, pred_sample$drift),
                                        function(x) x[sample(nrow(x), ntrial_train/Ndiff),]))
          diff <- sum((tempDat$fb - pred_sample$cj)^2)
          cost_conf[i,d] <- diff
        }
        setTxtProgressBar(bar,d)
      }
      temp_df <- data.frame(cost = as.vector(cost_conf),
                            Nrep = rep(1:nrepeat,length(v_s_all)),
                            Vs = rep(v_s_all,each=nrepeat),
                            sub = subs_1[s], fbcond = cond_1[cond])
      if (s==1 & cond==1) {
        cost_df <- temp_df
      }else{
        cost_df <- rbind(cost_df,temp_df)
      }
      cond <- cond + 1
    }
    s <- s + 1
    cond <- 1
  }
  write.csv(cost_df,file="Data/Aggregated/cost_vs_exp1.csv")
}else{
  cost_df <- read.csv("Data/Aggregated/cost_vs_exp1.csv")
}

means_fullconfRT <- with(cost_df,aggregate(cost,by=list(Vs=Vs,sub=sub,fbcond=fbcond),mean))
means_fullconfRT <- cast(means_fullconfRT,fbcond+sub~Vs)

# Smooth over subjective v_s_all
vs_smooth1 <- sapply(seq(nrow(means_fullconfRT)), function(i) {
  j <- as.numeric(means_fullconfRT[i,3:dim(means_fullconfRT)[2]])
  j <- lowess(j,f=.05)
  j <- which.min(j$y)
  c(v_s_all[j])
})

param_train_exp1 <- data.frame(Vs=vs_smooth1,bound = c(bound_train), ter = c(ter_train),Vo=c(v_train),
                 sub=rep(subs_1,Ncond_1),condition=rep(cond_1,each=Nsub_1))

## Generate model prediction ====
if (file.exists("Data/Aggregated/model_prediction_exp1.csv")) {
  Simuls <- read.csv("Data/Aggregated/model_prediction_exp1.csv")
}else{
  rm(Simuls)
  for(i in 1:Nsub_1){
    print(paste('simulating',i,'from',Nsub_1))
    for(c in 1:Ncond_1){
      temp_vs <- subset(param_train_exp1,condition==cond_1[c]&sub==subs_1[i])$Vs
      tempDat <- subset(Data_exp1,fbcond==cond_1[c]&sub==subs_1[i])
      hm_up <- build_hm(temp_vs)
      hm_low <- 1-hm_up
      hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
      temp <- quantile_optim_DDM_fullconfRT(c(bound[i,c],ter[i,c],0,nsim,.1,dt,1,v[i,c],v2[i,c],v3[i,c]),tempDat,0)
      
      #match to the heatmap
      temp$closest_evdnc2 <- match.closest(temp$evidence2,ev_mapping)
      
      temp$rt2[temp$rt2>5] <- 5 #heatmap doesn't go higher
      temp$rt2 <- temp$rt2*timesteps/5 #scale with the heatmap, between 0 and 2000
      temp[temp$resp==1,]$cj <- hmvec_up[(temp[temp$resp==1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==1,]$rt2)]
      temp[temp$resp==-1,]$cj <- hmvec_low[(temp[temp$resp==-1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==-1,]$rt2)]
      
      if(!exists('Simuls')){ Simuls <- cbind(temp,cond_1[c],subs_1[i])
      }else{ Simuls <- rbind(Simuls,cbind(temp,cond_1[c],subs_1[i]))
      }
    }
  }
  Simuls <- data.frame(Simuls);names(Simuls) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','closest_evdnc2','condition','sub')
  
  difflevels <- sort(unique(Data_exp1$trialdifflevel))
  Simuls$trialdifflevel <- 0
  for (i in 1:Nsub_1) {
    for(d in 1:length(difflevels)) Simuls$trialdifflevel[Simuls$sub==subs_1[i] & Simuls$drift %in% c(unique(subset(Simuls,sub==subs_1[i])$drift)[d],unique(subset(Simuls,sub==subs_1[i])$drift)[d+3],unique(subset(Simuls,sub==subs_1[i])$drift)[d+6])] <- difflevels[d] #recode drift to trialdifflevelerence
  }
  write.csv(Simuls,file = "Data/Aggregated/model_prediction_exp1.csv")
}
# EXP 2 -------------------------------------------------------------------
## Data Load ====
subs_2 <- sort(unique(Training_exp2$sub)); Nsub_2 <- length(subs_2)
cond_2 <- sort(unique(Training_exp2$traindiffcond)); Ncond_2 <- length(cond_2)
trialdifflevel <- sort(unique(Training_exp2$trialdifflevel));Ndiff <- length(trialdifflevel)


#Load fitted train DDM parameters
bound_train <- matrix(NA,Nsub_2,Ncond_2)
ter_train <- matrix(NA,Nsub_2,Ncond_2) 
v_train <- matrix(NA,Nsub_2,Ncond_2) # 1 difficulty per condition so only 1 drift rate
resid_train <- matrix(NA,Nsub_2,Ncond_2)

for (i in 1:Nsub_2) {
  tempAll <- subset(Training_exp2,sub==subs_2[i])
  for (c in 1:Ncond_2) {
    tempDat <- subset(tempAll,traindiffcond==cond_2[c])
    file_name <- paste0('Fits/Exp2/Train/trainfit',cond_2[c],subs_2[i],'.Rdata')
    if(file.exists(file_name)){
      load(file_name)
    }
    else{ #if not, fit the model
      optimal_params <- DEoptim(quantile_optim_DDM, # function to optimize
                                # a,ter,z,ntrials,sigma,dt,t2time,vratio,alpha,beta,v
                                lower = c( 0, 0, 0, 5000, .1, .001, 0,   1,0), 
                                upper = c(.2, 2, 0, 5000, .1, .001, 0, 1,.5),
                                observations = tempDat,returnFit = 1, binning = F,
                                control=c(itermax=1000,steptol=100,reltol=dt,NP=30))
      results <- summary(optimal_params)
      #save individual results
      save(results, file=file_name)
    }
    bound_train[i,c] <- results$optim$bestmem[1]
    ter_train[i,c] <- results$optim$bestmem[2]
    v_train[i,c] <- results$optim$bestmem[9]
    resid_train[i,c] <- results$optim$bestval
  }
}

# Load fitted test DDM parameters
bound <- matrix(NA,Nsub_2,Ncond_2)
ter <- matrix(NA,Nsub_2,Ncond_2)
v <- matrix(NA,Nsub_2,Ncond_2)
v2 <- matrix(NA,Nsub_2,Ncond_2)
v3 <- matrix(NA,Nsub_2,Ncond_2)
resid <- matrix(NA,Nsub_2,Ncond_2)
for(i in 1:Nsub_2){
  for(c in 1:Ncond_2){
    print(paste('Running participant',i,'from',Nsub_2,"condition",c))
    file_name <- paste0('Fits/Exp2/Test/testfit',cond_2[c],subs_2[i],'.Rdata')
    if(file.exists(file_name)){
      load(file_name)
    }
    else{ #if not, fit the model
      optimal_params <- DEoptim(quantile_optim_DDM, # function to optimize
                                # a,ter,z,ntrials,sigma,dt,t2time,vratio,alpha,beta,v
                                lower = c( 0, 0, 0, 5000, .1, .001, 0,   1,0,0,0), 
                                upper = c(.2, 2, 0, 5000, .1, .001, 0, 1,.5,.5,.5),
                                observations = tempDat,returnFit = 1, binning = F,
                                control=c(itermax=1000,steptol=100,reltol=dt,NP=50))
      results <- summary(optimal_params)
      #save individual results
      save(results, file=file_name)
    }
    bound[i,c] <- results$optim$bestmem[1]
    ter[i,c] <- results$optim$bestmem[2]
    v[i,c] <- results$optim$bestmem[9]
    v2[i,c] <-   results$optim$bestmem[10]
    v3[i,c] <-   results$optim$bestmem[11]
    resid[i,c] <- results$optim$bestval
  }
}

param_ddm_test_exp2 <- data.frame(drift = c(v,v2,v3),bound=rep(bound,Ndiff),ter=rep(ter,Ndiff),
                     sub=rep(subs_2,Ndiff*Ncond_2),
                     condition=rep(cond_2,each=Nsub_2,length.out=Nsub_2*Ncond_2*Ndiff),
                     difflevel=rep(trialdifflevel,each=Nsub_2*Ncond_2),exp=2,resid=rep(resid,Ndiff))
## Fit subjective drift  ====
Ndiff <- 1 # Only one difficulty level in the training phase
if (!(file.exists("Data/Aggregated/cost_vs_exp2.csv"))) {
  s <- 1; cond <- 1
  while (s <= Nsub_2) {
    while (cond <= Ncond_2) {
      print(paste("Running participant",s,"of",Nsub_2,"condition",cond))
      cost_conf <- matrix(NA,nrow=nrepeat,ncol=length(v_s_all))
      
      tempDat <- subset(Training_exp2,sub==subs_2[s]&traindiffcond==cond_2[cond])
      tempDat_test <- subset(Data_exp2,sub==subs_2[s]&traindiffcond==cond_2[cond])
      
      ntrial_train <- dim(tempDat)[1]
      
      temp_par <- c(bound_train[s,cond],ter_train[s,cond],0,nrepeat*2,.1,dt,1,v_train[s,cond],
                    1,1)
      
      # First generate trials from estimated DDM parameters in the training phase
      temp <- quantile_optim_DDM_fullconfRT(temp_par,observations=tempDat_test,returnFit=0)
      temp <- subset(temp,drift==v_train[s,cond])
      
      ntrial <- dim(temp)[1]/nrepeat
      
      # match to the heatmap
      temp$closest_evdnc2 <- match.closest(temp$evidence2,ev_mapping)
      temp$rt2[temp$rt2>5] <- 5 #heatmap doesn't go higher (5 seconds)
      temp$rt2 <- temp$rt2*timesteps/5 #scale with the heatmap
      
      # We compute the cost function several times to take noise into account
      temp$Nrep <- rep(1:nrepeat,each=ntrial/Ndiff,length.out=ntrial*nrepeat)        
      
      bar <- txtProgressBar(0,length(v_s_all),style=3,char="#")
      for (d in 1:length(v_s_all)) {
        hm_up <- build_hm(v_s_all[d])
        hm_low <- 1-hm_up
        hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
        
        temp[temp$resp==1,]$cj <- hmvec_up[(temp[temp$resp==1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==1,]$rt2)]
        temp[temp$resp==-1,]$cj <- hmvec_low[(temp[temp$resp==-1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==-1,]$rt2)]
        
        
        for (i in 1:nrepeat) {
          pred_sample <- subset(temp,Nrep==i)
          #' Randomly sample 120 trials to get the same number of trials as in the training
          pred_sample <- do.call(rbind,
                                 lapply(split(pred_sample, pred_sample$drift),
                                        function(x) x[sample(nrow(x), ntrial_train/Ndiff),]))
          diff <- sum((tempDat$cor - pred_sample$cj)^2)
          cost_conf[i,d] <- diff
        }
        setTxtProgressBar(bar,d)
      }
      temp_df_cor <- data.frame(cost = as.vector(cost_conf),
                                Nrep = rep(1:nrepeat,length(v_s_all)),
                                Vs = rep(v_s_all,each=nrepeat),
                                sub = subs_2[s], traindiffcond = cond_2[cond])
      if (s==1 & cond==1) {
        cost_df <- temp_df_cor
      }else{
        cost_df <- rbind(cost_df,temp_df_cor)
      }
      write.csv(cost_df,file="Data/Aggregated/cost_vs_exp2.csv")
      cond <- cond + 1
    }
    s <- s + 1
    cond <- 1
  }
  write.csv(cost_df,file="Data/Aggregated/cost_vs_exp2.csv")
}else{
  cost_df <- read.csv("Data/Aggregated/cost_vs_exp2.csv")
}

means_fullconfRT <- with(cost_df,aggregate(cost,by=list(Vs=Vs,sub=sub,traindiffcond=traindiffcond),mean))
means_fullconfRT <- cast(means_fullconfRT,traindiffcond+sub~Vs)

vs_smooth2 <- sapply(seq(nrow(means_fullconfRT)), function(i) {
  j <- as.numeric(means_fullconfRT[i,3:dim(means_fullconfRT)[2]])
  j <- lowess(j,f=.05)
  j <- which.min(j$y)
  c(v_s_all[j])
})

param_train_exp2 <- data.frame(Vs=vs_smooth2,bound = c(bound_train), ter = c(ter_train),
                  Vo=c(v_train),sub=rep(subs_2,Ncond_2),condition=rep(cond_2,each=Nsub_2))
## Generate model prediction ====
if (file.exists("Data/Aggregated/model_prediction_exp2.csv")) {
  Simuls2 <- read.csv("Data/Aggregated/model_prediction_exp2.csv")
}else{
  rm(Simuls2)
  for(i in 1:Nsub_2){
    print(paste('simulating',i,'from',Nsub_2))
    for(c in 1:Ncond_2){
      temp_vs <- subset(param_train_exp2,condition==cond_2[c]&sub==subs_2[i])$Vs
      tempDat <- subset(Data_exp2, sub==subs_2[i] & traindiffcond==cond_2[c])
      hm_up <- build_hm(temp_vs)
      hm_low <- 1-hm_up
      hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
      temp <- quantile_optim_DDM_fullconfRT(c(bound[i,c],ter[i,c],0,nsim,.1,
                                                dt,1,v[i,c],v2[i,c],v3[i,c]),tempDat,0)
      
      #match to the heatmap
      temp$closest_evdnc2 <- match.closest(temp$evidence2,ev_mapping)
      
      temp$rt2[temp$rt2>5] <- 5 #heatmap doesn't go higher
      temp$rt2 <- temp$rt2*timesteps/5 #scale with the heatmap, between 0 and 2000
      
      temp[temp$resp==1,]$cj <- hmvec_up[(temp[temp$resp==1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==1,]$rt2)]
      temp[temp$resp==-1,]$cj <- hmvec_low[(temp[temp$resp==-1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==-1,]$rt2)]
      
      if(!exists('Simuls2')){ Simuls2 <- cbind(temp,cond_2[c],subs_2[i])
      }else{ Simuls2 <- rbind(Simuls2,cbind(temp,cond_2[c],subs_2[i]))
      }
    }
  }
  Simuls2 <- data.frame(Simuls2);names(Simuls2) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','closest_evdnc2','condition','sub')
  
  difflevels <- sort(unique(Data_exp2$trialdifflevel))
  Simuls2$trialdifflevel <- 0
  for (i in 1:Nsub_2) {
    for(d in 1:length(difflevels)) Simuls2$trialdifflevel[Simuls2$sub==subs_2[i] & Simuls2$drift %in% c(unique(subset(Simuls2,sub==subs_2[i])$drift)[d],unique(subset(Simuls2,sub==subs_2[i])$drift)[d+3],unique(subset(Simuls2,sub==subs_2[i])$drift)[d+6])] <- difflevels[d] #recode drift to trialdifflevelerence
  }
  write.csv(Simuls2,file = "Data/Aggregated/model_prediction_exp2.csv")
}
# Stat tests ------------------------------------------------------------
# Exp1 ====
if (stat_test) {
  #DDM train
  df$sub <- as.factor(df$sub)
  m <- lmer(Vs ~ condition + (1|sub),data=df); anova(m);
  m <- lmer(Vo ~ condition + (1|sub), data = df); anova(m)
  m <- lmer(bound ~ condition + (1|sub),data=df); anova(m);
  m <- lmer(ter ~ condition + (1|sub),data=df); anova(m);
  
  #DDM test
  param_ddm_test_exp1$sub <- as.factor(param_ddm_test_exp1$sub)
  test_bound1 <- with(param_ddm_test_exp1,aggregate(bound,by=list(condition=condition,sub=sub),mean))
  test_ter1 <- with(param_ddm_test_exp1,aggregate(ter,by=list(condition=condition,sub=sub),mean))
  m <- lmer(x ~ condition + (1|sub),data=test_bound1); anova(m);
  m <- lmer(x ~ condition + (1|sub),data=test_ter1); anova(m);
  m <- lmer(drift ~ condition*difflevel + (1|sub),data=param_ddm_test_exp1); anova(m);
}
## Exp2 ====
if (stat_test) {
  #DDM train
  df2$sub <- as.factor(df2$sub)
  m <- lmer(Vs ~condition + (1|sub), data = df2); anova(m)
  m <- lmer(Vo ~ condition + (1|sub),data=df2); anova(m);
  m <- lmer(bound ~ condition + (1|sub),data=df2); anova(m);
  m <- lmer(ter ~ condition + (1|sub),data=df2); anova(m)
  
  #DDM test
  param_ddm_test_exp2$sub <- as.factor(param_ddm_test_exp2$sub)
  test_bound2 <- with(param_ddm_test_exp2,aggregate(bound,by=list(condition=condition,sub=sub),mean))
  test_ter2 <- with(param_ddm_test_exp2,aggregate(ter,by=list(condition=condition,sub=sub),mean))
  m <- lmer(x ~ condition + (1|sub),data=test_bound2); anova(m);
  m <- lmer(x ~ condition + (1|sub),data=test_ter2); anova(m);
  m <- lmer(drift ~ condition*difflevel + (condition|sub),data=param_ddm_test_exp2); anova(m);
  
  m <- lmer(cj ~ condition*trialdifflevel + (condition+trialdifflevel|sub), data = Simuls, 
            control = lmerControl(optimizer='bobyqa')); 
  anova(m)
  m2 <- lmer(cj ~ condition*trialdifflevel + (condition+trialdifflevel|sub), data = Simuls2, 
            control = lmerControl(optimizer='bobyqa')); 
  anova(m2)
  

}
# Save aggregated data for plotting purpose --------------------------------------------------------------------
save(param_ddm_test_exp1,file="Data/Aggregated/param_ddm_test_exp1.Rdata")
save(param_ddm_test_exp2,file="Data/Aggregated/param_ddm_test_exp2.Rdata")
save(param_train_exp1,file="Data/Aggregated/param_train_exp1.Rdata")
save(param_train_exp2,file="Data/Aggregated/param_train_exp2.Rdata")


