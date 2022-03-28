##' This script 
##' - loads the pre-processed behavioral data 
##' - loads the DDM fits of both the training and the testing phase
##' - Fits the prior belief parameter (Vs) to the feedback received in the training
##' - Generates model predictions of RT/accuracy and confidence in the testing phase
##' - Computes stat tests on the predictions and fitted parameters
##' - Plots the results

rm(list=ls())
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)
library(reshape)
library(effects)
library(lmerTest)
library(scales)
library(DEoptim)
library(prob)
library(car)
library(myPackage)
library(MALDIquant)
library(emmeans)
library(multcomp) # cld for post hoc test
# source("1_preprocessing.R")
## Transparent colors, Mark Gardener 2015, www.dataanalytics.org.uk
transp <- function(color, percent = 50, name = NULL) {
  #   color = color name
  #   percent = % transparency
  #   name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  transp <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                max = 255,
                alpha = (100 - percent) * 255 / 100,
                names = name)
  
  ## Save the color
  invisible(transp)
}

error.bar <- function(x, y, upper, lower=upper, length=0,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
cexkl <- 1.5;cexgr <- 2;lwdgr <- 3; lwddat <- 2
cex_lab <- 3; cex_legend <- 3; cex_title <- 2.5 
windowsFonts(A = windowsFont("Calibri")) 
par(family="A",font.main = 2, cex.main = cex_title)

# Global parameters --------------------------------------------------------------
## Heat map resolution
dt <- .001; ev_bound <- .5; ev_window <- dt*10; upperRT <- 5
ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
timesteps <- upperRT/dt

## List of heat maps
v_min <- .001; v_max <- .5; step <- .001
drifts <- seq(v_min,v_max,step)

nsim <- 500 # per drift/cond/participant
nrepeat <- 24 # Vs fitting
ntrials <- 5000 #model prediction
# EXP 1 -------------------------------------------------------------------
## Data load ====
go_to("results")
Data1 <- read.csv('data_exp1.csv')
Data2 <- read.csv('data_exp2.csv')
Data1_train <- read.csv("data_exp1_training.csv")
Data2_train <- read.csv("data_exp2_training.csv")


subs1 <- sort(unique(Data1_train$sub)); N1 <- length(subs1) 
cond_1 <- sort(unique(Data1_train$selfconf)); Ncond <- length(cond_1)
coh <- sort(unique(Data1_train$coh));Ndiff <- length(coh)


go_to("fits")
#Load fitted DDM parameters in the training phase + median confidence RT
bound_train <- matrix(NA,N1,Ncond);v_train <- matrix(NA,N1,Ncond);
ter_train <- matrix(NA,N1,Ncond); resid_train <- matrix(NA,N1,Ncond)
resid <- matrix(NA,N1,Ncond)
v2_train <- matrix(NA,N1,Ncond); v3_train <- matrix(NA,N1,Ncond)
conf_rt <- matrix(NA,N1,Ncond)
for (i in 1:N1) {
  tempAll <- subset(Data1_train,sub==subs1[i])
  for (cond in 1:Ncond) {
    print(paste('Running participant',i,'from',N1,"condition",cond))
    tempDat <- subset(tempAll,selfconf==cond_1[cond])
    file_name <- paste0('exp1/train/trainfit',cond_1[cond],subs1[i],'.Rdata')
    if(file.exists(file_name)){
      load(file_name)
    }
    else{ #if not, fit the model
      optimal_params <- DEoptim(chi_square_optim_DDM, # function to optimize
                                lower = c( 0, 0, 0, 5000, .1, .0025, 0,1,0,0,0), # a,ter,z,ntrials,sigma,dt,t2time,vratio,alpha,beta,v
                                upper = c(.2, 2, 0, 5000, .1, .0025, 0,1,.5,.5,.5), # a,ter,z,ntrials,sigma,dt,t2time,vratio,alpha,beta,v
                                observations = tempDat,
                                control=c(itermax=1000,steptol=100,reltol=.001,NP=50), 
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

#Load fitted DDM parameters in the testing phase 
bound <- matrix(NA,N1,Ncond);v <- matrix(NA,N1,Ncond);ter <- matrix(NA,N1,Ncond)
conf_rt <- matrix(NA,N1,Ncond); resid <- matrix(NA,N1,Ncond)
#Adjust the number of drift parameters to the model loaded
v2 <- matrix(NA,N1,Ncond);v3 <- matrix(NA,N1,Ncond) 
for(i in 1:N1){
  for(cond in 1:Ncond){
    print(paste('Running participant',i,'from',N1,"condition",cond))
    file_name <- paste0('exp1/test/testfit',cond_1[cond],subs1[i],'.Rdata')
    if(file.exists(file_name)){
      load(file_name)
    }
    else{ #if not, fit the model
      optimal_params <- DEoptim(chi_square_optim_DDM, # function to optimize
                                lower = c( 0, 0, 0, 5000, .1, .0025, 0,1,0,0,0), # a,ter,z,ntrials,sigma,dt,t2time,vratio,alpha,beta,v
                                upper = c(.2, 2, 0, 5000, .1, .0025, 0,1,.5,.5,.5), # a,ter,z,ntrials,sigma,dt,t2time,vratio,alpha,beta,v
                                observations = tempDat,
                                control=c(itermax=1000,steptol=100,reltol=.001,NP=50), 
                                returnFit = 1)
      results <- summary(optimal_params)
      #save individual results
      save(results, file=file_name)
    }
    bound[i,cond] <- results$optim$bestmem[1]
    ter[i,cond] <- results$optim$bestmem[2]
    conf_rt[i,cond] <- results$optim$bestmem[7]
    v[i,cond] <- results$optim$bestmem[9]
    v2[i,cond] <-   results$optim$bestmem[10]
    v3[i,cond] <-   results$optim$bestmem[11]
    resid[i,cond] <- results$optim$bestval
  }
}
param_1 <- data.frame(drift = c(v,v2,v3),bound=rep(bound,Ndiff),ter=rep(ter,Ndiff),
                      sub=rep(subs1,Ndiff*Ncond),
                      condition=rep(cond_1,each=N1,length.out=N1*Ncond*Ndiff),
                      difflevel=rep(coh,each=N1*Ncond),exp=1,resid=rep(resid,Ndiff))

##' Aggregate train and test
bounds1 <- data.frame(bound = c(bound_train,bound),
                     phase = rep(c("train","main"),each=length(bound)), 
                     sub = rep(subs1,Ncond*2),
                     condition = rep(cond_1,each = N1,length.out=N1*Ncond*2))
ters1 <- data.frame(ter = c(ter_train,ter),
                   phase = rep(c("train","main"),each=length(ter)), 
                   sub = rep(subs1,Ncond*2),
                   condition = rep(cond_1,each = N1,length.out=N1*Ncond*2))
vs <- data.frame(v = c(v_train,v2_train,v3_train,v,v2,v3),
                 phase = rep(c("train","main"),each=length(v)*Ndiff), 
                 sub = rep(subs1,Ncond*2*Ndiff),
                 condition = rep(cond_1,each = N1,length.out=N1*Ncond*2*Ndiff),
                 difficulty=rep(coh,each=N1*Ncond,length.out=N1*Ncond*2))

## Fit subjective drift ====
go_to("results")
if (!(file.exists("cost_vs_exp1.csv"))) {
  means <- matrix(NA,nrow=Ncond*N1,ncol=length(drifts)) 
  stds <- matrix(NA,nrow=Ncond*N1,ncol=length(drifts))
  s <- 1; cond <- 1
  while (s <= N1) {
    while (cond <= Ncond) {
      cost_conf <- matrix(NA,nrow=nrepeat,ncol=length(drifts))
      print(paste("Running participant",s,"of",N1,"condition",cond))
      tempDat <- subset(Data1_train,sub==subs1[s]&selfconf==cond_1[cond])
      tempDat_test <- subset(Data1,sub==subs1[s]&selfconf==cond_1[cond])
      ntrial <- dim(tempDat_test)[1]
      ntrial_train <- dim(tempDat)[1]
      temp_par <- c(bound_train[s,cond],ter_train[s,cond],0,nrepeat,
                    .1,.001,1,v_train[s,cond],v2_train[s,cond],v3_train[s,cond])
      
      # First generate trials from estimated DDM parameters in the training phase
      temp <- chi_square_optim_DDM_fullconfRT(temp_par,observations=tempDat_test,returnFit=0)
      
      # match to the heatmap
      temp$closest_evdnc2 <- match.closest(temp$evidence2,ev_mapping)
      temp$temprt2 <- temp$rt2;
      temp$temprt2[temp$temprt2>5] <- 5 #heatmap doesn't go higher (5 seconds)
      temp$temprt2 <- temp$temprt2*timesteps/5 #scale with the heatmap
      
      # We compute the cost function several times to take noise into account
      temp$Nrep <- rep(1:nrepeat,each=ntrial/Ndiff,length.out=ntrial*nrepeat)        
      
      bar <- txtProgressBar(0,length(drifts),style=3,char="#")
      for (d in 1:length(drifts)) {
        load(paste0("heatmaps/hm_",drifts[d],"_filled.Rdata"))
        hm_low <- output$lower; hm_up <- output$upper
        hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
        
        temp[temp$resp==1,]$cj <- hmvec_up[(temp[temp$resp==1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==1,]$temprt2)]
        temp[temp$resp==-1,]$cj <- hmvec_low[(temp[temp$resp==-1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==-1,]$temprt2)]
        
        
        for (i in 1:nrepeat) {
          pred_sample <- subset(temp,Nrep==i)
          #' Randomly sample 40 trials for each drift to get the same number of 
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
                            Nrep = rep(1:nrepeat,length(drifts)),
                            Vs = rep(drifts,each=nrepeat),
                            sub = subs1[s], selfconf = cond_1[cond])
      if (s==1 & cond==1) {
        cost_df <- temp_df
      }else{
        cost_df <- rbind(cost_df,temp_df)
      }
      means[s+N1*(cond-1),] <- colMeans(cost_conf)
      stds[s+N1*(cond-1),] <- colSds(cost_conf)
      cond <- cond + 1
    }
    s <- s + 1
    cond <- 1
  }
  # save(means,file="means_exp1_full.Rdata")
  # save(stds,file="stds_exp1_full.Rdata")
  write.csv(cost_df,file="cost_vs_exp1.csv")
}else{
  load("means_exp1_full.Rdata")
  load("stds_exp1_full.Rdata")
  cost_df <- read.csv("cost_vs_exp1.csv")
}

means_fullconfRT <- with(cost_df,aggregate(cost,by=list(Vs=Vs,sub=sub,selfconf=selfconf),mean))
means_fullconfRT <- cast(means_fullconfRT,selfconf+sub~Vs)

vs_smooth1 <- sapply(seq(nrow(means_fullconfRT)), function(i) {
  j <- as.numeric(means_fullconfRT[i,3:dim(means_fullconfRT)[2]])
  j <- lowess(j,f=.05)
  j <- which.min(j$y)
  c(drifts[j])
})

result <- sapply(seq(nrow(means)),function(i) {
  j <- which.min(means[i,])
  c(j)
})
Vs1 <- drifts[result]
Vs1_matrix <- matrix(Vs1,nrow=N1,ncol=Ncond)

df <- data.frame(Vs=vs_smooth1,bound = c(bound_train), ter = c(ter_train),Vo=c(v_train),
                 sub=rep(subs1,Ncond),condition=rep(cond_1,each=N1),Vs=Vs1)
## Generate model prediction ====
go_to("results")

if (file.exists("model_prediction_exp1.csv")) {
  Simuls <- read.table("model_prediction_exp1.csv")
}else{
  rm(Simuls)
  for(i in 1:N1){
    print(paste('simulating',i,'from',N1))
    for(c in 1:Ncond){
      load(paste0("heatmaps/hm_",Vs1_matrix[i,c],"_filled.Rdata"))
      hm_low <- output$lower; hm_up <- output$upper
      hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
      temp <- chi_square_optim_DDM(c(bound[i,c],ter[i,c],0,nsim,.1,.001,conf_rt[i,c],1,v[i,c],v2[i,c],v3[i,c]),NULL,0)
      
      #match to the heatmap
      temp$closest_evdnc2 <- match.closest(temp$evidence2,ev_mapping)
      temp$temprt2 <- temp$rt2;
      temp$temprt2[temp$temprt2>5] <- 5 #heatmap doesn't go higher
      temp$temprt2 <- temp$temprt2*timesteps/5 #scale with the heatmap, between 0 and 2000
      
      temp[temp$resp==1,]$cj <- hmvec_up[(temp[temp$resp==1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==1,]$temprt2)]
      temp[temp$resp==-1,]$cj <- hmvec_low[(temp[temp$resp==-1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==-1,]$temprt2)]
      
      temp$cj_cont <- temp$cj
      
      
      if(!exists('Simuls')){ Simuls <- cbind(temp,rep(cond_1[c],nsim),rep(subs1[i],nsim))
      }else{ Simuls <- rbind(Simuls,cbind(temp,rep(cond_1[c],nsim),rep(subs1[i],nsim)))
      }
    }
  }
  Simuls <- data.frame(Simuls);names(Simuls) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','closest_evdnc2',"temprt2",'cj_cont','condition','sub')
  
  coherences <- sort(unique(Data1$coh))
  Simuls$coh <- 0
  for (i in 1:N1) {
    for(d in 1:length(coherences)) Simuls$coh[Simuls$sub==subs1[i] & Simuls$drift %in% c(unique(subset(Simuls,sub==subs1[i])$drift)[d],unique(subset(Simuls,sub==subs1[i])$drift)[d+3],unique(subset(Simuls,sub==subs1[i])$drift)[d+6])] <- coherences[d] #recode drift to coherence
  }
  write.csv(Simuls,file = "model_prediction_exp1.csv")
}
# EXP 2 -------------------------------------------------------------------
## Data Load ====
subs_2 <- sort(unique(Data2_train$sub)); Nsub_2 <- length(subs_2)
cond_2 <- sort(unique(Data2_train$traindiffcond)); Ncond_2 <- length(cond_2)
coh <- sort(unique(Data2_train$coh));Ndiff <- length(coh)

#Load fitted train DDM parameters
bound_train <- matrix(NA,Nsub_2,Ncond_2);v_train <- matrix(NA,Nsub_2,Ncond_2);
ter_train <- matrix(NA,Nsub_2,Ncond_2); resid_train <- matrix(NA,Nsub_2,Ncond_2)
resid <- matrix(NA,Nsub_2,Ncond_2)
v2 <- matrix(NA,Nsub_2,Ncond_2); v3 <- matrix(NA,Nsub_2,Ncond_2)
conf_rt <- matrix(NA,Nsub_2,Ncond_2)

go_to("fits")
for (i in 1:Nsub_2) {
  tempAll <- subset(Data2_train,sub==subs_2[i])
  for (c in 1:Ncond_2) {
    tempDat <- subset(tempAll,traindiffcond==cond_2[c])
    file_name <- paste0('exp2/train/trainfit',cond_2[c],subs_2[i],'.Rdata')
    if(file.exists(file_name)){
      load(file_name)
    }
    else{ #if not, fit the model
      optimal_params <- DEoptim(chi_square_optim_DDM, # function to optimize
                                # a,ter,z,ntrials,sigma,dt,t2time,vratio,alpha,beta,v
                                lower = c( 0, 0, 0, 5000, .1, .0025, 0,   1,0), 
                                upper = c(.2, 2, 0, 5000, .1, .0025, 0, 1,.5),
                                observations = tempDat,returnFit = 1, binning = F,
                                control=c(itermax=1000,steptol=100,reltol=.001,NP=30))
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

bound <- matrix(NA,Nsub_2,Ncond_2);v <- matrix(NA,Nsub_2,Ncond_2);
ter <- matrix(NA,Nsub_2,Ncond_2);
conf_rt <- matrix(NA,Nsub_2,Ncond_2); resid <- matrix(NA,Nsub_2,Ncond_2)
#Adjust the number of drift parameters to the model loaded
v2 <- matrix(NA,Nsub_2,Ncond_2);v3 <- matrix(NA,Nsub_2,Ncond_2);
for(i in 1:Nsub_2){
  for(c in 1:Ncond_2){
    print(paste('Running participant',i,'from',Nsub_2,"condition",c))
    file_name <- paste0('exp2/test/testfit',cond_2[c],subs_2[i],'.Rdata')
    load(file_name)
    # plot(results$member$bestvalit[1:results$optim$iter],ylab='Goal function',ylim=c(-1,1),frame=F,type='l',main=paste('sub',i,'condition',cond_2[c]))
    bound[i,c] <- results$optim$bestmem[1]
    ter[i,c] <- results$optim$bestmem[2]
    conf_rt[i,c] <- results$optim$bestmem[7]
    v[i,c] <- results$optim$bestmem[9]
    v2[i,c] <-   results$optim$bestmem[10]
    v3[i,c] <-   results$optim$bestmem[11]
    resid[i,c] <- results$optim$bestval
  }
}

param2 <- data.frame(drift = c(v,v2,v3),bound=rep(bound,Ndiff),ter=rep(ter,Ndiff),
                     sub=rep(subs_2,Ndiff*Ncond_2),
                     condition=rep(cond_2,each=Nsub_2,length.out=Nsub_2*Ncond_2*Ndiff),
                     difflevel=rep(coh,each=Nsub_2*Ncond_2),exp=2,resid=rep(resid,Ndiff))

##' Aggregate train and test
bounds2 <- data.frame(bound = c(bound_train,bound),
                      phase = rep(c("train","main"),each=length(bound)), 
                      sub = rep(subs_2,Ncond*2),
                      condition = rep(cond_1,each = Nsub_2,length.out=Nsub_2*Ncond*2))
ters2 <- data.frame(ter = c(ter_train,ter),
                    phase = rep(c("train","main"),each=length(ter)), 
                    sub = rep(subs_2,Ncond*2),
                    condition = rep(cond_1,each = Nsub_2,length.out=Nsub_2*Ncond*2))


## Fit subjective drift ====
go_to("results")
Ndiff <- 1 #Only one difficulty level in the training phase
if (!(file.exists("cost_vs_exp2_cor.csv"))) {
  means <- matrix(NA,nrow=Ncond_2*Nsub_2,ncol=length(drifts)) 
  stds <- matrix(NA,nrow=Ncond_2*Nsub_2,ncol=length(drifts))
  s <- 1; cond <- 1
  while (s <= Nsub_2) {
    while (cond <= Ncond_2) {
      cost_conf_cor <- matrix(NA,nrow=nrepeat,ncol=length(drifts))
      cost_conf_fb <- matrix(NA,nrow=nrepeat,ncol=length(drifts))
      print(paste("Running participant",s,"of",Nsub_2,"condition",cond))
      tempDat <- subset(Data2_train,sub==subs_2[s]&traindiffcond==cond_2[cond])
      tempDat_test <- subset(Data2,sub==subs_2[s]&traindiffcond==cond_2[cond])
      ntrial <- dim(tempDat_test)[1]
      ntrial_train <- dim(tempDat)[1]
      temp_par <- c(bound_train[s,cond],ter_train[s,cond],0,nrepeat,
                    .1,.001,1,v_train[s,cond],v2_train[s,cond],v3_train[s,cond])
      
      # First generate trials from estimated DDM parameters in the training phase
      temp <- chi_square_optim_DDM_fullconfRT(temp_par,observations=tempDat_test,returnFit=0)
      
      # match to the heatmap
      temp$closest_evdnc2 <- match.closest(temp$evidence2,ev_mapping)
      temp$temprt2 <- temp$rt2;
      temp$temprt2[temp$temprt2>5] <- 5 #heatmap doesn't go higher (5 seconds)
      temp$temprt2 <- temp$temprt2*timesteps/5 #scale with the heatmap
      
      # We compute the cost function several times to take noise into account
      temp$Nrep <- rep(1:nrepeat,each=ntrial/Ndiff,length.out=ntrial*nrepeat)        
      
      bar <- txtProgressBar(0,length(drifts),style=3,char="#")
      for (d in 1:length(drifts)) {
        load(paste0("heatmaps/hm_",drifts[d],"_filled.Rdata"))
        hm_low <- output$lower; hm_up <- output$upper
        hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
        
        temp[temp$resp==1,]$cj <- hmvec_up[(temp[temp$resp==1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==1,]$temprt2)]
        temp[temp$resp==-1,]$cj <- hmvec_low[(temp[temp$resp==-1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==-1,]$temprt2)]
        
        
        for (i in 1:nrepeat) {
          pred_sample <- subset(temp,Nrep==i)
          #' Randomly sample 40 trials for each drift to get the same number of 
          #' trials as in the training
          pred_sample <- do.call(rbind,
                                 lapply(split(pred_sample, pred_sample$drift),
                                        function(x) x[sample(nrow(x), ntrial_train/Ndiff),]))
          diff <- sum((tempDat$cor - pred_sample$cj)^2)
          cost_conf_cor[i,d] <- diff
          diff <- sum((tempDat$fb - pred_sample$cj)^2)
          cost_conf_fb[i,d] <- diff
        }
        setTxtProgressBar(bar,d)
      }
      temp_df_cor <- data.frame(cost = as.vector(cost_conf_cor),
                            Nrep = rep(1:nrepeat,length(drifts)),
                            Vs = rep(drifts,each=nrepeat),
                            sub = subs_2[s], traindiffcond = cond_2[cond])
      temp_df_fb <- data.frame(cost = as.vector(cost_conf_fb),
                            Nrep = rep(1:nrepeat,length(drifts)),
                            Vs = rep(drifts,each=nrepeat),
                            sub = subs_2[s], traindiffcond = cond_2[cond])
      if (s==1 & cond==1) {
        cost_df_cor <- temp_df_cor
        cost_df_fb <- temp_df_fb
      }else{
        cost_df_cor <- rbind(cost_df_cor,temp_df_cor)
        cost_df_fb <- rbind(cost_df_fb,temp_df_fb)
      }
      means[s+Nsub_2*(cond-1),] <- colMeans(cost_conf)
      stds[s+Nsub_2*(cond-1),] <- colSds(cost_conf)
      cond <- cond + 1
    }
    s <- s + 1
    cond <- 1
  }
  # save(means,file="means_exp2_full.Rdata")
  # save(stds,file="stds_exp2_full.Rdata")
  write.csv(cost_df_cor,file="cost_vs_exp2_cor.csv")
  write.csv(cost_df_fb,file="cost_vs_exp2_fb.csv")
}else{
  load("means_exp2.Rdata")
  load("stds_exp2.Rdata")
  cost_df_cor <- read.csv("cost_vs_exp2_cor.csv")
  cost_df_fb <- read.csv("cost_vs_exp2_fb.csv")
}

means_fullconfRT <- with(cost_df_cor,aggregate(cost,by=list(Vs=Vs,sub=sub,traindiffcond=traindiffcond),mean))
means_fullconfRT <- cast(means_fullconfRT,traindiffcond+sub~Vs)

vs_smooth2 <- sapply(seq(nrow(means_fullconfRT)), function(i) {
  j <- as.numeric(means_fullconfRT[i,3:dim(means_fullconfRT)[2]])
  j <- lowess(j,f=.05)
  j <- which.min(j$y)
  c(drifts[j])
})

result <- sapply(seq(nrow(means)),function(i) {
  j <- which.min(means[i,])
  c(j)
})
Vs2 <- drifts[result]
Vs2_matrix <- matrix(Vs2,nrow=Nsub_2,ncol=Ncond_2)

df2 <- data.frame(Vs=vs_smooth2,bound = c(bound_train), ter = c(ter_train),
                  Vo=c(v_train),sub=rep(subs_2,Ncond_2),condition=rep(cond_2,each=Nsub_2),
                  Vs_raw=Vs2)

## Generate model prediction ====
go_to("results")
if (file.exists("model_prediction_exp2.csv")) {
  Simuls2 <- read.table("model_prediction_exp2.csv")
}else{
  rm(Simuls2)
  for(i in 1:Nsub_2){
    print(paste('simulating',i,'from',Nsub_2))
    for(c in 1:Ncond_2){
      tempdat <- subset(Data2, sub==subs_2[i] & traindiffcond==cond_2[c])
      load(paste0("heatmaps/hm_",Vs2_matrix[i,c],"_filled.Rdata"))
      hm_low <- output$lower; hm_up <- output$upper
      hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
      temp <- chi_square_optim_DDM(c(bound[i,c],ter[i,c],0,nsim,.1,.001,conf_rt[i,c],1,v[i,c],v2[i,c],v3[i,c]),NULL,0)
      
      #match to the heatmap
      temp$closest_evdnc2 <- match.closest(temp$evidence2,ev_mapping)
      temp$temprt2 <- temp$rt2;
      temp$temprt2[temp$temprt2>5] <- 5 #heatmap doesn't go higher
      temp$temprt2 <- temp$temprt2*timesteps/5 #scale with the heatmap, between 0 and 2000
      
      temp[temp$resp==1,]$cj <- hmvec_up[(temp[temp$resp==1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==1,]$temprt2)]
      temp[temp$resp==-1,]$cj <- hmvec_low[(temp[temp$resp==-1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==-1,]$temprt2)]
      
      temp$cj_cont <- temp$cj
      
      
      if(!exists('Simuls2')){ Simuls2 <- cbind(temp,rep(cond_2[c],nsim),rep(subs_2[i],nsim))
      }else{ Simuls2 <- rbind(Simuls2,cbind(temp,rep(cond_2[c],nsim),rep(subs_2[i],nsim)))
      }
    }
  }
  Simuls2 <- data.frame(Simuls2);names(Simuls2) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','closest_evdnc2',"temprt2",'cj_cont','condition','sub')
  
  coherences <- sort(unique(Data2$coh))
  Simuls2$coh <- 0
  for (i in 1:Nsub_2) {
    for(d in 1:length(coherences)) Simuls2$coh[Simuls2$sub==subs_2[i] & Simuls2$drift %in% c(unique(subset(Simuls2,sub==subs_2[i])$drift)[d],unique(subset(Simuls2,sub==subs_2[i])$drift)[d+3],unique(subset(Simuls2,sub==subs_2[i])$drift)[d+6])] <- coherences[d] #recode drift to coherence
  }
  write.csv(Simuls2,file = "model_prediction_exp2.csv")
}
# Vs fitting procedures comparison ----------------------------------------
#' Originally, the median confidence RT was considered for fitting Vs, with
#' the final value of Vs determined by the mean of 10 repetitions of the fitting
#' We proceeded to go instead with the full confidence RT distribution to preserve
#' more information.
#' We also explore different approaches to estimate Vs

means_fullconfRT <- with(cost_df,aggregate(cost,by=list(Vs=Vs,sub=sub,selfconf=selfconf),mean))
means_fullconfRT <- cast(means_fullconfRT,selfconf+sub~Vs)
result_mean <- means_fullconfRT[,3:dim(means_fullconfRT)[2]]
result_mean <- as.matrix(result_mean)
result_mean <- sapply(seq(nrow(result_mean)),function(i) {
  j <- which.min(result_mean[i,])
  c(j)
})
Vs_mean <- drifts[result_mean]

medians <- with(cost_df,aggregate(cost,by=list(Vs=Vs,sub=sub,selfconf=selfconf),median))
medians <- cast(medians,selfconf+sub~Vs)
result_median <- medians[,3:dim(medians)[2]]
result_median <- as.matrix(result_median)
result_median <- sapply(seq(nrow(result_median)),function(i) {
  j <- which.min(result_median[i,])
  c(j)
})
Vs_median <- drifts[result_median]

N_vs <- 4
confRT <- c("median",rep("dist",N_vs-1))
estimate_func <- c("mean","mean","median","mean")
type <- c("mean","mean_fullconfRT","median","smooth")
Vs_compare1 <- data.frame(Vs=c(Vs1,Vs_mean,Vs_median,vs_smooth1),sub=rep(subs1,Ncond*N_vs),
                         condition=rep(cond_1,each=N1,length.out=Ncond*N1*N_vs),
                         type=rep(type,each=Ncond*N1),estimate_func=rep(estimate_func,each=Ncond*N1),
                         confRT=rep(confRT,each=Ncond*N1))

means_fullconfRT <- with(cost_df_cor,aggregate(cost,by=list(Vs=Vs,sub=sub,traindiffcond=traindiffcond),mean))
means_fullconfRT <- cast(means_fullconfRT,traindiffcond+sub~Vs)
result_mean <- means_fullconfRT[,3:dim(means_fullconfRT)[2]]
result_mean <- as.matrix(result_mean)
result_mean <- sapply(seq(nrow(result_mean)),function(i) {
  j <- which.min(result_mean[i,])
  c(j)
})
Vs_mean <- drifts[result_mean]

medians <- with(cost_df_cor,aggregate(cost,by=list(Vs=Vs,sub=sub,traindiffcond=traindiffcond),median))
medians <- cast(medians,traindiffcond+sub~Vs)
result_median <- medians[,3:dim(medians)[2]]
result_median <- as.matrix(result_median)
result_median <- sapply(seq(nrow(result_median)),function(i) {
  j <- which.min(result_median[i,])
  c(j)
})
Vs_median <- drifts[result_median]

means_fullconfRT <- with(cost_df_fb,aggregate(cost,by=list(Vs=Vs,sub=sub,traindiffcond=traindiffcond),mean))
means_fullconfRT <- cast(means_fullconfRT,traindiffcond+sub~Vs)
result_mean <- means_fullconfRT[,3:dim(means_fullconfRT)[2]]
result_mean <- as.matrix(result_mean)
result_mean <- sapply(seq(nrow(result_mean)),function(i) {
  j <- which.min(result_mean[i,])
  c(j)
})
Vs_mean_fb <- drifts[result_mean]

medians <- with(cost_df_fb,aggregate(cost,by=list(Vs=Vs,sub=sub,traindiffcond=traindiffcond),median))
medians <- cast(medians,traindiffcond+sub~Vs)
result_median <- medians[,3:dim(medians)[2]]
result_median <- as.matrix(result_median)
result_median <- sapply(seq(nrow(result_median)),function(i) {
  j <- which.min(result_median[i,])
  c(j)
})
Vs_median_fb <- drifts[result_median]

N_vs <- 6
type <- c("mean","mean_fullconfRT_fb","median_fb","mean_fullconfRT","median","smooth")
confRT <- c("median",rep("dist",N_vs-1))
estimate_func <- c("mean","mean","median","mean","median","mean")
feedback <- c(rep("block",3),rep("trial",3))
Vs_compare2 <- data.frame(Vs=c(Vs2,Vs_mean_fb,Vs_median_fb,Vs_mean,Vs_median,vs_smooth2),
                          sub=rep(subs_2,Ncond*N_vs),
                         condition=rep(cond_2,each=Nsub_2,length.out=Ncond*Nsub_2*N_vs),
                         type=rep(type,each=Ncond*Nsub_2), confRT=rep(confRT,each=Ncond*Nsub_2),
                         feedback=rep(feedback,each=Ncond*Nsub_2), estimate_func=rep(estimate_func,each=Ncond*Nsub_2))

#' Questions : 
#' - Is there a difference in fitted Vs between median confRT and full distribution ?
#' - Which of the mean/median/min provide the most accurate estimate for Vs ?
#' - Exp2 : Does the trial-by-trial feedback give different results than the blockwise feedback ? 

## Exp 1: fake feedback
# Median confRT vs full distribution
m <- lmer(Vs~confRT*condition + (1|sub),data=subset(Vs_compare1,estimate_func=="mean"))
anova(m)
# Post-hoc test within each condition
emm <- emmeans(m, ~ confRT|condition) 
pairs(emm) # Slightly higher Vs with the full distribution in the positive FB condition 
with(subset(Vs_compare1,estimate_func=="mean"),aggregate(Vs,by=list(confRT,condition),mean)) # Show mean estimates

# Mean vs Median vs Min
m <- lmer(Vs~estimate_func*condition + (1|sub),data=subset(Vs_compare1,confRT=="dist"))
anova(m) # No difference between mean and median

# Smooth vs no Smooth
m <- lmer(Vs~type*condition + (1|sub),data=subset(Vs_compare1,type %in% c("mean_fullconfRT","smooth")))
anova(m) 

## Exp 2: Training difficulty
# Median confRT vs full distribution
m <- lmer(data = subset(Vs_compare2,feedback=="block"&estimate_func=="mean"),
          Vs~confRT*condition + (1|sub))
anova(m)
# Post-hoc test within each condition
emm <- emmeans(m, ~ confRT|condition)
pairs(emm) # Median confRT has lower Vs estimates in the easy condition

# Trial-by-trial FB vs block FB + aggregation function of the repetitions
m <- lmer(Vs~estimate_func*condition*feedback + (1|sub),data=subset(Vs_compare2,confRT=="dist"))
anova(m)
m <- lmer(Vs~condition*feedback + (1|sub),data=subset(Vs_compare2,confRT=="dist"&estimate_func=="mean"))
anova(m)
# Post-hoc test within each condition
emm <- emmeans(m, ~ estimate_func|condition) 
pairs(emm) # Higher Vs estimate using the mean in the easy condition 
with(Vs_compare2,aggregate(Vs,by=list(estimate_func,condition),mean)) # Show mean estimates

# Smooth vs no Smooth
m <- lmer(Vs~type*condition + (1|sub),data=subset(Vs_compare2,type %in% c("mean_fullconfRT","smooth")))
anova(m) 

par(mfrow=c(1,3))
for (i in 1:Nsub_2) {
  for (c in 1:Ncond_2) {
    tempmean <- as.numeric(means_fullconfRT[Nsub_2*(c-1)+i,])
    tempmean <- tempmean[complete.cases(tempmean)]
    smoothed <- lowess(tempmean[2:501],f=.05)
    plot(tempmean[2:501],main=paste(subs_2[i],cond_2[c],"mean"),
         xlab="Vs",ylab="Mean over repetitions",xaxt='n')
    lines(smoothed,col="green",lwd=2)
    axis(1,at=seq(0,500,100),labels = seq(0,.5,.1))
    abline(v=which.min(tempmean[2:501]),col="red")
    abline(v=which.min(smoothed$y),col="green")
  }
}

# Stat tests ------------------------------------------------------------
# Exp1 ====
#DDM train
df$sub <- as.factor(df$sub)
m <- lmer(Vs ~ condition + (1|sub),data=df); anova(m);
m <- aov(Vs ~ condition+ Error(sub/condition), data = df); summary(m) #Equivalent
m <- lmer(Vo ~condition + (1|sub), data = df); anova(m)
m <- lmer(bound ~ condition + (1|sub),data=df); anova(m);
m <- lmer(ter ~ condition + (1|sub),data=df); anova(m);

#DDM test
param_1$sub <- as.factor(param_1$sub)
test_bound1 <- with(param_1,aggregate(bound,by=list(condition=condition,sub=sub),mean))
test_ter1 <- with(param_1,aggregate(ter,by=list(condition=condition,sub=sub),mean))
m <- lmer(x ~ condition + (1|sub),data=test_bound1); anova(m);
m <- lmer(x ~ condition + (1|sub),data=test_ter1); anova(m);
m <- lmer(drift ~ condition*difflevel + (1|sub),data=param_1); anova(m);

#Train vs Test
bounds1$sub <- as.factor(bounds1$sub)
vs$sub <- as.factor(vs$sub)
ters1$sub <- as.factor(ters1$sub)
m <- lmer(bound ~ phase*condition + (condition|sub),data = bounds1); anova(m)
m <- lmer(v ~ phase*condition*difficulty + (difficulty|sub),data = vs); anova(m)
m <- lmer(ter ~ phase*condition + (condition|sub),data = ters1); anova(m)
## Exp2 ====
#DDM train
df2$sub <- as.factor(df2$sub)
m <- aov(Vs ~ condition+ Error(sub/condition), data = df2); summary(m)
m <- lmer(Vs ~condition + (1|sub), data = df2); anova(m)
m <- lmer(Vo ~ condition + (1|sub),data=df2); anova(m);
m <- lmer(bound ~ condition + (1|sub),data=df2); anova(m);
m <- lmer(ter ~ condition + (1|sub),data=df2); anova(m)

#DDM test
param2$sub <- as.factor(param2$sub)
test_bound2 <- with(param2,aggregate(bound,by=list(condition=condition,sub=sub),mean))
test_ter2 <- with(param2,aggregate(ter,by=list(condition=condition,sub=sub),mean))
m <- lmer(x ~ condition + (1|sub),data=test_bound2); anova(m);
m <- lmer(x ~ condition + (1|sub),data=test_ter2); anova(m);
m <- lmer(drift ~ condition*difflevel + (condition|sub),data=param2); anova(m);

#Train vs Test
bounds2$sub <- as.factor(bounds2$sub)
ters2$sub <- as.factor(ters2$sub)
m <- lmer(bound ~ phase*condition + (1|sub),data = bounds2); anova(m)
m <- lmer(ter ~ phase*condition + (condition|sub),data = ters2); anova(m)

sim_cj1 <- with(Simuls,aggregate(cj,by=list(condition=condition,coh=coh,sub=sub),mean))
sim_cj2 <- with(Simuls2,aggregate(cj,by=list(condition=condition,coh=coh,sub=sub),mean))
sim_cj1$sub <- as.factor(sim_cj1$sub)
sim_cj2$sub <- as.factor(sim_cj2$sub)
m <- lmer(x ~ condition*coh + (condition|sub), data = sim_cj1); anova(m)
m <- lmer(x ~ condition*coh + (condition|sub), data = sim_cj2); anova(m)
m <- aov(x ~ condition*coh+ Error(sub/condition), data = sim_cj1); summary(m)
m <- aov(x ~ condition*coh+ Error(sub/condition), data = sim_cj2); summary(m)
# Plot Layout -------------------------------------------------------------
go_to("plot")

jpeg(
  filename="results.jpeg",
  width=13,
  height=18,
  units="in",
  res=500)
# layout(matrix(c(1,3,7,9,11,12,1,5,7,9,11,12,2,4,8,10,11,13,2,6,8,10,11,13),ncol=4),heights = c(.4,1,1.5,1.5,.2,1.5))
layout(matrix(c(1,3,7,9,10,1,5,7,9,10,2,4,8,9,11,2,6,8,9,11),ncol=4),heights = c(.4,1,2,.2,2))

#' Add legend for both experiments on top
par(mar=c(0,0,0,0))
plot.new()
legend("top",legend=c("Negative","Average","Positive"),
       title = "Experiment 1: Feedback condition",pch=rep(16,3),bty = "n",inset=0,
       cex = cex_legend,col=c("brown3","cyan4","darkgoldenrod3"), horiz = T)

plot.new()
legend("top",legend=c("Difficult","Medium","Easy"),
       title = "Experiment 2: Training condition",pch=rep(16,3),bty = "n",inset=0, 
       cex = cex_legend,col=c("brown3","cyan4","darkgoldenrod3"), horiz = T)
par(mar=c(5,5,2,2)+0.1)
# Plot overlay accuracy ---------------------------------------------------

## Experiment 1
#Aggregate conf for data
corlow <- with(subset(Data1,selfconf=="lowSC"),aggregate(cor,by=list(sub,coh),mean));names(corlow) <- c('sub','coh','cor')
corlow <- cast(corlow,sub~coh,)
cormed <- with(subset(Data1,selfconf=="mediumSC"),aggregate(cor,by=list(sub,coh),mean));names(cormed) <- c('sub','coh','cor')
cormed <- cast(cormed,sub~coh)
corhigh <- with(subset(Data1,selfconf=="highSC"),aggregate(cor,by=list(sub,coh),mean));names(corhigh) <- c('sub','coh','cor')
corhigh <- cast(corhigh,sub~coh)


#aggregate cor for model
x <- corlow[,c(2:4)];xmed <- cormed[,c(2:4)];xhigh <- corhigh[,c(2:4)]
n <- length(x)

x <- x[,c("hard","average","easy")];
xmed <- xmed[,c("hard","average","easy")];
xhigh <- xhigh[,c("hard","average","easy")]

#Immediate condition
simDat_low <- subset(Simuls,condition=="lowSC");simDat_low <- simDat_low[c('sub','coh','cor')]
simDat_med <- subset(Simuls,condition=="mediumSC");simDat_med <- simDat_med[c('sub','coh','cor')]
simDat_high <- subset(Simuls,condition=="highSC");simDat_high <- simDat_high[c('sub','coh','cor')]

#aggregate cor for model
snrcorlowSim <- aggregate(.~sub+coh,simDat_low,mean)
snrcormedSim <- aggregate(.~sub+coh,simDat_med,mean)
snrcorhighSim <- aggregate(.~sub+coh,simDat_high,mean)
x_sim = cast(snrcorlowSim,sub~coh,value='cor');xmed_sim = cast(snrcormedSim,sub~coh,value='cor');xhigh_sim = cast(snrcorhighSim,sub~coh,value='cor')
x <- x[,c("hard","average","easy")];xmed <- xmed[,c("hard","average","easy")];x_sim <- x_sim[,c("hard","average","easy")];xmed_sim <- xmed_sim[,c("hard","average","easy")];
xhigh_sim <- xhigh_sim[,c("hard","average","easy")];xhigh <- xhigh[,c("hard","average","easy")]


stripchart(x, ylim=c(0.6,1), xlim=c(-.05,n-1), vertical = TRUE, col="white",
           frame=F,xaxt='n',yaxt='n',xlab="Trial difficulty",ylab="Accuracy",
           cex.lab=cex_lab/2)
mtext("A.", at = -.55, line = 1, cex = cex_title, font = 2)
axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);
axis(2,seq(0.6,1,.1),cex.axis=1.5)
# mtext("Accuracy",2,at=.8,line=2.5,cex=cex_lab);
# mtext("Trial difficulty",1,3,at=1,cex=cex_lab)
for(i in seq(.6,1,.1)) abline(h=i,col="lightgrey",lty = "dashed")
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N1)),
          (colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(1,0,0,.2))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xmed_sim,na.rm=T) + (colSds(as.matrix(xmed_sim))/sqrt(N1)),
          (colMeans(xmed_sim,na.rm=T) - colSds(as.matrix(xmed_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(1,.5,0,.2))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xhigh_sim,na.rm=T) + (colSds(as.matrix(xhigh_sim))/sqrt(N1)),
          (colMeans(xhigh_sim,na.rm=T) - colSds(as.matrix(xhigh_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(0,0,1,.2))
means <- sapply(x, mean);n<- length(x)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="brown3",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N1),lwd=lwdgr,col="brown3")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="cyan4",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N1),lwd=lwdgr,col="cyan4")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="darkgoldenrod3",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N1),lwd=lwdgr,col="darkgoldenrod3")

## Experiment 2
#Aggregate conf for Data2
corlow <- with(subset(Data2,traindiffcond=="hard"),aggregate(cor,by=list(sub,coh),mean));names(corlow) <- c('sub','coh','cor')
corlow <- cast(corlow,sub~coh,)
cormed <- with(subset(Data2,traindiffcond=="average"),aggregate(cor,by=list(sub,coh),mean));names(cormed) <- c('sub','coh','cor')
cormed <- cast(cormed,sub~coh)
corhigh <- with(subset(Data2,traindiffcond=="easy"),aggregate(cor,by=list(sub,coh),mean));names(corhigh) <- c('sub','coh','cor')
corhigh <- cast(corhigh,sub~coh)


#aggregate cor for model
x <- corlow[,c(2:4)];xmed <- cormed[,c(2:4)];xhigh <- corhigh[,c(2:4)]
n <- length(x)

x <- x[,c("hard","average","easy")];
xmed <- xmed[,c("hard","average","easy")];
xhigh <- xhigh[,c("hard","average","easy")]

#Immediate condition
simDat_low <- subset(Simuls2,condition=="hard");simDat_low <- simDat_low[c('sub','coh','cor')]
simDat_med <- subset(Simuls2,condition=="average");simDat_med <- simDat_med[c('sub','coh','cor')]
simDat_high <- subset(Simuls2,condition=="easy");simDat_high <- simDat_high[c('sub','coh','cor')]

#aggregate cor for model
snrcorlowSim <- aggregate(.~sub+coh,simDat_low,mean)
snrcormedSim <- aggregate(.~sub+coh,simDat_med,mean)
snrcorhighSim <- aggregate(.~sub+coh,simDat_high,mean)
x_sim = cast(snrcorlowSim,sub~coh,value='cor');xmed_sim = cast(snrcormedSim,sub~coh,value='cor');xhigh_sim = cast(snrcorhighSim,sub~coh,value='cor')
x <- x[,c("hard","average","easy")];xmed <- xmed[,c("hard","average","easy")];x_sim <- x_sim[,c("hard","average","easy")];xmed_sim <- xmed_sim[,c("hard","average","easy")];
xhigh_sim <- xhigh_sim[,c("hard","average","easy")];xhigh <- xhigh[,c("hard","average","easy")]


stripchart(x, ylim=c(0.6,1), xlim=c(-.05,n-1), vertical = TRUE, col="white",
           frame=F,xaxt='n',main=NULL,yaxt='n',xlab="Trial difficulty",ylab="Accuracy",cex.lab=cex_lab/2)
mtext("B.", at = -.55, line = 1, cex = cex_title, font = 2)
axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);
axis(2,seq(0.6,1,.1),cex.axis=1.5)
# mtext("Accuracy",2,at=.8,line=2.5,cex=cex_lab);
# mtext("Trial difficulty",1,3,at=1,cex=cex_lab)
for(i in seq(.6,1,.1)) abline(h=i,col="lightgrey",lty = "dashed")
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N1)),
          (colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(1,0,0,.2))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xmed_sim,na.rm=T) + (colSds(as.matrix(xmed_sim))/sqrt(N1)),
          (colMeans(xmed_sim,na.rm=T) - colSds(as.matrix(xmed_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(1,.5,0,.2))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xhigh_sim,na.rm=T) + (colSds(as.matrix(xhigh_sim))/sqrt(N1)),
          (colMeans(xhigh_sim,na.rm=T) - colSds(as.matrix(xhigh_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(0,0,1,.2))
means <- sapply(x, mean);n<- length(x)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="brown3",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N1),lwd=lwdgr,col="brown3")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="cyan4",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N1),lwd=lwdgr,col="cyan4")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="darkgoldenrod3",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N1),lwd=lwdgr,col="darkgoldenrod3")

# Plot overlay RT ---------------------------------------------------------
## Experiment 1
#Aggregate conf for data
rtlow <- with(subset(Data1,selfconf=="lowSC"),aggregate(rt,by=list(sub,coh),mean));
names(rtlow) <- c('sub','coh','rt')
rtlow <- cast(rtlow,sub~coh,)
rtmed <- with(subset(Data1,selfconf=="mediumSC"),aggregate(rt,by=list(sub,coh),mean));
names(rtmed) <- c('sub','coh','rt')
rtmed <- cast(rtmed,sub~coh)
rthigh <- with(subset(Data1,selfconf=="highSC"),aggregate(rt,by=list(sub,coh),mean));
names(rthigh) <- c('sub','coh','rt')
rthigh <- cast(rthigh,sub~coh)


#aggregate rt for model
x <- rtlow[,c(2:4)];xmed <- rtmed[,c(2:4)];xhigh <- rthigh[,c(2:4)]
n <- length(x)

x <- x[,c("hard","average","easy")];
xmed <- xmed[,c("hard","average","easy")];
xhigh <- xhigh[,c("hard","average","easy")]

#Immediate condition
simDat_low <- subset(Simuls,condition=="lowSC");
simDat_low <- simDat_low[c('sub','coh','rt')]
simDat_med <- subset(Simuls,condition=="mediumSC");
simDat_med <- simDat_med[c('sub','coh','rt')]
simDat_high <- subset(Simuls,condition=="highSC");
simDat_high <- simDat_high[c('sub','coh','rt')]

#aggregate rt for model
snrrtlowSim <- aggregate(.~sub+coh,simDat_low,mean)
snrrtmedSim <- aggregate(.~sub+coh,simDat_med,mean)
snrrthighSim <- aggregate(.~sub+coh,simDat_high,mean)
x_sim = cast(snrrtlowSim,sub~coh,value='rt');xmed_sim = cast(snrrtmedSim,sub~coh,value='rt');xhigh_sim = cast(snrrthighSim,sub~coh,value='rt')
x <- x[,c("hard","average","easy")];xmed <- xmed[,c("hard","average","easy")];x_sim <- x_sim[,c("hard","average","easy")];xmed_sim <- xmed_sim[,c("hard","average","easy")];
xhigh_sim <- xhigh_sim[,c("hard","average","easy")];xhigh <- xhigh[,c("hard","average","easy")]


stripchart(x, ylim=c(0.6,1.1), xlim=c(-.05,n-1), vertical = TRUE, col="white",
           frame=F,xaxt='n',main=NULL,yaxt='n',xlab="Trial difficulty",ylab="RT (s)",cex.lab=cex_lab/2)
axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);
axis(2,seq(0.6,1.1,.1),cex.axis=1.5)
# mtext("RT (s)",2,at=.85,line=2.5,cex=cex_lab);
# mtext("Trial difficulty",1,3,at=1,cex=cex_lab)
for(i in seq(.6,1.1,.1)) abline(h=i,col="lightgrey",lty = "dashed")

polygon(c(0:(n-1),(n-1):0),
        c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N1)),
          (colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(1,0,0,.2))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xmed_sim,na.rm=T) + (colSds(as.matrix(xmed_sim))/sqrt(N1)),
          (colMeans(xmed_sim,na.rm=T) - colSds(as.matrix(xmed_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(1,.5,0,.2))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xhigh_sim,na.rm=T) + (colSds(as.matrix(xhigh_sim))/sqrt(N1)),
          (colMeans(xhigh_sim,na.rm=T) - colSds(as.matrix(xhigh_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(0,0,1,.2))
means <- sapply(x, mean);n<- length(x)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="brown3",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N1),lwd=lwdgr,col="brown3")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="cyan4",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N1),lwd=lwdgr,col="cyan4")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="darkgoldenrod3",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N1),lwd=lwdgr,col="darkgoldenrod3")

## Experiment 2
#Aggregate conf for Data2
rtlow <- with(subset(Data2,traindiffcond=="hard"),aggregate(rt,by=list(sub,coh),mean));names(rtlow) <- c('sub','coh','rt')
rtlow <- cast(rtlow,sub~coh,)
rtmed <- with(subset(Data2,traindiffcond=="average"),aggregate(rt,by=list(sub,coh),mean));names(rtmed) <- c('sub','coh','rt')
rtmed <- cast(rtmed,sub~coh)
rthigh <- with(subset(Data2,traindiffcond=="easy"),aggregate(rt,by=list(sub,coh),mean));names(rthigh) <- c('sub','coh','rt')
rthigh <- cast(rthigh,sub~coh)


#aggregate rt for model
x <- rtlow[,c(2:4)];xmed <- rtmed[,c(2:4)];xhigh <- rthigh[,c(2:4)]
n <- length(x)

x <- x[,c("hard","average","easy")];
xmed <- xmed[,c("hard","average","easy")];
xhigh <- xhigh[,c("hard","average","easy")]

#Immediate condition
simDat_low <- subset(Simuls2,condition=="hard");simDat_low <- simDat_low[c('sub','coh','rt')]
simDat_med <- subset(Simuls2,condition=="average");simDat_med <- simDat_med[c('sub','coh','rt')]
simDat_high <- subset(Simuls2,condition=="easy");simDat_high <- simDat_high[c('sub','coh','rt')]

#aggregate rt for model
snrrtlowSim <- aggregate(.~sub+coh,simDat_low,mean)
snrrtmedSim <- aggregate(.~sub+coh,simDat_med,mean)
snrrthighSim <- aggregate(.~sub+coh,simDat_high,mean)
x_sim = cast(snrrtlowSim,sub~coh,value='rt');xmed_sim = cast(snrrtmedSim,sub~coh,value='rt');xhigh_sim = cast(snrrthighSim,sub~coh,value='rt')
x <- x[,c("hard","average","easy")];xmed <- xmed[,c("hard","average","easy")];x_sim <- x_sim[,c("hard","average","easy")];xmed_sim <- xmed_sim[,c("hard","average","easy")];
xhigh_sim <- xhigh_sim[,c("hard","average","easy")];xhigh <- xhigh[,c("hard","average","easy")]


stripchart(x, ylim=c(0.6,1.1), xlim=c(-.05,n-1), vertical = TRUE, col="white",
           frame=F,xaxt='n',main=NULL,yaxt='n',xlab="Trial difficulty",ylab="RT (s)",cex.lab=cex_lab/2)
axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);
axis(2,seq(0.6,1.1,.1),cex.axis=1.5)
# mtext("RT (s)",2,at=.85,line=2.5,cex=cex_lab);
# mtext("Trial difficulty",1,3,at=1,cex=cex_lab)
for(i in seq(.6,1.1,.1)) abline(h=i,col="lightgrey",lty = "dashed")
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N1)),
          (colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(1,0,0,.2))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xmed_sim,na.rm=T) + (colSds(as.matrix(xmed_sim))/sqrt(N1)),
          (colMeans(xmed_sim,na.rm=T) - colSds(as.matrix(xmed_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(1,.5,0,.2))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xhigh_sim,na.rm=T) + (colSds(as.matrix(xhigh_sim))/sqrt(N1)),
          (colMeans(xhigh_sim,na.rm=T) - colSds(as.matrix(xhigh_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(0,0,1,.2))
means <- sapply(x, mean);n<- length(x)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="brown3",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N1),lwd=lwdgr,col="brown3")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="cyan4",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N1),lwd=lwdgr,col="cyan4")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="darkgoldenrod3",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N1),lwd=lwdgr,col="darkgoldenrod3")

# Plot Confidence - Empirical Data ----------------------------------------
##' Experiment 1
CJ_SC_diff_data <- with(Data1,aggregate(cj,by=list(sub=sub,selfconf=selfconf, coh=coh),mean));
CJ_SC_diff_data <- cast(CJ_SC_diff_data,sub~selfconf+coh)
average_CJ_SC_diff_data <- with(Data1,aggregate(cj,by=list(selfconf=selfconf,coh=coh),mean));
average_CJ_SC_diff_data <- cast(average_CJ_SC_diff_data,selfconf~coh)

# use family to adjust the font and cex. to adjust font size
CJ_SC_diff_plot = plot(as.numeric(average_CJ_SC_diff_data[1,]),type='n',frame=F,
                       main=NULL,
                       ylab="Confidence",
                       xlab="Trial difficulty",
                       xaxt='n',
                       xlim=c(1,3.3),ylim=c(3,6),
                       cex.axis = cex_lab-1, 
                       cex.lab = cex_lab,
                       family="A")
mtext("C.", at=.7, line = 1, cex = cex_title, font = 2)
axis(1,at=1.1:3.1,labels=c("hard","average","easy"),cex.axis=cex_lab-1,family="A")
abline(h = seq(3,6,0.5), col = "lightgrey", lty = "dashed")

# High SC
for(i in 1:N1) points(jitter(1:3,0.1),CJ_SC_diff_data[i,c(4,2,3)],lty=i,type='p',pch=21,col='white',bg=transp('darkgoldenrod2'))
lines(1:3,average_CJ_SC_diff_data[1,c(4,2,3)],lty=2,type='b',pch=21,
      col='darkgoldenrod3',bg='darkgoldenrod2',lwd=lwddat,cex=cexkl)
# Medium SC
for(i in 1:N1) points(jitter(1.1:3.1,0.1),CJ_SC_diff_data[i,c(10,8,9)],lty=i,type='p',pch=24,col='white',bg=transp('cyan3'))
lines(1.1:3.1,average_CJ_SC_diff_data[3,c(4,2,3)],lty=2,type='b',pch=24,
      col='cyan4',bg='cyan3',lwd=lwddat,cex=cexkl)
# Low SC
for(i in 1:N1) points(jitter(1.2:3.2,0.1),CJ_SC_diff_data[i,c(7,5,6)],lty=i,type='p',pch=22,col='white',bg=transp('brown2'))
lines(1.2:3.2,average_CJ_SC_diff_data[2,c(4,2,3)],lty=2,type='b',pch=22,
      col='brown3',bg="brown2",lwd=lwddat,cex=cexkl)

# plot error bars
error.bar(1:3,colMeans(CJ_SC_diff_data[,c(4,2,3)]),colSds(as.matrix(CJ_SC_diff_data[,c(4,2,3)])/sqrt(N1)),
          length=0,lwd=lwdgr, col='darkgoldenrod3')
error.bar(1.1:3.1,colMeans(CJ_SC_diff_data[,c(10,8,9)]),colSds(as.matrix(CJ_SC_diff_data[,c(10,8,9)])/sqrt(N1)),
          length=0,lwd=lwdgr, col='cyan4')
error.bar(1.2:3.2,colMeans(CJ_SC_diff_data[,c(7,5,6)]),colSds(as.matrix(CJ_SC_diff_data[,c(7,5,6)])/sqrt(N1)),
          length=0,lwd=lwdgr, col='brown3')


## Experiment 2
CJ_SC_diff_data <- with(Data2,aggregate(cj,by=list(sub=sub,traindiffcond=traindiffcond, coh=coh),mean));
CJ_SC_diff_data <- cast(CJ_SC_diff_data,sub~traindiffcond+coh)
average_CJ_SC_diff_data <- with(Data2,aggregate(cj,by=list(traindiffcond=traindiffcond,coh=coh),mean));
average_CJ_SC_diff_data <- cast(average_CJ_SC_diff_data,traindiffcond~coh)

# use family to adjust the font and cex. to adjust font size
CJ_SC_diff_plot = plot(as.numeric(average_CJ_SC_diff_data[1,]),type='n',frame=F,
                       main=NULL,
                       ylab="Confidence",
                       xlab="Trial difficulty",
                       xaxt='n',
                       xlim=c(1,3.3),ylim=c(3,6),
                       cex.axis = cex_lab-1, 
                       cex.lab = cex_lab,
                       family="A")
mtext("D.", at = .7, line = 1, cex = cex_title, font = 2)
axis(1,at=1.1:3.1,labels=c("hard","average","easy"),cex.axis=cex_lab-1,family="A")
abline(h = seq(3,6,0.5), col = "lightgrey", lty = "dashed")

# High SC
for(i in 1:Nsub_2) points(jitter(1:3,0.1),CJ_SC_diff_data[i,c(7,5,6)],lty=i,type='p',pch=21,col='white',bg=transp('darkgoldenrod2'))
lines(1:3,average_CJ_SC_diff_data[2,c(4,2,3)],lty=2,type='b',pch=21,
      col='darkgoldenrod3',bg='darkgoldenrod2',lwd=lwddat,cex=cexkl)
# Medium SC
for(i in 1:Nsub_2) points(jitter(1.1:3.1,0.1),CJ_SC_diff_data[i,c(4,2,3)],lty=i,type='p',pch=24,col='white',bg=transp('cyan3'))
lines(1.1:3.1,average_CJ_SC_diff_data[1,c(4,2,3)],lty=2,type='b',pch=24,
      col='cyan4',bg='cyan3',lwd=lwddat,cex=cexkl)
# Low SC
for(i in 1:Nsub_2) points(jitter(1.2:3.2,0.1),CJ_SC_diff_data[i,c(10,8,9)],lty=i,type='p',pch=22,col='white',bg=transp('brown2'))
lines(1.2:3.2,average_CJ_SC_diff_data[3,c(4,2,3)],lty=2,type='b',pch=22,
      col='brown3',bg="brown2",lwd=lwddat,cex=cexkl)

# plot error bars
error.bar(1:3,colMeans(CJ_SC_diff_data[,c(7,5,6)]),colSds(as.matrix(CJ_SC_diff_data[,c(7,5,6)])/sqrt(Nsub_2)),
          length=0,lwd=lwdgr, col='darkgoldenrod3')
error.bar(1.1:3.1,colMeans(CJ_SC_diff_data[,c(4,2,3)]),colSds(as.matrix(CJ_SC_diff_data[,c(4,2,3)])/sqrt(Nsub_2)),
          length=0,lwd=lwdgr, col='cyan4')
error.bar(1.2:3.2,colMeans(CJ_SC_diff_data[,c(10,8,9)]),colSds(as.matrix(CJ_SC_diff_data[,c(10,8,9)])/sqrt(Nsub_2)),
          length=0,lwd=lwdgr, col='brown3')


# # Plot Confidence per block ----------------------------------------
# ##' Experiment 1
# CJ_SC_diff_data <- with(Data1,aggregate(cj,by=list(sub=sub,selfconf=selfconf, block=block),mean));
# CJ_SC_diff_data <- cast(CJ_SC_diff_data,sub~selfconf+block)
# average_CJ_SC_diff_data <- with(Data1,aggregate(cj,by=list(selfconf=selfconf,block=block),mean));
# average_CJ_SC_diff_data <- cast(average_CJ_SC_diff_data,selfconf~block)
# 
# # use family to adjust the font and cex. to adjust font size
# CJ_SC_diff_plot = plot(as.numeric(average_CJ_SC_diff_data[1,]),type='n',frame=F,
#                        main=NULL,
#                        ylab="Confidence",
#                        xlab="Block",
#                        xaxt='n',
#                        xlim=c(1,3.3),ylim=c(3,6),
#                        cex.axis = cex_lab-1, 
#                        cex.lab = cex_lab,
#                        family="A")
# mtext("E.", at=.7, line = 1, cex = cex_title, font = 2)
# axis(1,at=1.1:3.1,labels=c("block 1","block 2","block 3"),cex.axis=cex_lab-1,family="A")
# abline(h = seq(3,6,0.5), col = "lightgrey", lty = "dashed")
# 
# # High SC
# for(i in 1:N1) points(jitter(1:3,0.1),CJ_SC_diff_data[i,c(2,3,4)],lty=i,type='p',pch=21,col='white',bg=transp('darkgoldenrod2'))
# lines(1:3,average_CJ_SC_diff_data[1,c(2,3,4)],lty=2,type='b',pch=21,
#       col='darkgoldenrod3',bg='darkgoldenrod2',lwd=lwddat,cex=cexkl)
# # Medium SC
# for(i in 1:N1) points(jitter(1.1:3.1,0.1),CJ_SC_diff_data[i,c(8,9,10)],lty=i,type='p',pch=24,col='white',bg=transp('cyan3'))
# lines(1.1:3.1,average_CJ_SC_diff_data[3,c(2,3,4)],lty=2,type='b',pch=24,
#       col='cyan4',bg='cyan3',lwd=lwddat,cex=cexkl)
# # Low SC
# for(i in 1:N1) points(jitter(1.2:3.2,0.1),CJ_SC_diff_data[i,c(5,6,7)],lty=i,type='p',pch=22,col='white',bg=transp('brown2'))
# lines(1.2:3.2,average_CJ_SC_diff_data[2,c(2,3,4)],lty=2,type='b',pch=22,
#       col='brown3',bg="brown2",lwd=lwddat,cex=cexkl)
# 
# # plot error bars
# error.bar(1:3,colMeans(CJ_SC_diff_data[,c(2,3,4)]),colSds(as.matrix(CJ_SC_diff_data[,c(2,3,4)])/sqrt(N1)),
#           length=0,lwd=lwdgr, col='darkgoldenrod3')
# error.bar(1.1:3.1,colMeans(CJ_SC_diff_data[,c(8,9,10)]),colSds(as.matrix(CJ_SC_diff_data[,c(8,9,10)])/sqrt(N1)),
#           length=0,lwd=lwdgr, col='cyan4')
# error.bar(1.2:3.2,colMeans(CJ_SC_diff_data[,c(5,6,7)]),colSds(as.matrix(CJ_SC_diff_data[,c(5,6,7)])/sqrt(N1)),
#           length=0,lwd=lwdgr, col='brown3')
# 
# 
# ## Experiment 2
# CJ_SC_diff_data <- with(Data2,aggregate(cj,by=list(sub=sub,traindiffcond=traindiffcond, block=block),mean));
# CJ_SC_diff_data <- cast(CJ_SC_diff_data,sub~traindiffcond+block)
# average_CJ_SC_diff_data <- with(Data2,aggregate(cj,by=list(traindiffcond=traindiffcond,block=block),mean));
# average_CJ_SC_diff_data <- cast(average_CJ_SC_diff_data,traindiffcond~block)
# 
# # use family to adjust the font and cex. to adjust font size
# CJ_SC_diff_plot = plot(as.numeric(average_CJ_SC_diff_data[1,]),type='n',frame=F,
#                        main=NULL,
#                        ylab="Confidence",
#                        xlab="Block",
#                        xaxt='n',
#                        xlim=c(1,3.3),ylim=c(3,6),
#                        cex.axis = cex_lab-1, 
#                        cex.lab = cex_lab,
#                        family="A")
# mtext("F.", at = .7, line = 1, cex = cex_title, font = 2)
# axis(1,at=1.1:3.1,labels=c("block 1","block 2","block 3"),cex.axis=cex_lab-1,family="A")
# abline(h = seq(3,6,0.5), col = "lightgrey", lty = "dashed")
# 
# # High SC
# for(i in 1:Nsub_2) points(jitter(1:3,0.1),CJ_SC_diff_data[i,c(5,6,7)],lty=i,type='p',pch=21,col='white',bg=transp('darkgoldenrod2'))
# lines(1:3,average_CJ_SC_diff_data[2,c(2,3,4)],lty=2,type='b',pch=21,
#       col='darkgoldenrod3',bg='darkgoldenrod2',lwd=lwddat,cex=cexkl)
# # Medium SC
# for(i in 1:Nsub_2) points(jitter(1.1:3.1,0.1),CJ_SC_diff_data[i,c(2,3,4)],lty=i,type='p',pch=24,col='white',bg=transp('cyan3'))
# lines(1.1:3.1,average_CJ_SC_diff_data[1,c(2,3,4)],lty=2,type='b',pch=24,
#       col='cyan4',bg='cyan3',lwd=lwddat,cex=cexkl)
# # Low SC
# for(i in 1:Nsub_2) points(jitter(1.2:3.2,0.1),CJ_SC_diff_data[i,c(8,9,10)],lty=i,type='p',pch=22,col='white',bg=transp('brown2'))
# lines(1.2:3.2,average_CJ_SC_diff_data[3,c(2,3,4)],lty=2,type='b',pch=22,
#       col='brown3',bg="brown2",lwd=lwddat,cex=cexkl)
# 
# # plot error bars
# error.bar(1:3,colMeans(CJ_SC_diff_data[,c(5,6,7)]),colSds(as.matrix(CJ_SC_diff_data[,c(5,6,7)])/sqrt(Nsub_2)),
#           length=0,lwd=lwdgr, col='darkgoldenrod3')
# error.bar(1.1:3.1,colMeans(CJ_SC_diff_data[,c(2,3,4)]),colSds(as.matrix(CJ_SC_diff_data[,c(2,3,4)])/sqrt(Nsub_2)),
#           length=0,lwd=lwdgr, col='cyan4')
# error.bar(1.2:3.2,colMeans(CJ_SC_diff_data[,c(8,9,10)]),colSds(as.matrix(CJ_SC_diff_data[,c(8,9,10)])/sqrt(Nsub_2)),
#           length=0,lwd=lwdgr, col='brown3')
# 
# 
# 
# Plot confidence prediction ----------------------------------------------
par(mar=c(0,0,0,0))
plot.new()
text(.5,.75, labels="G. Model Predictions",cex = cex_legend+.5,font=2)
par(mar=c(5,5,0,2)+0.1)

Simuls$cj <- Simuls$cj_cont
## Experiment 1
#Aggregate conf for data
CJlow <- with(subset(Simuls,condition=="lowSC"),aggregate(cj,by=list(sub,coh),mean));names(CJlow) <- c('sub','coh','cj')
CJlow <- cast(CJlow,sub~coh,)
CJmed <- with(subset(Simuls,condition=="mediumSC"),aggregate(cj,by=list(sub,coh),mean));names(CJmed) <- c('sub','coh','cj')
CJmed <- cast(CJmed,sub~coh)
CJhigh <- with(subset(Simuls,condition=="highSC"),aggregate(cj,by=list(sub,coh),mean));names(CJhigh) <- c('sub','coh','cj')
CJhigh <- cast(CJhigh,sub~coh)


#aggregate cj for model
x <- CJlow[,c(2:4)];xmed <- CJmed[,c(2:4)];xhigh <- CJhigh[,c(2:4)]
n <- length(x)

x <- x[,c("hard","average","easy")];
xmed <- xmed[,c("hard","average","easy")];
xhigh <- xhigh[,c("hard","average","easy")]

stripchart(x,ylim=c(.45,1), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main=NULL, yaxt = 'n',family="A",xlab="Trial difficulty",ylab = "Confidence",cex.lab=cex_lab )
# mtext("Confidence",2,at=.75,line=2.5,cex=cex_lab);
# mtext("Trial difficulty",1,3,at=1,cex=cex_lab)
axis(1,at=0:(n-1),labels=names(x), cex.axis=cex_lab-1);
axis(2, seq(.5,1,.1), cex.axis=cex_lab-1)
means <- sapply(x, mean);n<- length(x)
for(i in seq(.5,1,length.out = 5)) abline(h=i,col="lightgrey",lty = "dashed")
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="brown3",lwd=lwddat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N1),lwd=lwdgr,col="brown3")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="cyan4",lwd=lwddat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N1),lwd=lwdgr,col="cyan4")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="darkgoldenrod3",lwd=lwddat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N1),lwd=lwdgr,col="darkgoldenrod3")


## Experiment 2
Simuls2$cj <- Simuls2$cj_cont
#Aggregate conf for Simuls2
CJlow <- with(subset(Simuls2,condition=="hard"),aggregate(cj,by=list(sub,coh),mean));names(CJlow) <- c('sub','coh','cj')
CJlow <- cast(CJlow,sub~coh,)
CJmed <- with(subset(Simuls2,condition=="average"),aggregate(cj,by=list(sub,coh),mean));names(CJmed) <- c('sub','coh','cj')
CJmed <- cast(CJmed,sub~coh)
CJhigh <- with(subset(Simuls2,condition=="easy"),aggregate(cj,by=list(sub,coh),mean));names(CJhigh) <- c('sub','coh','cj')
CJhigh <- cast(CJhigh,sub~coh)


#aggregate cj for model
x <- CJlow[,c(2:4)];xmed <- CJmed[,c(2:4)];xhigh <- CJhigh[,c(2:4)]
n <- length(x)

x <- x[,c("hard","average","easy")];
xmed <- xmed[,c("hard","average","easy")];
xhigh <- xhigh[,c("hard","average","easy")]

stripchart(x,ylim=c(.45,1), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main=NULL, yaxt = 'n',family="A",ylab="Confidence",xlab="Trial difficulty",cex.lab=cex_lab)
# mtext("Confidence",2,at=.75,line=2.5,cex=cex_lab);
# mtext("Trial difficulty",1,3,at=1,cex=cex_lab)
axis(1,at=0:(n-1),labels=names(x), cex.axis=cex_lab-1);
axis(2, seq(.5,1,.1), cex.axis=cex_lab-1)
means <- sapply(x, mean);n<- length(x)
for(i in seq(.5,1,.1)) abline(h=i,col="lightgrey",lty = "dashed")
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="brown3",lwd=lwddat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N1),lwd=lwdgr,col="brown3")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="cyan4",lwd=lwddat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N1),lwd=lwdgr,col="cyan4")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="darkgoldenrod3",lwd=lwddat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N1),lwd=lwdgr,col="darkgoldenrod3")

#'save
dev.off()
# EXP 1 Estimated parameters plot ----------------------------------------------
# Layout ------------------------------------------------------------------
jpeg(
  filename="fit_results1.jpeg",
  width=13,
  height=11,
  units="in",
  res=500)
# layout(matrix(c(1,2,1,3,1,4),ncol=3),heights = c(1.5,1))
layout(matrix(c(1,2,3,1,5,4),ncol=2),heights = c(.2,1,1))
par(mar=c(0,0,0,0))
plot.new()
text(.5,.75, labels="Experiment 1: fake feedback",cex = cex_legend+.5,font=2)
par(mar=c(5,5,4,0)+0.1)
# Plot Subjective drift -----------------------------------------------------
##Exp1
plot_drift <- with(df,aggregate(Vs,by=list(sub=sub,condition=condition),mean))
plot_drift <- cast(plot_drift,sub~condition)
plot_drift <- plot_drift[,c(3,4,2)] #Reorder columns to have hard -> easy
plot(colMeans(plot_drift),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylab='Subjective drift',ylim=c(0,.16),
     xlab="Feedback condition",xaxt='n', yaxt='n')
mtext("A.", at = .55, line = 3, cex = cex_title, font = 2)
segments(y0 = seq(0,.16,.04),y1 = seq(0,.16,.04),x0 = 0, x1 = Ncond, col = "lightgrey", lty = "dotted")
axis(1,1:Ncond,c("Negative","Average","Positive"),cex.axis=1.75)
axis(2,at=seq(0,.16,.04),cex.axis=1.75)
for(i in 1:Nsub_2){
  x <- jitter(1:Ncond,.2)
  lines(x,plot_drift[i,1:Ncond],lty=2,col=transp('grey'))
  points(x,plot_drift[i,1:Ncond],col="white", bg = transp('grey'),pch=21, cex = 2)
} 
points(colMeans(plot_drift),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_drift),colSds(plot_drift,na.rm=T)/sqrt(N1),lwd=3)

# Plot DDM parameters test phase EXP1 ------------------------------------------
par(mar=c(5,5,2,0)+0.1)
##Non-decision time
plot_ter <- with(param_1,aggregate(ter,by=list(sub=sub,condition=condition),mean))
plot_ter <- cast(plot_ter,sub~condition)
plot_ter <- plot_ter[,c(3,4,2)] #Reorder columns to have easy -> hard
plot(colMeans(plot_ter),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylim=c(.15,.8),ylab="Non-decision time",
     xlab="Feedback condition",xaxt='n',main="", cex.main = 2);
mtext("C.", at = .55, line = 1, cex = cex_title, font = 2)
axis(1,1:Ncond,c("Negative","Average","Positive"),cex.axis=1.75)
segments(y0 = seq(.2,.8,.1),y1 = seq(.2,.8,.1),x0 = 0, x1 = Ncond, col = "lightgrey", lty = "dotted")
for(i in 1:Nsub_2){
  x <- jitter(1:Ncond,.2)
  lines(x,plot_ter[i,1:Ncond],lty=2,col=transp('grey'))
  points(x,plot_ter[i,1:Ncond],col="white", bg = transp('grey'),pch=21, cex = 2)
} 
points(colMeans(plot_ter),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_ter),colSds(plot_ter,na.rm=T)/sqrt(N1),lwd=3,length=0)

##Bound
plot_bound <- with(param_1,aggregate(bound,by=list(sub=sub,condition=condition),mean))
plot_bound <- cast(plot_bound,sub~condition)
plot_bound <- plot_bound[,c(3,4,2)] #Reorder columns to have easy -> hard
plot(colMeans(plot_bound),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylim=c(0,.2),ylab="Bound",
     xlab="Feedback condition",xaxt='n',main="", cex.main = 2);
mtext("D.", at = .55, line = 1, cex = cex_title, font = 2)
segments(y0 = seq(0,.2,.05),y1 = seq(0,.2,.05),x0 = 0, x1 = Ncond, col = "lightgrey", lty = "dotted")
axis(1,1:Ncond,c("Negative","Average","Positive"),cex.axis=1.75)
for(i in 1:Nsub_2){
  x <- jitter(1:Ncond,.2)
  lines(x,plot_bound[i,1:Ncond],lty=2,col=transp('grey'))
  points(x,plot_bound[i,1:Ncond],col="white", bg = transp('grey'),pch=21, cex = 2)
} 
points(colMeans(plot_bound),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_bound),colSds(plot_bound,na.rm=T)/sqrt(N1),lwd=3,length=0)

##Drift interaction
plot_drift_minus <- with(subset(param_1,difflevel=="hard"),
                         aggregate(drift,by=list(sub=sub,condition=condition),mean))
plot_drift_minus <- cast(plot_drift_minus,sub~condition)
plot_drift_minus <- plot_drift_minus[,c(3,4,2)] #Reorder columns to have easy -> hard
plot(colMeans(plot_drift_minus),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,xlim=c(.8,Ncond+.2),
     ylim=c(min(plot_drift_minus),.27),ylab="Drift rate",xlab="Feedback condition",xaxt='n');
axis(1,1:Ncond,c("Negative","Average","Positive"),cex.axis=1.75)
mtext("B.", at = .55, line = 1, cex = cex_title, font = 2)
segments(y0 = seq(0,.25,.05),y1 = seq(0,.25,.05),x0 = 0, x1 = Ncond, col = "lightgrey", lty = "dotted")
points(colMeans(plot_drift_minus),type='b',lwd=5,col="darkolivegreen",lty="dashed")
error.bar(1:Ncond,colMeans(plot_drift_minus),
          colSds(plot_drift_minus,na.rm=T)/sqrt(N1),lwd=3,length=0,col="darkolivegreen")

plot_drift_control <- with(subset(param_1,difflevel=="average"),
                           aggregate(drift,by=list(sub=sub,condition=condition),mean))
plot_drift_control <- cast(plot_drift_control,sub~condition)
plot_drift_control <- plot_drift_control[,c(3,4,2)] #Reorder columns to have easy -> hard
points(colMeans(plot_drift_control),type='b',lwd=5,col="darkolivegreen3",lty="dotdash")
error.bar(1:Ncond,colMeans(plot_drift_control),
          colSds(plot_drift_control,na.rm=T)/sqrt(N1),lwd=3,length=0,col="darkolivegreen3")

plot_drift_plus <- with(subset(param_1,difflevel=="easy"),
                        aggregate(drift,by=list(sub=sub,condition=condition),mean))
plot_drift_plus <- cast(plot_drift_plus,sub~condition)
plot_drift_plus <- plot_drift_plus[,c(3,4,2)] #Reorder columns to have easy -> hard
points(colMeans(plot_drift_plus),type='b',lwd=5,col="darkolivegreen1",)
error.bar(1:Ncond,colMeans(plot_drift_plus),
          colSds(plot_drift_plus,na.rm=T)/sqrt(N1),lwd=3,length=0,col="darkolivegreen1")
par(xpd=T)
legend("top",border=F,legend=c("Hard","Average","Easy"),lwd=3, horiz = T, inset = c(0,-.065),
       col=c("darkolivegreen","darkolivegreen3","darkolivegreen1"),bty="n",cex=2,
       title = "Trial difficulty", lty = c("dashed","dotdash","solid"))
par(xpd=F)

#'save
dev.off()
# EXP 2 Estimated parameters plot -----------------------------------------
# Layout ------------------------------------------------------------------
jpeg(
  filename="fit_results2.jpeg",
  width=13,
  height=11,
  units="in",
  res=500)
layout(matrix(c(1,2,3,1,5,4),ncol=2),heights = c(.2,1,1))
par(mar=c(0,0,0,0))
plot.new()
text(.5,.75, labels="Experiment 2: training condition",cex = cex_legend+.5,font=2)
par(mar=c(5,5,4,0)+0.1)
# Plot Subjective drift -----------------------------------------------------
##Exp1
plot_drift <- with(df2,aggregate(Vs,by=list(sub=sub,condition=condition),mean))
plot_drift <- cast(plot_drift,sub~condition)
plot_drift <- plot_drift[,c(4,2,3)] #Reorder columns to have hard -> easy
plot(colMeans(plot_drift),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylab='Subjective drift',ylim=c(0,.5),
     xlab="Training condition",xaxt='n', yaxt='n')
mtext("A.", at = .55, line = 3, cex = cex_title, font = 2)
segments(y0 = seq(0,.5,.1),y1 = seq(0,.5,.1),x0 = 0, x1 = Ncond, col = "lightgrey", lty = "dotted")
axis(1,1:Ncond,c("Difficult","Medium","Easy"),cex.axis=1.75)
axis(2,at=seq(0,.5,.1),cex.axis=1.75)
for(i in 1:Nsub_2){
  x <- jitter(1:Ncond,.2)
  lines(x,plot_drift[i,1:Ncond],lty=2,col=transp('grey'))
  points(x,plot_drift[i,1:Ncond],col="white", bg = transp('grey'),pch=21, cex = 2)
} 
points(colMeans(plot_drift),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_drift),colSds(plot_drift,na.rm=T)/sqrt(Nsub_2),lwd=3)

# Plot DDM parameters test phase EXP2 ------------------------------------------
par(mar=c(5,5,2,0)+0.1)
##Non-decision time
plot_ter <- with(param2,aggregate(ter,by=list(sub=sub,condition=condition),mean))
plot_ter <- cast(plot_ter,sub~condition)
plot_ter <- plot_ter[,c(4,2,3)] #Reorder columns to have easy -> hard
plot(colMeans(plot_ter),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylim=c(.15,.8),ylab="Non-decision time",
     xlab="Training condition",xaxt='n',main="", cex.main = 2);
mtext("C.", at = .55, line = 1, cex = cex_title, font = 2)
axis(1,1:Ncond,c("Difficult","Medium","Easy"),cex.axis=1.75)
segments(y0 = seq(.2,.8,.1),y1 = seq(.2,.8,.1),x0 = 0, x1 = Ncond, col = "lightgrey", lty = "dotted")
for(i in 1:Nsub_2){
  x <- jitter(1:Ncond,.2)
  lines(x,plot_ter[i,1:Ncond],lty=2,col=transp('grey'))
  points(x,plot_ter[i,1:Ncond],col="white", bg = transp('grey'),pch=21, cex = 2)
} 
points(colMeans(plot_ter),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_ter),colSds(plot_ter,na.rm=T)/sqrt(Nsub_2),lwd=3,length=0)

##Bound
plot_bound <- with(param2,aggregate(bound,by=list(sub=sub,condition=condition),mean))
plot_bound <- cast(plot_bound,sub~condition)
plot_bound <- plot_bound[,c(4,2,3)] #Reorder columns to have easy -> hard
plot(colMeans(plot_bound),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylim=c(.04,.12),ylab="Bound",
     xlab="Training condition",xaxt='n',main="", cex.main = 2);
mtext("D.", at = .55, line = 1, cex = cex_title, font = 2)
segments(y0 = seq(.04,.12,.02),y1 = seq(.04,.12,.02),x0 = 0, x1 = Ncond, col = "lightgrey", lty = "dotted")
axis(1,1:Ncond,c("Difficult","Medium","Easy"),cex.axis=1.75)
for(i in 1:Nsub_2){
  x <- jitter(1:Ncond,.2)
  lines(x,plot_bound[i,1:Ncond],lty=2,col=transp('grey'))
  points(x,plot_bound[i,1:Ncond],col="white", bg = transp('grey'),pch=21, cex = 2)
} 
points(colMeans(plot_bound),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_bound),colSds(plot_bound,na.rm=T)/sqrt(Nsub_2),lwd=3,length=0)

##Drift interaction
plot_drift_minus <- with(subset(param2,difflevel=="hard"),
                         aggregate(drift,by=list(sub=sub,condition=condition),mean))
plot_drift_minus <- cast(plot_drift_minus,sub~condition)
plot_drift_minus <- plot_drift_minus[,c(4,2,3)] #Reorder columns to have easy -> hard
plot(colMeans(plot_drift_minus),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,xlim=c(.8,Ncond+.2),
     ylim=c(min(plot_drift_minus),.27),ylab="Drift rate",xlab="Training condition",xaxt='n',yaxt='n');
axis(1,1:Ncond,c("Difficult","Medium","Easy"),cex.axis=1.75)
axis(2,seq(0,.25,.05),cex.axis=1.75)
mtext("B.", at = .55, line = 1, cex = cex_title, font = 2)
segments(y0 = seq(0,.25,.05),y1 = seq(0,.25,.05),x0 = 0, x1 = Ncond, col = "lightgrey", lty = "dotted")
points(colMeans(plot_drift_minus),type='b',lwd=5,col="darkolivegreen",lty="dashed")
error.bar(1:Ncond,colMeans(plot_drift_minus),
          colSds(plot_drift_minus,na.rm=T)/sqrt(Nsub_2),lwd=3,length=0,col="darkolivegreen")

plot_drift_control <- with(subset(param2,difflevel=="average"),
                           aggregate(drift,by=list(sub=sub,condition=condition),mean))
plot_drift_control <- cast(plot_drift_control,sub~condition)
plot_drift_control <- plot_drift_control[,c(4,2,3)] #Reorder columns to have easy -> hard
points(colMeans(plot_drift_control),type='b',lwd=5,col="darkolivegreen3",lty="dotdash")
error.bar(1:Ncond,colMeans(plot_drift_control),
          colSds(plot_drift_control,na.rm=T)/sqrt(Nsub_2),lwd=3,length=0,col="darkolivegreen3")

plot_drift_plus <- with(subset(param2,difflevel=="easy"),
                        aggregate(drift,by=list(sub=sub,condition=condition),mean))
plot_drift_plus <- cast(plot_drift_plus,sub~condition)
plot_drift_plus <- plot_drift_plus[,c(4,2,3)] #Reorder columns to have easy -> hard
points(colMeans(plot_drift_plus),type='b',lwd=5,col="darkolivegreen1",)
error.bar(1:Ncond,colMeans(plot_drift_plus),
          colSds(plot_drift_plus,na.rm=T)/sqrt(Nsub_2),lwd=3,length=0,col="darkolivegreen1")
par(xpd=T)
legend("top",border=F,legend=c("Hard","Average","Easy"),lwd=3, horiz = T, inset = c(0,-.07),
       col=c("darkolivegreen","darkolivegreen3","darkolivegreen1"),bty="n",cex=2,
       title = "Trial difficulty", lty = c("dashed","dotdash","solid"))
par(xpd=F)

#'save
dev.off()
## Correlation Objective/Subjective drift ====
par(mfrow=c(1,3))
cond_ordered_1 <- c("lowSC","mediumSC","highSC")
cond_ordered_2 <- c("hard","average","easy")
#' Exp 1
for (c in 1:Ncond) {
obj_drift <- subset(df,condition==cond_ordered_1[c])$Vo
subj_drift <- subset(df,condition==cond_ordered_1[c])$Vs
drift_range <- c(min(min(subj_drift),min(obj_drift)),max(max(subj_drift),max(obj_drift)))
plot(obj_drift~subj_drift,cex.axis=1.75,cex.lab=1.75,frame=F,pch=19, main= "Experiment 1",
     xlim=drift_range, ylim =drift_range );
print(cor.test(obj_drift,subj_drift));
abline(lm(obj_drift~subj_drift),lty=2)
mtext(paste(cond_ordered_1[c],'r = ',round(cor(obj_drift,subj_drift),3)))
}
#' Exp 2
for (c in 1:Ncond) {
  obj_drift <- subset(df2,condition==cond_ordered_2[c])$Vo
  subj_drift <- subset(df2,condition==cond_ordered_2[c])$Vs
  drift_range <- c(min(min(subj_drift),min(obj_drift)),max(max(subj_drift),max(obj_drift)))
  plot(obj_drift~subj_drift,cex.axis=1.75,cex.lab=1.75,frame=F,pch=19, main = "Experiment 2", 
       xlim=drift_range, ylim =drift_range );
  print(cor.test(obj_drift,subj_drift));
  abline(lm(obj_drift~subj_drift),lty=2)
  mtext(paste(cond_ordered_2[c],'r = ',round(cor(obj_drift,subj_drift),3)))
}
#Comparing both experiments
par(mfrow=c(1,2))
obj_drift <- df$Vo
subj_drift <- df$Vs
drift_range <- c(min(min(subj_drift),min(obj_drift)),max(max(subj_drift),max(obj_drift)))
plot(obj_drift~subj_drift,cex.axis=1.75,cex.lab=1.75,frame=F,pch=19, xlim=drift_range, ylim =drift_range );
print(cor.test(obj_drift,subj_drift));
abline(lm(obj_drift~subj_drift),lty=2)
mtext(paste('Experiment 1: r = ',round(cor(obj_drift,subj_drift),3)))
obj_drift <- df2$Vo
subj_drift <- df2$Vs
drift_range <- c(min(min(subj_drift),min(obj_drift)),max(max(subj_drift),max(obj_drift)))
plot(obj_drift~subj_drift,cex.axis=1.75,cex.lab=1.75,frame=F,pch=19, xlim=drift_range, ylim =drift_range );
print(cor.test(obj_drift,subj_drift));
abline(lm(obj_drift~subj_drift),lty=2)
mtext(paste('Experiment 2: r = ',round(cor(obj_drift,subj_drift),3)))
