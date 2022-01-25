rm(list=ls())
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
wd <- dirname(curdir)
datadir <- paste0(curdir,"/data")
library(reshape);library(effects);library(lmerTest);library(scales)
setwd(paste0(wd,"/functions"))
library(Rcpp);library(DEoptim)
library(prob);library(car)
sourceCpp("DDM_with_confidence_slow.cpp")
source('fastmerge.R')
source("quantilefit_function_DDMonly.R")
source("rw_createHM_driftsign.R")

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
cexkl <- 2.5;cexgr <- 2;lwdgr <- 3;
# Global Parameters --------------------------------------------------------------
## Heat map resolution
dt <- .001; ev_bound <- .5; ev_window <- dt*10; upperRT <- 5
ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
timesteps <- upperRT/dt

#List of heat maps
v_min <- .001; v_max <- .5; step <- .001
drifts <- seq(v_min,v_max,step)

nsim <- 500 #per drift/cond/participant
ntrial <- 120; nrepeat <- 20 #Vs fitting
# EXP 1 -------------------------------------------------------------------
## Pre-process ====
if (!(file.exists(paste0(datadir,'/dataexp1_helene_full.csv')))) {
  setwd(paste0(wd,"/realdata_fit/"))
  for(i in 1:50){ #Load raw data
    if(i == 1){
      Data <- read.csv(paste0('RealData_1A/selfconfidence1A_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
    }else{
      temp <- read.csv(paste0('RealData_1A/selfconfidence1A_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
      Data <- rbind(Data,temp)
    }
  }
  
  Data <- subset(Data,rt>150)
  
  Training <- subset(Data, running == "training")
  Data['response'] <- 0
  Data$response[Data$resp == "['n']"] <- 1
  Data <- subset(Data,running == "main")
  
  
  ## Diagnostic plot per participant and task + chance performance testing
  N <- length(unique(Data$sub)); subs <- unique(Data$sub); exclusion <- c()
  tasks <- unique(Data$task)
  par(mfrow=c(2,2))
  for(i in 1:N){
    for (t in tasks) {
      tempDat <- subset(Data,sub==subs[i]&task==t)
      acc_block <- with(tempDat,aggregate(cor,by=list(block=block),mean))
      bias_block <- with(tempDat,aggregate(response,by=list(block=block),mean))
      plot(acc_block,ylab="Acc (.) and bias (x)",frame=F,ylim=c(0,1));abline(h=.5,lty=2,col="grey")
      points(bias_block,pch=4)
      plot(tempDat$rt/1000,frame=F,col=c("black"),main=paste('subject',i,"task :",t),ylab="RT",ylim=c(0,5))
      plot(tempDat$cj,frame=F,col=c("black"),ylim=c(1,6),ylab="conf")
      plot(tempDat$RTconf,frame=F,col=c("black"),ylab="RT_conf")
      test <- binom.test(length(tempDat$cor[tempDat$cor==1]),n=length(tempDat$cor),alternative = "greater")
      print(paste("In t0, sub",subs[i],"p =", round(test$p.value,3),"compared to chance"))
      if (test$p.value > .05) {
        exclusion <- c(exclusion,subs[i])
      }
    }
  }
  
  # /!\ Need to check with Hélène about exclusion of participants
  # Data <- subset(Data,!(sub %in% exclusion))
  
  Data$response[Data$response==0] <- -1
  
  df <- Data[,c("sub","task","selfconf","difflevel","rt","response","cor","cj","RTconf")]
  names(df) <- c("sub","task","selfconf","coh","rt","resp","cor","cj","RTconf")
  df$rt <- df$rt/1000;df$RTconf <- df$RTconf/1000
  # write.csv(df,"dataexp1_helene_full.csv",row.names = FALSE)
  
  Training <- subset(Training,!(sub %in% exclusion))
  Training['response'] <- -1
  Training$response[Training$resp == "['n']"] <- 1
  df_train <- Training[,c("sub","task","selfconf","difflevel","rt","response","cor","cj","RTconf")]
  names(df_train) <- c("sub","task","selfconf","coh","rt","resp","cor","cj","RTconf")
  df_train$rt <- df_train$rt/1000;df_train$RTconf <- df_train$RTconf/1000
  setwd(datadir)
  # write.csv(df_train,"dataexp1_helene_training.csv",row.names = FALSE)
  
}

## Data load ====
setwd(datadir)
Data1 <- read.csv('dataexp1_helene_full.csv')
Data1_train <- read.csv("dataexp1_helene_training_full.csv")


subs1 <- sort(unique(Data1_train$sub)); N1 <- length(subs1) 
cond <- sort(unique(Data1_train$selfconf)); Ncond <- length(cond)
coh <- sort(unique(Data1_train$coh));Ndiff <- length(coh)


#Retrieve feedback given /!\ Should be within the pre-process
len_block <- 24
Data1_train$fb <- -99
block <- 1
for (i in seq(1,dim(Data1_train)[1],len_block)) {
  if (Data1_train[i,"selfconf"]=="lowSC") {
    if (block%%5==1) {
      Data1_train[i:(i+len_block-1),]$fb <- round(runif(1,min=.66,max=.69),2)
    }else{
      Data1_train[i:(i+len_block-1),]$fb <- round(runif(1,min=.53,max=.66),2)
    }
  }
  
  if (Data1_train[i,"selfconf"]=="mediumSC") {
    Data1_train[i:(i+len_block-1),]$fb <- round(runif(1,min=.69,max=.82),2)
  }
  
  if (Data1_train[i,"selfconf"]=="highSC") {
    if (block%%5==2) {
      Data1_train[i:(i+len_block-1),]$fb <- round(runif(1,min=.82,max=.84),2)
    }else{
      Data1_train[i:(i+len_block-1),]$fb <- round(runif(1,min=.85,max=.98),2)
    }
  }
  block <- block + 1
}

#Load fitted DDM parameters in the training phase + median confidence RT
setwd(wd)
bound_train <- matrix(NA,N1,Ncond);v_train <- matrix(NA,N1,Ncond);
ter_train <- matrix(NA,N1,Ncond); resid_train <- matrix(NA,N1,Ncond)
resid <- matrix(NA,N1,Ncond)
v2_train <- matrix(NA,N1,Ncond); v3_train <- matrix(NA,N1,Ncond)
conf_rt <- matrix(NA,N1,Ncond)
for (i in 1:N1) {
  tempAll <- subset(Data1_train,sub==subs1[i])
  for (c in 1:Ncond) {
    tempDat <- subset(tempAll,selfconf==cond[c])
    if(file.exists(paste0('realdata_fit/RealData_1A/fits/trainfit/trainfit',cond[c],i,'.Rdata'))){
      load(paste0('realdata_fit/RealData_1A/fits/trainfit/trainfit',cond[c],i,'.Rdata'))
    }
    else{ #if not, fit the model
      optimal_params <- DEoptim(chi_square_optim, # function to optimize
                                lower = c( 0, 0, 0, 5000, .1, .0025, 0,   1,0,0,0,0,0), # a,ter,z,ntrials,sigma,dt,t2time,vratio,alpha,beta,v
                                upper = c(.2, 2, 0, 5000, .1, .0025, 0, 1,0,0,.5,.5,.5), # a,ter,z,ntrials,sigma,dt,t2time,vratio,alpha,beta,v
                                observations = tempDat,control=c(itermax=1000,steptol=100,reltol=.001,NP=50), returnFit = 1, binning = F)
      results <- summary(optimal_params)
      #save individual results
      save(results, file=paste0('realdata_fit/RealData_1A/fits/trainfit/trainfit',cond[c],i,'.Rdata'))
    }
    bound_train[i,c] <- results$optim$bestmem[1]
    ter_train[i,c] <- results$optim$bestmem[2]
    v_train[i,c] <- results$optim$bestmem[11]
    v2_train[i,c] <- results$optim$bestmem[12]
    v3_train[i,c] <- results$optim$bestmem[13]
    resid_train[i,c] <- results$optim$bestval

    tempDat <- subset(Data1,sub==subs1[i]&selfconf==cond[c])
    conf_rt[i,c] <- median(tempDat$RTconf)
  }  
}

#Load fitted DDM parameters in the testing phase 
#/!\ Current parameters come from a fitting with alpha/beta 
setwd(paste0(wd,"/realdata_fit/RealData_1A"))
bound <- matrix(NA,N1,Ncond);v <- matrix(NA,N1,Ncond);ter <- matrix(NA,N1,Ncond)
conf_rt <- matrix(NA,N1,Ncond); resid <- matrix(NA,N1,Ncond)
#Adjust the number of drift parameters to the model loaded
v2 <- matrix(NA,N1,Ncond);v3 <- matrix(NA,N1,Ncond) 
for(i in 1:N1){
  for(c in 1:Ncond){
    print(paste('Running participant',i,'from',N1,"condition",c))
    # load(paste0('fits/results_sub_',subs[i],'_',cond[c],'.Rdata'))
    load(paste0('fits/exp1_vratio/results_sub_',subs1[i],'_',cond[c],'.Rdata'))
    bound[i,c] <- results$optim$bestmem[1]
    ter[i,c] <- results$optim$bestmem[2]
    conf_rt[i,c] <- results$optim$bestmem[7]
    v[i,c] <- results$optim$bestmem[11]
    v2[i,c] <-   results$optim$bestmem[12]
    v3[i,c] <-   results$optim$bestmem[13]
    resid[i,c] <- results$optim$bestval
  }
}
param_1 <- data.frame(drift = c(v,v2,v3),bound=rep(bound,3),ter=rep(ter,3),
                      sub=rep(subs1,3*Ncond),
                      condition=rep(cond,each=N1,length.out=N1*Ncond*3),
                      difflevel=rep(coh,each=N1*Ncond),exp=2,resid=rep(resid,3))

#Aggregate train and test
bounds <- data.frame(bound = c(bound_train,bound),phase = rep(c("train","main"),each=length(bound)), sub = rep(subs1,Ncond*2),condition = rep(cond,each = N1,length.out=N1*Ncond*2))
ters <- data.frame(ter = c(ter_train,ter),phase = rep(c("train","main"),each=length(ter)), sub = rep(subs1,Ncond*2),condition = rep(cond,each = N1,length.out=N1*Ncond*2))
vs <- data.frame(v = c(v_train,v2_train,v3_train,v,v2,v3),phase = rep(c("train","main"),each=length(v)*Ndiff), sub = rep(subs1,Ncond*2*Ndiff),condition = rep(cond,each = N1,length.out=N1*Ncond*2*Ndiff),difficulty=rep(coh,each=N1*Ncond,length.out=N1*Ncond*2))

## Fit subjective drift ====
setwd(wd)

cost_conf <- matrix(NA,nrow=nrepeat,ncol=length(drifts))
if (!(file.exists(paste0(datadir,"/means_exp1_full.Rdata")))) {
  means <- matrix(NA,nrow=Ncond*N1,ncol=length(drifts)) 
  stds <- matrix(NA,nrow=Ncond*N1,ncol=length(drifts))
  for (s in 1:N1) {
    print(paste("Running participant",s,"of",N1))
    for (c in 1:Ncond) {
      tempDat <- subset(Data1_train,sub==subs1[s]&selfconf==cond[c])
      temp <- chi_square_optim(c(bound_train[s,c],ter_train[s,c],0,ntrial*nrepeat/Ndiff,.1,.0025,conf_rt[s,c],1,0,0,v_train[s,c],v2_train[s,c],v3_train[s,c]),NULL,0,binning=F)
      
      #match to the heatmap
      temp$closest_evdnc2 <- match.closest(temp$evidence2,ev_mapping)
      temp$temprt2 <- temp$rt2;
      temp$temprt2[temp$temprt2>5] <- 5 #heatmap doesn't go higher
      temp$temprt2 <- temp$temprt2*timesteps/5 #scale with the heatmap, between 0 and 2000
      
      temp$Nrep <- rep(1:nrepeat,each=ntrial/Ndiff,length.out=ntrial*nrepeat)        
      
      for (d in 1:length(drifts)) {
        print(paste("drift",d))
        load(paste0("heatmaps/hm_",drifts[d],"_filled.Rdata"))
        hm_low <- output$lower; hm_up <- output$upper
        hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)

        temp[temp$resp==1,]$cj <- hmvec_up[(temp[temp$resp==1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==1,]$temprt2)]
        temp[temp$resp==-1,]$cj <- hmvec_low[(temp[temp$resp==-1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==-1,]$temprt2)]
        

        for (i in 1:nrepeat) {
          pred_sample <- subset(temp,Nrep==i)
          diff <- sum((tempDat$fb - pred_sample$cj)^2)
          cost_conf[i,d] <- diff
        }
      }
      means[s+N*(c-1),] <- colMeans(cost_conf)
      stds[s+N*(c-1),] <- colSds(cost_conf)
    }
  }
  setwd(datadir)
  # save(means,file="means_exp1_full.Rdata")
  # save(stds,file="stds_exp1_full.Rdata")
}else{
  setwd(datadir)
  load("means_exp1_full.Rdata")
  load("stds_exp1_full.Rdata")
}


result <- sapply(seq(nrow(means)),function(i) {
  j <- which.min(means[i,])
  c(j)
})


Vs <- drifts[result]
Vs_matrix <- matrix(Vs,nrow=N1,ncol=Ncond)
df <- data.frame(Vs=Vs,bound = c(bound_train), ter = c(ter_train),Vo=c(v_train),sub=rep(subs1,Ncond),condition=rep(cond,each=N1))


## Plot cost for each drift ====
# 
# 
# n <- ncol(means)
# N <- length(subs)
# for (i in 1:nrow(means)) {
#   plot(0:(n-1),means[i,],type='b',lty=2,cex=1,lwd=lwdgr,pch=16,col=rgb(0,0,0,.5), xlim=c(-.05,n+-1),frame=F,xaxt="n",ylab="",xlab="",cex.axis=1.25,main=paste("Participant",subs[i%%N],"Condition",cond[(i%/%N)+1]))
#   polygon(c(0:(n-1),(n-1):0),c(means[i,] + (stds[i,]/sqrt(nrepeat)),(means[i,] - stds[i,]/sqrt(nrepeat))[n:1]),
#           border=F,col=rgb(0,0,0,.2))
#   mtext("Cost value",2,at=mean(c(max(means[i,]),min(means[i,]))),line=3,cex=1.5);axis(1,at=round(seq(0,n-1,length.out = 6)),labels=round(drifts[round(seq(0,n,length.out = 6))+1],2), cex.axis=1.25);mtext("Confidence",1,2.5,at=n/2,cex=1.5)
# }
# 
## Generate model prediction ====
setwd(wd)
rm(Simuls)
for(i in 1:N1){
  print(paste('simulating',i,'from',N1))
  for(c in 1:Ncond){
    load(paste0("heatmaps/hm_",Vs_matrix[i,c],"_filled.Rdata"))
    hm_low <- output$lower; hm_up <- output$upper
    hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
    temp <- chi_square_optim(c(bound[i,c],ter[i,c],0,nsim,.1,.0025,conf_rt[i,c],1,0,0,v[i,c],v2[i,c],v3[i,c]),NULL,0)

    #match to the heatmap
    temp$closest_evdnc2 <- match.closest(temp$evidence2,ev_mapping)
    temp$temprt2 <- temp$rt2;
    temp$temprt2[temp$temprt2>5] <- 5 #heatmap doesn't go higher
    temp$temprt2 <- temp$temprt2*timesteps/5 #scale with the heatmap, between 0 and 2000

    temp[temp$resp==1,]$cj <- hmvec_up[(temp[temp$resp==1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==1,]$temprt2)]
    temp[temp$resp==-1,]$cj <- hmvec_low[(temp[temp$resp==-1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==-1,]$temprt2)]
    
    temp$cj_cont <- temp$cj
    
    #/!\ What is this ?
    tempDat <- subset(Data1,sub==subs1[i]&selfconf==cond[c])
    for(k in 1:6){
      temp$cj[temp$cj < quantile(temp$cj,probs=sum(tempDat$cj==k)/dim(tempDat)[1])] <- k
    }
    temp$cj[temp$cj<1] <- 6 #put the extremes to six
    

    if(!exists('Simuls')){ Simuls <- cbind(temp,rep(cond[c],nsim),rep(subs1[i],nsim))
    }else{ Simuls <- rbind(Simuls,cbind(temp,rep(cond[c],nsim),rep(subs1[i],nsim)))
    }
  }
}
Simuls <- data.frame(Simuls);names(Simuls) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','closest_evdnc2',"temprt2",'cj_cont','condition','sub')

#Linear scaling into confidence judgments (deprecated)
Simuls$cj_bin <- as.numeric(cut(Simuls$cj_cont,breaks=seq(0,1,length.out = 7),include.lowest = TRUE))
Simuls$cj_scaled <- Simuls$cj #/!\ Delete ?

Simuls$cj <- Simuls$cj_bin

coherences <- sort(unique(Data1$coh))
Simuls$coh <- 0
for (i in 1:N1) {
  for(d in 1:length(coherences)) Simuls$coh[Simuls$sub==subs1[i] & Simuls$drift %in% c(unique(subset(Simuls,sub==subs1[i])$drift)[d],unique(subset(Simuls,sub==subs1[i])$drift)[d+3],unique(subset(Simuls,sub==subs1[i])$drift)[d+6])] <- coherences[d] #recode drift to coherence
}

#/!\ Is this the linear regression scaling ?
# Conf1_p <- with(Simuls,aggregate(cj,by=list(condition),mean))
# names(Conf1_p) <- c('condition','sub','conf')
# Conf1_p <- cast(Conf1_p,condition~sub)
# for(i in 1:N) Conf1_p_all[,i] <- predict(lm(as.numeric(CJ0[i,1:10])~Conf1_p_all[,i]))
# EXP 2 -------------------------------------------------------------------
## Data Load ====
setwd(paste0(wd,"/realdata_fit/RealData_1B"))
Data2 <- read.csv('dataexp2_helene.csv')

setwd(datadir)
Data2_train <- read.csv("dataexp2_helene_training.csv")

# exclusion <- c(50) #Fit didn't work
# Data2 <- subset(Data2,!(sub %in% exclusion))
# Data2_train <- subset(Data2_train,!(sub %in% exclusion))

#Retrieve feedback given
len_block <- 24
Data2_train$fb <- -99
for (i in seq(1,dim(Data2_train)[1],len_block)) {
  Data2_train[i:(i+len_block-1),]$fb <- round(mean(Data2_train[i:(i+len_block-1),]$cor),2)
}
Data2_train[Data2_train$fb<.5,]$fb <- .5 #Actual feedback was "lower than 50%"

subs <- sort(unique(Data2_train$sub)); N <- length(subs)
cond <- sort(unique(Data2_train$traindiffcond)); Ncond <- length(cond)
coh <- sort(unique(Data2_train$coh));Ndiff <- length(coh)

#Load fitted train DDM parameters
bound_train <- matrix(NA,N,Ncond);v_train <- matrix(NA,N,Ncond);ter_train <- matrix(NA,N,Ncond); resid_train <- matrix(NA,N,Ncond)
resid <- matrix(NA,N,Ncond)
v2 <- matrix(NA,N,Ncond); v3 <- matrix(NA,N,Ncond)
conf_rt <- matrix(NA,N,Ncond)

setwd(wd)
for (i in 1:N) {
  tempAll <- subset(Data2_train,sub==subs[i])
  for (c in 1:Ncond) {
    tempDat <- subset(tempAll,traindiffcond==cond[c])
    if(file.exists(paste0('realdata_fit/RealData_1B/fits/trainfit/trainfit',cond[c],i,'.Rdata'))){
      load(paste0('realdata_fit/RealData_1B/fits/trainfit/trainfit',cond[c],i,'.Rdata'))
    }
    else{ #if not, fit the model
      optimal_params <- DEoptim(chi_square_optim, # function to optimize
                                 # a,ter,z,ntrials,sigma,dt,t2time,vratio,alpha,beta,v
                                 lower = c( 0, 0, 0, 5000, .1, .0025, 0,   1,0,0,0), 
                                 upper = c(.2, 2, 0, 5000, .1, .0025, 0, 1,0,0,.5),
                                 observations = tempDat,returnFit = 1, binning = F,
                                 control=c(itermax=1000,steptol=100,reltol=.001,NP=30))
      results <- summary(optimal_params)
      #save individual results
      save(results, file=paste0('realdata_fit/RealData_1B/fits/trainfit/trainfit',cond[c],i,'.Rdata'))
    }
    bound_train[i,c] <- results$optim$bestmem[1]
    ter_train[i,c] <- results$optim$bestmem[2]
    v_train[i,c] <- results$optim$bestmem[11]
    resid_train[i,c] <- results$optim$bestval
    
    #Fit DDM param2eters for the testing phase
    tempDat <- subset(Data2,sub==subs[i]&traindiffcond==cond[c])
    conf_rt[i,c] <- median(tempDat$RTconf)
  }
}

setwd(paste0(wd,"/realdata_fit/RealData_1B"))
bound <- matrix(NA,N,Ncond);v <- matrix(NA,N,Ncond);ter <- matrix(NA,N,Ncond);
conf_rt <- matrix(NA,N,Ncond); resid <- matrix(NA,N,Ncond)
#Adjust the number of drift parameters to the model loaded
v2 <- matrix(NA,N,Ncond);v3 <- matrix(NA,N,Ncond);
for(i in 1:N){
  for(c in 1:Ncond){
    print(paste('Running participant',i,'from',N,"condition",c))
    # load(paste0('fits/results_sub_',subs[i],'_',cond[c],'.Rdata'))
    load(paste0('fits/exp2_vratio/results_sub_',subs[i],'_',cond[c],'.Rdata'))
    # plot(results$member$bestvalit[1:results$optim$iter],ylab='Goal function',ylim=c(-1,1),frame=F,type='l',main=paste('sub',i,'condition',cond[c]))
    bound[i,c] <- results$optim$bestmem[1]
    ter[i,c] <- results$optim$bestmem[2]
    conf_rt[i,c] <- results$optim$bestmem[7]
    v[i,c] <- results$optim$bestmem[11]
    v2[i,c] <-   results$optim$bestmem[12]
    v3[i,c] <-   results$optim$bestmem[13]
    resid[i,c] <- results$optim$bestval
  }
}

param2 <- data.frame(drift = c(v,v2,v3),bound=rep(bound,3),ter=rep(ter,3),
                     sub=rep(subs,3*Ncond),
                     condition=rep(cond,each=N,length.out=N*Ncond*3),
                     difflevel=rep(coh,each=N*Ncond),exp=2,resid=rep(resid,3))

## Fit subjective drift ====
setwd(wd)
#Generate model simulations
Ndiff <- 1 #Only one difficulty level in the training phase
Ncond <- length(cond)
cost_conf <- matrix(NA,nrow=nrepeat,ncol=length(drifts))
means <- matrix(NA,nrow=Ncond*N,ncol=length(drifts)); stds <- matrix(NA,nrow=Ncond*N,ncol=length(drifts))
if (!(file.exists(paste0(datadir,"/means_exp2.Rdata")))){
  for (s in 1:N) {
    print(paste("Running participant",s,"of",N))
    for (c in 1:Ncond) {
      tempDat <- subset(Data2_train,sub==subs[s]&traindiffcond==cond[c])
      
      temp <- chi_square_optim(c(bound_train[s,c],ter_train[s,c],0,ntrial*nrepeat/Ndiff,.1,.0025,conf_rt[s,c],1,0,0,v_train[s,c]),NULL,0,binning=F)
      
      #match to the heatmap
      temp$closest_evdnc2 <- match.closest(temp$evidence2,ev_mapping)
      temp$temprt2 <- temp$rt2;
      temp$temprt2[temp$temprt2>5] <- 5 #heatmap doesn't go higher
      temp$temprt2 <- temp$temprt2*timesteps/5 #scale with the heatmap, between 0 and 2000
      
      temp$Nrep <- rep(1:nrepeat,each=ntrial/Ndiff, length.out = ntrial*nrepeat*Ndiff)
      
      for (d in 1:length(drifts)) {
        print(paste("drift",d))
        load(paste0("heatmaps/hm_",drifts[d],"_filled.Rdata"))
        hm_low <- output$lower; hm_up <- output$upper
        hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
        
        temp[temp$resp==1,]$cj <- hmvec_up[(temp[temp$resp==1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==1,]$temprt2)]
        temp[temp$resp==-1,]$cj <- hmvec_low[(temp[temp$resp==-1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==-1,]$temprt2)]
        
        for (i in 1:nrepeat) {
          pred_sample <- subset(temp,Nrep==i)
          diff <- sum((tempDat$fb - pred_sample$cj)^2)
          cost_conf[i,d] <- diff
        }
        
        means[s+N*(c-1),] <- colMeans(cost_conf)
        stds[s+N*(c-1),] <- colSds(cost_conf)
      }
    }
  }
  setwd(datadir)
  # save(means,file="means_exp2.Rdata")
  # save(stds,file="stds_exp2.Rdata")
}else{
  setwd(datadir)
  load("means_exp2.Rdata")
  load("stds_exp2.Rdata")
}


result <- sapply(seq(nrow(means)),function(i) {
  j <- which.min(means[i,])
  # c(paste(i, j, sep='/'), means[i,j])
  c(j)
})


Vs <- drifts[result]
Vs_matrix <- matrix(Vs,nrow=N,ncol=Ncond)
df2 <- data.frame(drift=Vs,bound = c(bound_train), ter = c(ter_train),v=c(v_train),sub=rep(subs,3),condition=rep(cond,each=N))
## Plot cost for each drift ====
# error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
#   if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
#     stop("vectors must be same length")
#   arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
# }
# cexkl <- 2.5;cexgr <- 2;lwdgr <- 3;
# 
# 
# n <- ncol(means)
# N <- length(subs)
# for (i in 1:nrow(means)) {
#   plot(0:(n-1),means[i,],type='b',lty=2,cex=1,lwd=lwdgr,pch=16,col=rgb(0,0,0,.5), xlim=c(-.05,n+-1),frame=F,xaxt="n",ylab="",xlab="",cex.axis=1.25,main=paste("Participant",subs[i%%N],"Condition",cond[(i%/%N)+1]))
#   polygon(c(0:(n-1),(n-1):0),c(means[i,] + (stds[i,]/sqrt(nrepeat)),(means[i,] - stds[i,]/sqrt(nrepeat))[n:1]),
#           border=F,col=rgb(0,0,0,.2))
#   mtext("Cost value",2,at=max(means[i,]/2),line=3,cex=1.5);axis(1,at=round(seq(0,n-1,length.out = 6)),labels=round(drifts[round(seq(0,n,length.out = 6))+1],2), cex.axis=1.25);mtext("Confidence",1,2.5,at=n/2,cex=1.5)
# }
# # stripchart(x, ylim=c(0,1), xlim=c(-.05,n+-1), vertical = TRUE, col="white",frame=F,xaxt='n',
#            # main="",cex.axis=1.25)
# # means <- sapply(x, mean, na.rm=T);n<- length(x)
# plot(x=0:(n-1),y=means,type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=16,col=rgb(0,0,0,.5),ylim=c(0,1), xlim=c(-.05,n+-1),frame=F,xaxt="n")
# polygon(c(0:(n-1),(n-1):0),c(means + (stds/sqrt(N)),(means - stds/sqrt(N))[n:1]),
#         border=F,col=rgb(0,0,0,.2))
# mtext("Cost value",2,at=.5,line=3,cex=1.5);axis(1,at=0:(n-1),labels=names(x), cex.axis=1.25);mtext("Confidence",1,2.5,at=2.5,cex=1.5)
# # lines(0:(n-1),means,type='p',pch=16,cex=cexkl)
# # error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05)
# 
# 
## Generate model prediction ====
setwd(wd)
rm(Simuls2)
for(i in 1:N){
  print(paste('simulating',i,'from',N))
  for(c in 1:Ncond){
    tempdat <- subset(Data2, sub==subs[i] & traindiffcond==cond[c])
    load(paste0("heatmaps/hm_",Vs_matrix[i,c],"_filled.Rdata"))
    hm_low <- output$lower; hm_up <- output$upper
    hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
    temp <- chi_square_optim(c(bound[i,c],ter[i,c],0,nsim,.1,.0025,conf_rt[i,c],1,0,0,v[i,c],v2[i,c],v3[i,c]),NULL,0)
    
    #match to the heatmap
    temp$closest_evdnc2 <- match.closest(temp$evidence2,ev_mapping)
    temp$temprt2 <- temp$rt2;
    temp$temprt2[temp$temprt2>5] <- 5 #heatmap doesn't go higher
    temp$temprt2 <- temp$temprt2*timesteps/5 #scale with the heatmap, between 0 and 2000
    
    temp[temp$resp==1,]$cj <- hmvec_up[(temp[temp$resp==1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==1,]$temprt2)]
    temp[temp$resp==-1,]$cj <- hmvec_low[(temp[temp$resp==-1,]$closest_evdnc2-1)*timesteps+round(temp[temp$resp==-1,]$temprt2)]
    
    temp$cj_cont <- temp$cj
    
    #/!\ What is this ?
    tempDat <- subset(Data2,sub==subs[i]&traindiffcond==cond[c])
    for(k in 1:6){
      temp$cj[temp$cj < quantile(temp$cj,probs=sum(tempDat$cj==k)/dim(tempDat)[1])] <- k
    }
    temp$cj[temp$cj<1] <- 6 #put the extremes to six
    
    if(!exists('Simuls2')){ Simuls2 <- cbind(temp,rep(cond[c],nsim),rep(subs[i],nsim))
    }else{ Simuls2 <- rbind(Simuls2,cbind(temp,rep(cond[c],nsim),rep(subs[i],nsim)))
    }
  }
}
Simuls2 <- data.frame(Simuls2);names(Simuls2) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','closest_evdnc2',"temprt2",'cj_cont','condition','sub')

#Linear scaling into confidence ratings (deprecated)
Simuls2$cj_bin <- as.numeric(cut(Simuls2$cj_cont,breaks=seq(0,1,length.out = 7),include.lowest = TRUE))
Simuls2$cj_scaled <- Simuls2$cj #/!\ Delete ?

coherences <- sort(unique(Data2$coh))
Simuls2$coh <- 0
for (i in 1:N) {
  for(d in 1:length(coherences)) Simuls2$coh[Simuls2$sub==subs[i] & Simuls2$drift %in% c(unique(subset(Simuls2,sub==subs[i])$drift)[d],unique(subset(Simuls2,sub==subs[i])$drift)[d+3],unique(subset(Simuls2,sub==subs[i])$drift)[d+6])] <- coherences[d] #recode drift to coherence
}

# Remove ? ----------------------------------------------------------------
## Correlation Objective/Subjective drift ====
for (c in 1:Ncond) {
  obj_drift <- subset(df2,condition==cond[c])$v
  subj_drift <- subset(df2,condition==cond[c])$drift
  drift_range <- c(min(min(subj_drift),min(obj_drift)),max(max(subj_drift),max(obj_drift)))
  plot(obj_drift~subj_drift,cex.axis=1.75,cex.lab=1.75,frame=F,pch=19, xlim=drift_range, ylim =drift_range );print(cor.test(obj_drift,subj_drift));abline(lm(obj_drift~subj_drift),lty=2);mtext(paste(cond[c],'r = ',round(cor(obj_drift,subj_drift),3)))
}
## Scale heatmap conf to reports ====
# 
# Conf_dat_all <- matrix(NA,nrow=N,ncol=9)
# Conf_sim_all <- matrix(NA,nrow=N,ncol=9)
# N1 <- length(subs1)
# for (s in 1:N1) {
#   temp_sim <- subset(Simuls,sub==subs1[s])
#   temp_dat <- subset(Data1,sub==subs1[s])
#   
#   cj_dat <- with(temp_dat,aggregate(cj,by=list(selfconf,coh),mean))
#   names(cj_dat) <- c('condition','coh','cj')
#   cj_dat <- cast(cj_dat,coh~condition)
#   cj_dat <- cj_dat[c(3,1,2),]
#   
#   cj_sim <- with(temp_sim,aggregate(cj_cont,by=list(condition,coh),mean))
#   names(cj_sim) <- c('condition','coh','cj')
#   cj_sim <- cast(cj_sim,coh~condition)
#   cj_sim <- cj_sim[c(3,1,2),]
#   
#   Conf_dat_all[s,] <- c(cj_dat[,3],cj_dat[,4],cj_dat[,2])
#   Conf_sim_all[s,] <- c(cj_sim[,3],cj_sim[,4],cj_sim[,2])
#   
#   Conf_sim_all[s,] <- predict(lm(Conf_dat_all[s,]~Conf_sim_all[s,]))
# }
# 
# Conf_dat_all2 <- matrix(NA,nrow=N,ncol=9)
# Conf_sim_all2 <- matrix(NA,nrow=N,ncol=9)
# 
# for (s in 1:N) {
#   temp_sim <- subset(Simuls2,sub==subs[s])
#   temp_dat <- subset(Data2,sub==subs[s])
#   
#   cj_dat <- with(temp_dat,aggregate(cj,by=list(traindiffcond,coh),mean))
#   names(cj_dat) <- c('condition','coh','cj')
#   cj_dat <- cast(cj_dat,coh~condition)
#   cj_dat <- cj_dat[c(3,1,2),]
#   
#   cj_sim <- with(temp_sim,aggregate(cj_cont,by=list(condition,coh),mean))
#   names(cj_sim) <- c('condition','coh','cj')
#   cj_sim <- cast(cj_sim,coh~condition)
#   cj_sim <- cj_sim[c(3,1,2),]
#   
#   Conf_dat_all2[s,] <- c(cj_dat[,4],cj_dat[,2],cj_dat[,3])
#   Conf_sim_all2[s,] <- c(cj_sim[,4],cj_sim[,2],cj_sim[,3])
#   
#   Conf_sim_all2[s,] <- predict(lm(Conf_dat_all2[s,]~Conf_sim_all2[s,]))
# }
# 
# cj_dat <- with(Data1,aggregate(cj,by=list(sub,selfconf),mean))
# names(cj_dat) <- c('sub','condition','cj')
# cj_dat <- cast(cj_dat,sub~condition)
# cj_dat <- cj_dat[,c(3,4,2)]
# cj_dat <- as.matrix(cj_dat)
# cj_sim <- with(Simuls,aggregate(cj_cont,by=list(sub,condition),mean))
# names(cj_sim) <- c('sub','condition','cj')
# cj_sim <- cast(cj_sim,sub~condition)
# cj_sim <- cj_sim[,c(3,4,2)]; 
# cj_sim <- as.matrix(cj_sim)
# for (s in 1:N1) {
#   cj_sim[s,] <- predict(lm(cj_dat[s,]~cj_sim[s,]))
# }
# plot(colMeans(cj_dat),pch=16,cex=cexkl,ylim=c(4.2,5),xaxt='n',frame=F,ylab='',xlab='')
# polygon(c(1:3,3:1),c(colMeans(cj_sim,na.rm=T) + (colSds(cj_sim)/sqrt(N1)),(colMeans(cj_sim,na.rm=T) - colSds(cj_sim)/sqrt(N1))[3:1]),
#         border=F,col=rgb(0,0,0,.2))
# error.bar(1:3,colMeans(cj_dat),colSds(cj_dat,na.rm=T)/sqrt(N1),lwd=lwdgr,length=.05)
# mtext("Confidence",2,at=4.6,line=3,cex=1.75);axis(1,at=1:3,labels=c('Negative','Average','Positive'), cex.axis=1.5);mtext("Fake feedback condition",1,3,at=2,cex=1.75)
# 
# 
# cj_dat <- with(Data2,aggregate(cj,by=list(sub,traindiffcond),mean))
# names(cj_dat) <- c('sub','condition','cj')
# cj_dat <- cast(cj_dat,sub~condition)
# cj_dat <- cj_dat[,c(4,2,3)]
# cj_dat <- as.matrix(cj_dat)
# cj_sim <- with(Simuls2,aggregate(cj_cont,by=list(sub,condition),mean))
# names(cj_sim) <- c('sub','condition','cj')
# cj_sim <- cast(cj_sim,sub~condition)
# cj_sim <- cj_sim[,c(4,2,3)]; 
# cj_sim <- as.matrix(cj_sim)
# for (s in 1:N1) {
#   cj_sim[s,] <- predict(lm(cj_dat[s,]~cj_sim[s,]))
# }
# plot(colMeans(cj_dat),pch=16,cex=cexkl,ylim=c(4.2,5),xaxt='n',frame=F,ylab='',xlab='')
# polygon(c(1:3,3:1),c(colMeans(cj_sim,na.rm=T) + (colSds(cj_sim)/sqrt(N)),(colMeans(cj_sim,na.rm=T) - colSds(cj_sim)/sqrt(N))[3:1]),
#         border=F,col=rgb(0,0,0,.2))
# error.bar(1:3,colMeans(cj_dat),colSds(cj_dat,na.rm=T)/sqrt(N),lwd=lwdgr,length=.05)
# mtext("Confidence",2,at=4.6,line=3,cex=1.75);axis(1,at=1:3,labels=c('Hard','Average','Easy'), cex.axis=1.5);mtext("Training difficulty condition",1,3,at=2,cex=1.75)
# 
# 
# Plot Vs ~ Condition -----------------------------------------------------
##Exp1
plot_drift <- with(df,aggregate(Vs,by=list(sub=sub,condition=condition),mean))
plot_drift <- cast(plot_drift,sub~condition)
plot_drift <- plot_drift[,c(2,4,3)] #Reorder columns to have easy -> hard
plot(colMeans(plot_drift),frame=F,type='n',cex.lab=2,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylab='',ylim=c(min(plot_drift),max(plot_drift)),
     xlab="Fake Feedback",xaxt='n')
axis(1,1:Ncond,c("Positive","Average","Negative"),cex.axis=1.75)
mtext("Subjective drift",side = 2, line = 2.5, cex = 2)
for(i in 1:N) lines(1:Ncond,plot_drift[i,1:Ncond],type='b',lty=2,col="grey",pch=19)
points(colMeans(plot_drift),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_drift),colSds(plot_drift,na.rm=T)/sqrt(N),lwd=3,length=.05)

##Exp2
plot_drift <- with(df2,aggregate(Vs,by=list(sub=sub,condition=condition),mean));plot_drift <- cast(plot_drift,sub~condition)
plot_drift <- plot_drift[,c(3,2,4)] #Reorder columns to have easy -> hard
plot(colMeans(plot_drift),frame=F,type='n',cex.lab=2,cex.axis=1.75,xlim=c(.8,Ncond+.2),ylim=c(min(plot_drift),max(plot_drift)),ylab="",xlab="Training Difficulty",xaxt='n');axis(1,1:Ncond,c("Easy","Average","Hard"),cex.axis=1.75)
mtext("Subjective drift",side = 2, line = 2.5, cex = 2)
for(i in 1:N) lines(1:Ncond,plot_drift[i,1:Ncond],type='b',lty=2,col="grey",pch=19)
points(colMeans(plot_drift),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_drift),colSds(plot_drift,na.rm=T)/sqrt(N),lwd=3,length=.05)

# Mixed models ------------------------------------------------------------
## Exp1 ====
# #DDM train
# m <- lmer(Vs ~ cond + (1|sub),data=df); anova(m); with(df,aggregate(drift,by=list(cond=cond),mean))
# m <- lmer(Vo ~ cond + (1|sub),data=df); anova(m); with(df,aggregate(v,by=list(cond=cond),mean))
# m <- lmer(bound ~ cond + (1|sub),data=df); anova(m); with(df,aggregate(bound,by=list(cond=cond),mean))
# m <- lmer(ter ~ cond + (1|sub),data=df); anova(m); with(df,aggregate(ter,by=list(cond=cond),mean))
# 
# #Train vs Test
# m <- lmer(bound ~ phase*condition + (1|sub),data = bounds); anova(m)
# m <- lmer(v ~ phase*condition*difficulty + (1|sub),data = vs); anova(m)
# m <- lmer(ter ~ phase*condition + (1|sub),data = ters); anova(m)
# 
# ## Exp2 ====
# #DDM train
# m <- lmer(Vs ~ condition + (1|sub),data=df2); anova(m); with(df2,aggregate(drift,by=list(condition),mean))
# m <- lmer(Vo ~ condition + (1|sub),data=df2); anova(m); with(df2,aggregate(v,by=list(condition),mean))
# m <- lmer(bound ~ condition + (1|sub),data=df2); anova(m); with(df2,aggregate(bound,by=list(condition),mean))
# m <- lmer(ter ~ condition + (1|sub),data=df2); anova(m); with(df2,aggregate(ter,by=list(condition),mean))