##' 
##' 

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

error.bar <- function(x, y, upper, lower=upper, length=0,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
cexkl <- 1.5;cexgr <- 2;lwdgr <- 3; lwddat <- 2
# Global parameters --------------------------------------------------------------
## Heat map resolution
dt <- .001; ev_bound <- .5; ev_window <- dt*10; upperRT <- 5
ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
timesteps <- upperRT/dt

## List of heat maps
v_min <- .001; v_max <- .5; step <- .001
drifts <- seq(v_min,v_max,step)

nsim <- 500 # per drift/cond/participant
ntrial <- 120; nrepeat <- 20 # Vs fitting
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
cost_conf <- matrix(NA,nrow=nrepeat,ncol=length(drifts))
if (!(file.exists("means_exp1_full.Rdata"))) {
  means <- matrix(NA,nrow=Ncond*N1,ncol=length(drifts)) 
  stds <- matrix(NA,nrow=Ncond*N1,ncol=length(drifts))
  for (s in 1:N1) {
    for (cond in 1:Ncond) {
      print(paste("Running participant",s,"of",N1,"condition",cond))
      tempDat <- subset(Data1_train,sub==subs1[s]&selfconf==cond_1[cond])
      temp_par <- c(bound_train[s,cond],ter_train[s,cond],0,ntrial*nrepeat/Ndiff,
                    .1,.001,conf_rt[s,cond],1,v_train[s,cond],v2_train[s,cond],v3_train[s,cond])
      temp <- chi_square_optim_DDM(temp_par,observations=NULL,returnFit=0)
      
      # match to the heatmap
      temp$closest_evdnc2 <- match.closest(temp$evidence2,ev_mapping)
      temp$temprt2 <- temp$rt2;
      temp$temprt2[temp$temprt2>5] <- 5 #heatmap doesn't go higher (5 seconds)
      temp$temprt2 <- temp$temprt2*timesteps/5 #scale with the heatmap
      
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
          diff <- sum((tempDat$fb - pred_sample$cj)^2)
          cost_conf[i,d] <- diff
        }
        setTxtProgressBar(bar,d)
      }
      means[s+N1*(cond-1),] <- colMeans(cost_conf)
      stds[s+N1*(cond-1),] <- colSds(cost_conf)
    }
  }
  save(means,file="means_exp1_full.Rdata")
  save(stds,file="stds_exp1_full.Rdata")
}else{
  load("means_exp1_full.Rdata")
  load("stds_exp1_full.Rdata")
}


result <- sapply(seq(nrow(means)),function(i) {
  j <- which.min(means[i,])
  c(j)
})


Vs1 <- drifts[result]
Vs1_matrix <- matrix(Vs1,nrow=N1,ncol=Ncond)
df <- data.frame(Vs=Vs1,bound = c(bound_train), ter = c(ter_train),Vo=c(v_train),
                 sub=rep(subs1,Ncond),condition=rep(cond_1,each=N1))


## Generate model prediction ====
go_to("results")
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
#Generate model simulations
Ndiff <- 1 #Only one difficulty level in the training phase
cost_conf <- matrix(NA,nrow=nrepeat,ncol=length(drifts))
means <- matrix(NA,nrow=Ncond_2*Nsub_2,ncol=length(drifts))
stds <- matrix(NA,nrow=Ncond_2*Nsub_2,ncol=length(drifts))
if (!(file.exists("means_exp2.Rdata"))){
  for (s in 1:Nsub_2) {
    for (c in 1:Ncond_2) {
      print(paste("Running participant",s,"of",Nsub_2,"condition",c))
      tempDat <- subset(Data2_train,sub==subs_2[s]&traindiffcond==cond_2[c])
      
      temp <- chi_square_optim_DDM(c(bound_train[s,c],ter_train[s,c],0,ntrial*nrepeat/Ndiff,.1,.001,conf_rt[s,c],1,v_train[s,c]),
                                   observations = NULL, returnFit = 0)
      
      #match to the heatmap
      temp$closest_evdnc2 <- match.closest(temp$evidence2,ev_mapping)
      temp$temprt2 <- temp$rt2;
      temp$temprt2[temp$temprt2>5] <- 5 #heatmap doesn't go higher
      temp$temprt2 <- temp$temprt2*timesteps/5 #scale with the heatmap, between 0 and 2000
      
      temp$Nrep <- rep(1:nrepeat,each=ntrial/Ndiff, length.out = ntrial*nrepeat*Ndiff)
      
      bar <- txtProgressBar(0,length(drifts),style=3,char="#")
      for (d in 1:length(drifts)) {
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
        
        means[s+Nsub_2*(c-1),] <- colMeans(cost_conf)
        stds[s+Nsub_2*(c-1),] <- colSds(cost_conf)
        setTxtProgressBar(bar,d)
      }
    }
  }
  save(means,file="means_exp2.Rdata")
  save(stds,file="stds_exp2.Rdata")
}else{
  load("means_exp2.Rdata")
  load("stds_exp2.Rdata")
}


result <- sapply(seq(nrow(means)),function(i) {
  j <- which.min(means[i,])
  # c(paste(i, j, sep='/'), means[i,j])
  c(j)
})


Vs2 <- drifts[result]
Vs2_matrix <- matrix(Vs2,nrow=Nsub_2,ncol=Ncond_2)
df2 <- data.frame(Vs=Vs2,bound = c(bound_train), ter = c(ter_train),Vo=c(v_train),
                  sub=rep(subs_2,Ndiff),condition=rep(cond_2,each=Nsub_2))
## Generate model prediction ====
go_to("results")
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

# Export ------------------------------------------------------------------
go_to("results")
write.csv(Simuls,file = "model_prediction_exp1.csv")
write.csv(Simuls2,file = "model_prediction_exp2.csv")
# Stat tests ------------------------------------------------------------
# Exp1 ====
#DDM train
df$sub <- as.factor(df$sub)
m <- lm(Vs ~ condition,data=df); Anova(m);
m <- lmer(Vs ~ condition + (1|sub),data=df); Anova(m);
m <- aov(Vs ~ condition+ Error(sub/condition), data = df); summary(m)
m <- lmer(Vs ~condition + (1|sub), data = df); anova(m)
m <- lm(Vo ~ condition,data=df); Anova(m);
m <- lm(bound ~ condition,data=df); Anova(m);
m <- lm(ter ~ condition,data=df); Anova(m);

#DDM test
m <- lm(bound ~ condition,data=param_1); Anova(m);
m <- lm(ter ~ condition,data=param_1); Anova(m);
m <- lm(drift ~ condition*difflevel,data=param_1); Anova(m);

#Train vs Test
m <- lm(bound ~ phase*condition,data = bounds1); Anova(m)
m <- lm(v ~ phase*condition*difficulty,data = vs); Anova(m)
m <- lm(ter ~ phase*condition,data = ters1); Anova(m)
## Exp2 ====
#DDM train
df2$sub <- as.factor(df2$sub)
m <- aov(Vs ~ condition+ Error(sub/condition), data = df2); summary(m)
m <- lmer(Vs ~condition + (1|sub), data = df2); anova(m)
m <- lm(Vo ~ condition,data=df2); Anova(m);
m <- lm(bound ~ condition,data=df2); Anova(m);
m <- lm(ter ~ condition,data=df2); Anova(m)

#DDM test
m <- lm(bound ~ condition,data=param2); Anova(m);
m <- lm(ter ~ condition,data=param2); Anova(m);
m <- lm(drift ~ condition*difflevel,data=param2); Anova(m);

#Train vs Test
m <- lm(bound ~ phase*condition,data = bounds2); Anova(m)
m <- lm(ter ~ phase*condition,data = ters2); Anova(m)

Simuls$sub <- as.factor(Simuls$sub)
Simuls2$sub <- as.factor(Simuls2$sub)
m <- lmer(cj ~ condition*coh + (condition|sub), data = Simuls); anova(m)
m <- lmer(cj ~ condition*coh + (coh|sub), data = Simuls); anova(m)
m <- aov(cj ~ condition*coh+ Error(sub/condition), data = Simuls); summary(m)
m <- aov(cj ~ condition+ Error(sub/condition), data = Simuls2); summary(m)
# Plot Layout -------------------------------------------------------------
go_to("plot")

windowsFonts(A = windowsFont("Calibri")) 
par(family="A")

cex_axis <- 1.5; cex_legend <- 3 
jpeg(
  filename="results4.jpeg",
  width=13,
  height=18,
  units="in",
  res=500)
layout(matrix(c(1,3,7,9,10,1,5,7,9,10,2,4,8,9,11,2,6,8,9,11),ncol=4),heights = c(.4,1,2,.2,2))

#' Add legend for both experiments on top
par(mar=c(0,0,0,0))
plot.new()
legend("top",legend=c("Negative","Average","Positive"),
       title = "Feedback condition",pch=rep(16,3),bty = "n",inset=0, 
       cex = cex_legend,col=c("brown3","cyan4","darkgoldenrod3"), horiz = T)

plot.new()
legend("top",legend=c("Difficult","Medium","Easy"),
       title = "Training condition",pch=rep(16,3),bty = "n",inset=0, 
       cex = cex_legend,col=c("brown3","cyan4","darkgoldenrod3"), horiz = T)
par(mar=c(5,5,0,2)+0.1)
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
           frame=F,xaxt='n',main=NULL,yaxt='n')
axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);
axis(2,seq(0.6,1,.1),cex.axis=1.5)
mtext("Accuracy",2,at=.8,line=2.5,cex=cex_axis);
mtext("Trial difficulty",1,3,at=1,cex=cex_axis)
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
           frame=F,xaxt='n',main=NULL,yaxt='n')
axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);
axis(2,seq(0.6,1,.1),cex.axis=1.5)
mtext("Accuracy",2,at=.8,line=2.5,cex=cex_axis);
mtext("Trial difficulty",1,3,at=1,cex=cex_axis)
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
           frame=F,xaxt='n',main=NULL,yaxt='n')
axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);
axis(2,seq(0.6,1.1,.1),cex.axis=1.5)
mtext("RT (s)",2,at=.85,line=2.5,cex=cex_axis);
mtext("Trial difficulty",1,3,at=1,cex=cex_axis)
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
           frame=F,xaxt='n',main=NULL,yaxt='n')
axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);
axis(2,seq(0.6,1.1,.1),cex.axis=1.5)
mtext("RT (s)",2,at=.85,line=2.5,cex=cex_axis);
mtext("Trial difficulty",1,3,at=1,cex=cex_axis)
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
##Experiment 1
#Aggregate conf for data
CJlow <- with(subset(Data1,selfconf=="lowSC"),aggregate(cj,by=list(sub,coh),mean));names(CJlow) <- c('sub','coh','cj')
CJlow <- cast(CJlow,sub~coh,)
CJmed <- with(subset(Data1,selfconf=="mediumSC"),aggregate(cj,by=list(sub,coh),mean));names(CJmed) <- c('sub','coh','cj')
CJmed <- cast(CJmed,sub~coh)
CJhigh <- with(subset(Data1,selfconf=="highSC"),aggregate(cj,by=list(sub,coh),mean));names(CJhigh) <- c('sub','coh','cj')
CJhigh <- cast(CJhigh,sub~coh)


#aggregate cj for model
x <- CJlow[,c(2:4)];xmed <- CJmed[,c(2:4)];xhigh <- CJhigh[,c(2:4)]
n <- length(x)

x <- x[,c("hard","average","easy")];
xmed <- xmed[,c("hard","average","easy")];
xhigh <- xhigh[,c("hard","average","easy")]

stripchart(x,ylim=c(3.75,5.5), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main=NULL, yaxt = 'n',family="A")
axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);
axis(2, seq(4,5.5,.5), cex.axis=1.5)
mtext("Confidence",2,at=4.75,line=2.5,cex=cex_axis);
means <- sapply(x, mean);n<- length(x)
for(i in seq(4,5.5,.5)) abline(h=i,col="lightgrey",lty = "dashed")
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="brown3",lwd=lwddat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N1),lwd=lwdgr,col="brown3")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="cyan4",lwd=lwddat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N1),lwd=lwdgr,col="cyan4")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="darkgoldenrod3",lwd=lwddat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N1),lwd=lwdgr,col="darkgoldenrod3")

## Experiment 2
#Aggregate conf for data
CJlow <- with(subset(Data2,traindiffcond=="hard"),aggregate(cj,by=list(sub,coh),mean));names(CJlow) <- c('sub','coh','cj')
CJlow <- cast(CJlow,sub~coh,)
CJmed <- with(subset(Data2,traindiffcond=="average"),aggregate(cj,by=list(sub,coh),mean));names(CJmed) <- c('sub','coh','cj')
CJmed <- cast(CJmed,sub~coh)
CJhigh <- with(subset(Data2,traindiffcond=="easy"),aggregate(cj,by=list(sub,coh),mean));names(CJhigh) <- c('sub','coh','cj')
CJhigh <- cast(CJhigh,sub~coh)


#aggregate cj for model
x <- CJlow[,c(2:4)];xmed <- CJmed[,c(2:4)];xhigh <- CJhigh[,c(2:4)]
n <- length(x)

x <- x[,c("hard","average","easy")];
xmed <- xmed[,c("hard","average","easy")];
xhigh <- xhigh[,c("hard","average","easy")]

stripchart(x,ylim=c(3.75,5.5), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main=NULL, yaxt = 'n',family="A")
axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);
axis(2, seq(4,5.5,.5), cex.axis=1.5)
mtext("Confidence",2,at=4.75,line=2.5,cex=cex_axis);
means <- sapply(x, mean);n<- length(x)
for(i in seq(4,5.5,.5)) abline(h=i,col="lightgrey",lty = "dashed")
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="brown3",lwd=lwddat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N1),lwd=lwdgr,col="brown3")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="cyan4",lwd=lwddat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N1),lwd=lwdgr,col="cyan4")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="darkgoldenrod3",lwd=lwddat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N1),lwd=lwdgr,col="darkgoldenrod3")

# Plot confidence prediction ----------------------------------------------
par(mar=c(0,0,0,0))
plot.new()
text(.5,.75, labels="Model Prediction",cex = cex_legend+.5)
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
           main=NULL, yaxt = 'n',family="A")
mtext("Confidence",2,at=.75,line=2.5,cex=cex_axis);
mtext("Trial difficulty",1,3,at=1,cex=cex_axis)
axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);
axis(2, seq(.5,1,.1), cex.axis=1.5)
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
           main=NULL, yaxt = 'n',family="A")
mtext("Confidence",2,at=.75,line=2.5,cex=cex_axis);
axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);
mtext("Trial difficulty",1,3,at=1,cex=cex_axis)
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
axis(2, seq(.5,1,.1), cex.axis=1.5)

#'save
dev.off()
## Correlation Objective/Subjective drift ====
for (c in 1:Ncond) {
  obj_drift <- subset(df2,condition==cond_2[c])$Vo
  subj_drift <- subset(df2,condition==cond_2[c])$Vs
  drift_range <- c(min(min(subj_drift),min(obj_drift)),max(max(subj_drift),max(obj_drift)))
  plot(obj_drift~subj_drift,cex.axis=1.75,cex.lab=1.75,frame=F,pch=19, xlim=drift_range, ylim =drift_range );
  print(cor.test(obj_drift,subj_drift));
  abline(lm(obj_drift~subj_drift),lty=2)
  ;mtext(paste(cond_2[c],'r = ',round(cor(obj_drift,subj_drift),3)))
}
# Plot Vs ~ Condition -----------------------------------------------------
##Exp1
plot_drift <- with(df,aggregate(Vs,by=list(sub=sub,condition=condition),mean))
plot_drift <- cast(plot_drift,sub~condition)
plot_drift <- plot_drift[,c(3,4,2)] #Reorder columns to have hard -> easy
plot(colMeans(plot_drift),frame=F,type='n',cex.lab=2,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylab='',ylim=c(min(plot_drift),max(plot_drift)),
     xlab="Feedback",xaxt='n')
axis(1,1:Ncond,c("Negative","Average","Positive"),cex.axis=1.75)
mtext("Subjective drift",side = 2, line = 2.5, cex = 2)
for(i in 1:N1) lines(1:Ncond,plot_drift[i,1:Ncond],type='b',lty=2,col="grey",pch=19)
points(colMeans(plot_drift),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_drift),colSds(plot_drift,na.rm=T)/sqrt(N1),lwd=3)

##Exp2
plot_drift <- with(df2,aggregate(Vs,by=list(sub=sub,condition=condition),mean));
plot_drift <- cast(plot_drift,sub~condition)
plot_drift <- plot_drift[,c(4,2,3)] #Reorder columns to have hard -> easy
plot(colMeans(plot_drift),frame=F,type='n',cex.lab=2,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylim=c(min(plot_drift),max(plot_drift)),
     ylab="",xlab="Training",xaxt='n');
axis(1,1:Ncond,c("Hard","Average","Easy"),cex.axis=1.75)
mtext("Subjective drift",side = 2, line = 2.5, cex = 2)
for(i in 1:Nsub_2) lines(1:Ncond,plot_drift[i,1:Ncond],type='b',lty=2,col="grey",pch=19)
points(colMeans(plot_drift),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_drift),colSds(plot_drift,na.rm=T)/sqrt(Nsub_2),lwd=3)


# Plot DDM parameters test phase EXP1 ------------------------------------------
##Non-decision time
par(mfrow=c(1,2))
plot_ter <- with(param_1,aggregate(ter,by=list(sub=sub,condition=condition),mean))
plot_ter <- cast(plot_ter,sub~condition)
plot_ter <- plot_ter[,c(3,4,2)] #Reorder columns to have easy -> hard
plot(colMeans(plot_ter),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylim=c(min(plot_ter),max(plot_ter)),ylab="",
     xlab="Condition",xaxt='n',main="Non-decision time", cex.main = 2);
axis(1,1:Ncond,c("Negative","Average","Positive"),cex.axis=1.75)
for(i in 1:N1) lines(1:Ncond,plot_ter[i,1:Ncond],type='b',lty=2,col="grey",pch=19)
points(colMeans(plot_ter),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_ter),colSds(plot_ter,na.rm=T)/sqrt(N1),lwd=3,length=0)

##Bound
plot_bound <- with(param_1,aggregate(bound,by=list(sub=sub,condition=condition),mean))
plot_bound <- cast(plot_bound,sub~condition)
plot_bound <- plot_bound[,c(3,4,2)] #Reorder columns to have easy -> hard
plot(colMeans(plot_bound),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylim=c(min(plot_bound),max(plot_bound)),ylab="",
     xlab="Condition",xaxt='n',main="Bound", cex.main = 2);
axis(1,1:Ncond,c("Negative","Average","Positive"),cex.axis=1.75)
for(i in 1:N1) lines(1:Ncond,plot_bound[i,1:Ncond],type='b',lty=2,col="grey",pch=19)
points(colMeans(plot_bound),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_bound),colSds(plot_bound,na.rm=T)/sqrt(N1),lwd=3,length=0)

##Drift interaction
par(mfrow=c(1,1))
plot_drift_minus <- with(subset(param_1,condition=="lowSC"),
                         aggregate(drift,by=list(sub=sub,difflevel=difflevel),mean))
plot_drift_minus <- cast(plot_drift_minus,sub~difflevel)
plot_drift_minus <- plot_drift_minus[,c(4,2,3)] #Reorder columns to have easy -> hard
plot(colMeans(plot_drift_minus),frame=F,type='n',cex.lab=2,cex.axis=1.75,xlim=c(.8,Ncond+.2),
     ylim=c(min(plot_drift_minus),.25),ylab="Drift rate",xlab="Trial Difficulty",xaxt='n');
axis(1,1:Ncond,c("Hard","Medium","Easy"),cex.axis=1.75)
points(colMeans(plot_drift_minus),type='b',lwd=5,col="brown3")
error.bar(1:Ncond,colMeans(plot_drift_minus),
          colSds(plot_drift_minus,na.rm=T)/sqrt(N1),lwd=3,length=0,col="brown3")

plot_drift_control <- with(subset(param_1,condition=="mediumSC"),
                           aggregate(drift,by=list(sub=sub,difflevel=difflevel),mean))
plot_drift_control <- cast(plot_drift_control,sub~difflevel)
plot_drift_control <- plot_drift_control[,c(4,2,3)] #Reorder columns to have easy -> hard
points(colMeans(plot_drift_control),type='b',lwd=5,col="cyan4")
error.bar(1:Ncond,colMeans(plot_drift_control),
          colSds(plot_drift_control,na.rm=T)/sqrt(N1),lwd=3,length=0,col="cyan4")

plot_drift_plus <- with(subset(param_1,condition=="highSC"),
                        aggregate(drift,by=list(sub=sub,difflevel=difflevel),mean))
plot_drift_plus <- cast(plot_drift_plus,sub~difflevel)
plot_drift_plus <- plot_drift_plus[,c(4,2,3)] #Reorder columns to have easy -> hard
points(colMeans(plot_drift_plus),type='b',lwd=5,col="darkgoldenrod3")
error.bar(1:Ncond,colMeans(plot_drift_plus),
          colSds(plot_drift_plus,na.rm=T)/sqrt(N1),lwd=3,length=0,col="darkgoldenrod3")
legend("topleft",border=F,legend=c("Negative","Average","Positive"),lwd=2,
       col=c("brown3","cyan4","darkgoldenrod3"),bty="n",cex=1.5,title = "Condition")


# Plot DDM parameters test phase EXP2 -------------------------------------
##Non-decision time
par(mfrow=c(1,2))
plot_ter <- with(param2,aggregate(ter,by=list(sub=sub,condition=condition),mean))
plot_ter <- cast(plot_ter,sub~condition)
plot_ter <- plot_ter[,c(4,2,3)] #Reorder columns to have easy -> hard
plot(colMeans(plot_ter),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylim=c(min(plot_ter),max(plot_ter)),ylab="",
     xlab="Condition",xaxt='n',main="Non-decision time", cex.main = 2);
axis(1,1:Ncond,c("Hard","Medium","Easy"),cex.axis=1.75)
for(i in 1:Nsub_2) lines(1:Ncond,plot_ter[i,1:Ncond],type='b',lty=2,col="grey",pch=19)
points(colMeans(plot_ter),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_ter),colSds(plot_ter,na.rm=T)/sqrt(Nsub_2),lwd=3,length=0)

##Bound
plot_bound <- with(param2,aggregate(bound,by=list(sub=sub,condition=condition),mean))
plot_bound <- cast(plot_bound,sub~condition)
plot_bound <- plot_bound[,c(4,2,3)] #Reorder columns to have easy -> hard
plot(colMeans(plot_bound),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylim=c(min(plot_bound),max(plot_bound)),ylab="",
     xlab="Condition",xaxt='n',main="Bound", cex.main = 2);
axis(1,1:Ncond,c("Hard","Medium","Easy"),cex.axis=1.75)
for(i in 1:Nsub_2) lines(1:Ncond,plot_bound[i,1:Ncond],type='b',lty=2,col="grey",pch=19)
points(colMeans(plot_bound),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_bound),colSds(plot_bound,na.rm=T)/sqrt(Nsub_2),lwd=3,length=0)

##Drift interaction
par(mfrow=c(1,1))
plot_drift_minus <- with(subset(param2,condition=="hard"),
                         aggregate(drift,by=list(sub=sub,difflevel=difflevel),mean))
plot_drift_minus <- cast(plot_drift_minus,sub~difflevel)
plot_drift_minus <- plot_drift_minus[,c(4,2,3)] #Reorder columns to have easy -> hard
plot(colMeans(plot_drift_minus),frame=F,type='n',cex.lab=2,cex.axis=1.75,xlim=c(.8,Ncond+.2),
     ylim=c(min(plot_drift_minus),.25),ylab="Drift rate",xlab="Trial Difficulty",xaxt='n');
axis(1,1:Ncond,c("Hard","Medium","Easy"),cex.axis=1.75)
points(colMeans(plot_drift_minus),type='b',lwd=5,col="brown3")
error.bar(1:Ncond,colMeans(plot_drift_minus),
          colSds(plot_drift_minus,na.rm=T)/sqrt(Nsub_2),lwd=3,length=0,col="brown3")

plot_drift_control <- with(subset(param2,condition=="average"),
                           aggregate(drift,by=list(sub=sub,difflevel=difflevel),mean))
plot_drift_control <- cast(plot_drift_control,sub~difflevel)
plot_drift_control <- plot_drift_control[,c(4,2,3)] #Reorder columns to have easy -> hard
points(colMeans(plot_drift_control),type='b',lwd=5,col="cyan4")
error.bar(1:Ncond,colMeans(plot_drift_control),
          colSds(plot_drift_control,na.rm=T)/sqrt(Nsub_2),lwd=3,length=0,col="cyan4")

plot_drift_plus <- with(subset(param2,condition=="easy"),
                        aggregate(drift,by=list(sub=sub,difflevel=difflevel),mean))
plot_drift_plus <- cast(plot_drift_plus,sub~difflevel)
plot_drift_plus <- plot_drift_plus[,c(4,2,3)] #Reorder columns to have easy -> hard
points(colMeans(plot_drift_plus),type='b',lwd=5,col="darkgoldenrod3")
error.bar(1:Ncond,colMeans(plot_drift_plus),
          colSds(plot_drift_plus,na.rm=T)/sqrt(Nsub_2),lwd=3,length=0,col="darkgoldenrod3")
legend("topleft",border=F,legend=c("Hard","Medium","Easy"),lwd=2,
       col=c("brown3","cyan4","darkgoldenrod3"),bty="n",cex=1.5,title = "Condition")



# Standardizing -- DEAD END -----------------------------------------------

cj2_pred <- with(Simuls2,aggregate(cj,by=list(sub=sub,condition=condition,coh=coh),mean))
cj2_pred$x <- scale(cj2_pred$x)

cj1_pred <- with(Simuls,aggregate(cj,by=list(sub=sub,condition=condition,coh=coh),mean))
cj1_pred$x <- scale(cj1_pred$x)

cj1_obs <- with(Data1,aggregate(cj,by=list(sub=sub,condition=selfconf,coh=coh),mean))
cj1_obs$x <- scale(cj1_obs$x)

cj2_obs <- with(Data2,aggregate(cj,by=list(sub=sub,condition=traindiffcond,coh=coh),mean))
cj2_obs$x <- scale(cj2_obs$x)

# Plot confidence prediction ----------------------------------------------
Simuls$cj <- Simuls$cj_cont
## Experiment 1
#Aggregate conf for data
CJlow <- subset(cj1_obs, condition =="lowSC");
CJlow <- CJlow[,c("sub","coh","x")];names(CJlow) <- c('sub','coh','cj')
CJlow <- cast(CJlow,sub~coh,)
CJmed <- subset(cj1_obs, condition =="mediumSC");
CJmed <- CJmed[,c("sub","coh","x")];names(CJmed) <- c('sub','coh','cj')
CJmed <- cast(CJmed,sub~coh,)
CJhigh <- subset(cj1_obs, condition =="highSC");
CJhigh <- CJhigh[,c("sub","coh","x")];names(CJhigh) <- c('sub','coh','cj')
CJhigh <- cast(CJhigh,sub~coh,)

#aggregate CJ for model
x <- CJlow[,c(2:4)];xmed <- CJmed[,c(2:4)];xhigh <- CJhigh[,c(2:4)]
n <- length(x)

x <- x[,c("hard","average","easy")];
xmed <- xmed[,c("hard","average","easy")];
xhigh <- xhigh[,c("hard","average","easy")]

xlow_sim = subset(cj1_pred, condition =="lowSC");
xlow_sim <- xlow_sim[,c("sub","coh","x")];names(xlow_sim) <- c('sub','coh','cj')
xlow_sim <- cast(xlow_sim,sub~coh,)
xmed_sim = subset(cj1_pred, condition =="mediumSC");
xmed_sim <- xmed_sim[,c("sub","coh","x")];names(xmed_sim) <- c('sub','coh','cj')
xmed_sim <- cast(xmed_sim,sub~coh,)
xhigh_sim = subset(cj1_pred, condition =="highSC");
xhigh_sim <- xhigh_sim[,c("sub","coh","x")];names(xhigh_sim) <- c('sub','coh','cj')
xhigh_sim <- cast(xhigh_sim,sub~coh,)
xlow_sim <- xlow_sim[,c("hard","average","easy")];
xmed_sim <- xmed_sim[,c("hard","average","easy")];
xhigh_sim <- xhigh_sim[,c("hard","average","easy")];

par(mfrow=c(1,1))
stripchart(x, ylim=c(-1.5,1.5), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main="",cex.axis=1.25)
mtext("Z Confidence",2,at=0,line=2,cex=1.75);axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);mtext("Trial difficulty",1,3,at=1,cex=1.75)
means <- sapply(x, mean);n<- length(x)
# lines(0:(n-1),colMeans(xlow_sim,na.rm=T),type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=16,col=rgb(1,0,0,.5))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xlow_sim,na.rm=T) + (colSds(as.matrix(xlow_sim))/sqrt(N1)),
                             (colMeans(xlow_sim,na.rm=T) - colSds(as.matrix(xlow_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(1,0,0,.2))
# lines(0:(n-1),colMeans(xmed_sim,na.rm=T),type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=16,col=rgb(1,.5,0,.5))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xmed_sim,na.rm=T) + (colSds(as.matrix(xmed_sim),na.rm=T)/sqrt(N1)),
                             (colMeans(xmed_sim,na.rm=T) - colSds(as.matrix(xmed_sim),na.rm=T)/sqrt(N1))[3:1]),
        border=F,col=rgb(1,.5,0,.2))
# lines(0:(n-1),colMeans(xhigh_sim,na.rm=T),type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=16,col=rgb(0,0,1,.5))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xhigh_sim,na.rm=T) + (colSds(as.matrix(xhigh_sim),na.rm=T)/sqrt(N1)),
                             (colMeans(xhigh_sim,na.rm=T) - colSds(as.matrix(xhigh_sim),na.rm=T)/sqrt(N1))[3:1]),
        border=F,col=rgb(0,0,1,.2))
legend(0,1.75,legend=c("Negative","Average","Positive"),
       title = "Feedback condition",pch=rep(16,3),bty = "n",inset=.1, 
       cex = 1.25,col=c("brown3","cyan4","darkgoldenrod3"))
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
#Aggregate conf for data
CJlow <- subset(cj2_obs, condition =="hard");
CJlow <- CJlow[,c("sub","coh","x")];names(CJlow) <- c('sub','coh','cj')
CJlow <- cast(CJlow,sub~coh,)
CJmed <- subset(cj2_obs, condition =="average");
CJmed <- CJmed[,c("sub","coh","x")];names(CJmed) <- c('sub','coh','cj')
CJmed <- cast(CJmed,sub~coh,)
CJhigh <- subset(cj2_obs, condition =="easy");
CJhigh <- CJhigh[,c("sub","coh","x")];names(CJhigh) <- c('sub','coh','cj')
CJhigh <- cast(CJhigh,sub~coh,)

#aggregate CJ for model
x <- CJlow[,c(2:4)];xmed <- CJmed[,c(2:4)];xhigh <- CJhigh[,c(2:4)]
n <- length(x)

x <- x[,c("hard","average","easy")];
xmed <- xmed[,c("hard","average","easy")];
xhigh <- xhigh[,c("hard","average","easy")]

xlow_sim = subset(cj2_pred, condition =="hard");
xlow_sim <- xlow_sim[,c("sub","coh","x")];names(xlow_sim) <- c('sub','coh','cj')
xlow_sim <- cast(xlow_sim,sub~coh,)
xmed_sim = subset(cj2_pred, condition =="average");
xmed_sim <- xmed_sim[,c("sub","coh","x")];names(xmed_sim) <- c('sub','coh','cj')
xmed_sim <- cast(xmed_sim,sub~coh,)
xhigh_sim = subset(cj2_pred, condition =="easy");
xhigh_sim <- xhigh_sim[,c("sub","coh","x")];names(xhigh_sim) <- c('sub','coh','cj')
xhigh_sim <- cast(xhigh_sim,sub~coh,)
xlow_sim <- xlow_sim[,c("hard","average","easy")];
xmed_sim <- xmed_sim[,c("hard","average","easy")];
xhigh_sim <- xhigh_sim[,c("hard","average","easy")];

stripchart(x, ylim=c(-1.5,1.5), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main="",cex.axis=1.25)
mtext("Z Confidence",2,at=0,line=2,cex=1.75);axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);mtext("Trial difficulty",1,3,at=1,cex=1.75)
means <- sapply(x, mean);n<- length(x)
# lines(0:(n-1),colMeans(xlow_sim,na.rm=T),type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=16,col=rgb(1,0,0,.5))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xlow_sim,na.rm=T) + (colSds(as.matrix(xlow_sim))/sqrt(Nsub_2)),(colMeans(xlow_sim,na.rm=T) - colSds(as.matrix(xlow_sim))/sqrt(Nsub_2))[3:1]),
        border=F,col=rgb(1,0,0,.2))
# lines(0:(n-1),colMeans(xmed_sim,na.rm=T),type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=16,col=rgb(1,.5,0,.5))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xmed_sim,na.rm=T) + (colSds(as.matrix(xmed_sim),na.rm=T)/sqrt(Nsub_2)),
                             (colMeans(xmed_sim,na.rm=T) - colSds(as.matrix(xmed_sim),na.rm=T)/sqrt(Nsub_2))[3:1]),
        border=F,col=rgb(1,.5,0,.2))
# lines(0:(n-1),colMeans(xhigh_sim,na.rm=T),type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=16,col=rgb(0,0,1,.5))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xhigh_sim,na.rm=T) + (colSds(as.matrix(xhigh_sim),na.rm=T)/sqrt(Nsub_2)),(colMeans(xhigh_sim,na.rm=T) - colSds(as.matrix(xhigh_sim),na.rm=T)/sqrt(Nsub_2))[3:1]),
        border=F,col=rgb(0,0,1,.2))
legend(0,1.75,legend=c("Hard","Average","Easy"),
       title = "Training condition",pch=rep(16,3),bty = "n",inset=.1,
       cex = 1.25,col=c("brown3","cyan4","darkgoldenrod3"))
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="brown3",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(Nsub_2),lwd=lwdgr,col="brown3")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="cyan4",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(Nsub_2),lwd=lwdgr,col="cyan4")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="darkgoldenrod3",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(Nsub_2),lwd=lwdgr,col="darkgoldenrod3")

