rm(list=ls())
library(prob)
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir) ## change our current working directory

# Experiment 1 ------------------------------------------------------------
for(i in 1:50){
  if(i == 1){
    Data <- read.csv(paste0('RealData_1A/selfconfidence1A_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
  }else{
    temp <- read.csv(paste0('RealData_1A/selfconfidence1A_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
    Data <- rbind(Data,temp)
  }
}

head(Data)

Data <- subset(Data, running == "main")

Data['response'] <- 0
Data$response[Data$resp == "['n']"] <- 1

Data <- subset(Data,rt>.2) #There was no trial below 200ms

## Diagnostic plot per participant and task + chance performance testing
N <- length(unique(Data$sub)); subs <- unique(Data$sub); exclusion <- c()
tasks <- unique(Data$task)
par(mfrow=c(2,2))
for(i in 1:N){
  for (t in tasks) {
    tempDat <- subset(Data,sub==subs[i]&task==t)
    acc_block <- with(tempDat,aggregate(cor,by=list(block=block),mean))
    bias_block <- with(tempDat,aggregate(response,by=list(block=block),mean))
    test <- binom.test(length(tempDat$cor[tempDat$cor==1]),n=length(tempDat$cor),alternative = "greater")
    print(paste("In t0, sub",subs[i],"p =", round(test$p.value,3),"compared to chance"))
    if (test$p.value > .05) {
      exclusion <- c(exclusion,subs[i])
      plot(acc_block,ylab="Acc (.) and bias (x)",frame=F,ylim=c(0,1));abline(h=.5,lty=2,col="grey")
      points(bias_block,pch=4)
      plot(tempDat$rt/1000,frame=F,col=c("black"),main=paste('subject',i,"task :",t),ylab="RT",ylim=c(0,5))
      plot(tempDat$cj,frame=F,col=c("black"),ylim=c(1,6),ylab="conf")
      plot(tempDat$RTconf,frame=F,col=c("black"),ylab="RT_conf")
    }
  }
}

Data <- subset(Data,!(sub %in% exclusion)) #Sub 10 and 49 removed
# Experiment 2 ------------------------------------------------------------

for(i in 1:50){ #123=tinne
  if(i == 1){ 
    Data <- read.csv(paste0('RealData_1B/selfconfidence1B_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
  }else{
    temp <- read.csv(paste0('RealData_1B/selfconfidence1B_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
    Data <- rbind(Data,temp)
  }
}

head(Data)

Data <- subset(Data, running == "main")

Data['response'] <- 0
Data$response[Data$resp == "['n']"] <- 1

Data <- subset(Data,rt>200) #There was no trial below 200ms

## Diagnostic plot per participant and task + chance performance testing
N <- length(unique(Data$sub)); subs <- unique(Data$sub); exclusion <- c()
tasks <- unique(Data$task)
par(mfrow=c(2,2))
# at_chance <- matrix(NA,nrow=N,ncol=2)
for(i in 1:N){
  for (t in tasks) {
    tempDat <- subset(Data,sub==subs[i]&task==t)
    acc_block <- with(tempDat,aggregate(cor,by=list(block=block),mean))
    # at_chance[i,] <- c(subs[i], t.test(acc_block$x,mu=.5)$p.value)
    bias_block <- with(tempDat,aggregate(response,by=list(block=block),mean))
    test <- binom.test(length(tempDat$cor[tempDat$cor==1]),n=length(tempDat$cor),alternative = "greater")
    print(paste("In t0, sub",subs[i],"p =", round(test$p.value,3),"compared to chance"))
    if (test$p.value > .05) {
      exclusion <- c(exclusion,subs[i])
      plot(acc_block,ylab="Acc (.) and bias (x)",frame=F,ylim=c(0,1));abline(h=.5,lty=2,col="grey")
      points(bias_block,pch=4)
      plot(tempDat$rt/1000,frame=F,col=c("black"),main=paste('subject',i,"task :",t),ylab="RT",ylim=c(0,5))
      plot(tempDat$cj,frame=F,col=c("black"),ylim=c(1,6),ylab="conf")
      plot(tempDat$RTconf,frame=F,col=c("black"),ylab="RT_conf")
    }
  }
}

Data <- subset(Data,!(sub %in% exclusion)) #Sub 12, 32 and 50 removed
