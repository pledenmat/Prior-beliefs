##' Preprocessing of both experiments
##' 
##' - Loads raw data files
##' - Aggregates data into one dataframe per experiment
##' - Trims 200 ms < RT < 5000 ms
##' - Removes participants that did not exceed chance level performance in the test phase
##' - Writes csv files of aggregated data for each experiment
##' 
##' Note : exp1 = fake feedback, exp2 = training difficulty

rm(list=ls())
library(reshape)
library(effects)
library(lmerTest)
library(scales)
library(prob)
library(myPackage)
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir) ## change our current working directory

plot <- F

# Experiment 1 ------------------------------------------------------------
N <- 50
go_to("data")
for(i in 1:N){
  if(i == 1){
    Data_exp1 <- read.csv(paste0('RealData_1A/selfconfidence1A_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
  }else{
    temp <- read.csv(paste0('RealData_1A/selfconfidence1A_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
    Data_exp1 <- rbind(Data_exp1,temp)
  }
}

head(Data_exp1)

Data_exp1 <- subset(Data_exp1, running == "main")

Data_exp1['response'] <- 0
Data_exp1$response[Data_exp1$resp == "['n']"] <- 1

Data_exp1 <- subset(Data_exp1,rt>200) #There was no trial below 200ms
Data_exp1 <- subset(Data_exp1,rt<5000)

## Diagnostic plot per participant and task + chance performance testing
N <- length(unique(Data_exp1$sub))
sub_list <- unique(Data_exp1$sub)
exclusion <- c()
tasks <- unique(Data_exp1$task)
par(mfrow=c(2,2))
for(i in 1:N){
  for (t in tasks) {
    tempDat <- subset(Data_exp1, sub == sub_list[i] & task == t)
    acc_block <- with(tempDat,aggregate(cor,by=list(block=block),mean))
    bias_block <- with(tempDat,aggregate(response,by=list(block=block),mean))
    test <- binom.test(length(tempDat$cor[tempDat$cor==1]),n=length(tempDat$cor),alternative = "greater")
    print(paste("In t0, sub",sub_list[i],"p =", round(test$p.value,3),"compared to chance"))
    if (plot) {
      plot(acc_block,ylab="Acc (.) and bias (x)",frame=F,ylim=c(0,1));abline(h=.5,lty=2,col="grey")
      points(bias_block,pch=4)
      plot(tempDat$rt/1000,frame=F,col=c("black"),main=paste('subject',i,"task :",t),ylab="RT")
      plot(tempDat$cj,frame=F,col=c("black"),ylim=c(1,6),ylab="conf")
      plot(tempDat$RTconf,frame=F,col=c("black"),ylab="RT_conf")
      
    }
    if (test$p.value > .05) {
      exclusion <- c(exclusion,sub_list[i])
    }
  }
}

Data_exp1 <- subset(Data_exp1,!(sub %in% exclusion)) #Sub 10 and 49 removed

Data_exp1$response[Data_exp1$response==0] <- -1

df_1 <- Data_exp1[,c("sub","task","selfconf","difflevel","rt","response","cor","cj","RTconf","block")]
names(df_1) <- c("sub","task","selfconf","coh","rt","resp","cor","cj","RTconf","block")
## Convert into seconds
df_1$rt <- df_1$rt/1000
df_1$RTconf <- df_1$RTconf/1000
# Experiment 2 ------------------------------------------------------------
go_to("data")
N <- 50
for(i in 1:N){ 
  if(i == 1){ 
    Data_exp2 <- read.csv(paste0('RealData_1B/selfconfidence1B_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
  }else{
    temp <- read.csv(paste0('RealData_1B/selfconfidence1B_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
    Data_exp2 <- rbind(Data_exp2,temp)
  }
}

head(Data_exp2)

Data_exp2 <- subset(Data_exp2, running == "main")

Data_exp2['response'] <- 0
Data_exp2$response[Data_exp2$resp == "['n']"] <- 1

Data_exp2 <- subset(Data_exp2,rt>200) #There was no trial below 200ms
Data_exp2 <- subset(Data_exp2,rt<5000)

## Diagnostic plot per participant and task + chance performance testing
N <- length(unique(Data_exp2$sub)) 
sub_list <- unique(Data_exp2$sub)
exclusion <- c()
tasks <- unique(Data_exp2$task)
par(mfrow=c(2,2))
for(i in 1:N){
  for (t in tasks) {
    tempDat <- subset(Data_exp2,sub==sub_list[i]&task==t)
    acc_block <- with(tempDat,aggregate(cor,by=list(block=block),mean))
    bias_block <- with(tempDat,aggregate(response,by=list(block=block),mean))
    test <- binom.test(length(tempDat$cor[tempDat$cor==1]),n=length(tempDat$cor),alternative = "greater")
    print(paste("In t0, sub",sub_list[i],"p =", round(test$p.value,3),"compared to chance"))
    if (test$p.value > .05) {
      exclusion <- c(exclusion,sub_list[i])
    }
    if (plot) {
      plot(acc_block,ylab="Acc (.) and bias (x)",frame=F,ylim=c(0,1));abline(h=.5,lty=2,col="grey")
      points(bias_block,pch=4)
      plot(tempDat$rt/1000,frame=F,col=c("black"),main=paste('subject',i,"task :",t),ylab="RT",ylim=c(0,5))
      plot(tempDat$cj,frame=F,col=c("black"),ylim=c(1,6),ylab="conf")
      plot(tempDat$RTconf,frame=F,col=c("black"),ylab="RT_conf")
    }
  }
}

Data_exp2 <- subset(Data_exp2,!(sub %in% exclusion)) #Sub 12, 32 and 50 removed

Data_exp2$response[Data_exp2$response==0] <- -1

df_2 <- Data_exp2[,c("sub","task","traindiffcond","trialdifflevel","rt","response","cor","cj","RTconf","block")]
names(df_2) <- c("sub","task","traindiffcond","coh","rt","resp","cor","cj","RTconf","block")
## Convert into seconds
df_2$rt <- df_2$rt/1000
df_2$RTconf <- df_2$RTconf/1000


# Export aggregated data in csv files -------------------------------------
go_to("results")
write.csv(df_1,"data_exp1.csv",row.names = FALSE)
write.csv(df_2,"data_exp2.csv",row.names = FALSE)

