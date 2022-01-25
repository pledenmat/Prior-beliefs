setwd(paste0(curdir,"/realdata_fit/"))
for(i in 1:50){ #123=tinne
  if(i == 1){
    Data <- read.csv(paste0('RealData_1B/selfconfidence1B_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
  }else{
    temp <- read.csv(paste0('RealData_1B/selfconfidence1B_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
    Data <- rbind(Data,temp)
  }
}


Training <- subset(Data, running == "training")
Data['response'] <- 0
Data$response[Data$resp == "['n']"] <- 1
Data <- subset(Data,running == "main")

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

exclusion <- c(12) #Hélène's list
Data <- subset(Data,!(sub %in% exclusion)) 

Data$response[Data$response==0] <- -1

df <- Data[,c("sub","task","traindiffcond","trialdifflevel","rt","response","cor","cj","RTconf")]
names(df) <- c("sub","task","traindiffcond","coh","rt","resp","cor","cj","RTconf")
df$rt <- df$rt/1000;df$RTconf <- df$RTconf/1000
# write.csv(df,"dataexp2_helene_full.csv",row.names = FALSE)

Training <- subset(Training,!(sub %in% exclusion))
Training['response'] <- -1
Training$response[Training$resp == "['n']"] <- 1
df_train <- Training[,c("sub","task","traindiffcond","trialdifflevel","rt","response","cor","cj","RTconf")]
names(df_train) <- c("sub","task","traindiffcond","coh","rt","resp","cor","cj","RTconf")
df_train$rt <- df_train$rt/1000;df_train$RTconf <- df_train$RTconf/1000
setwd(curdir)
# write.csv(df_train,"dataexp2_helene_training_full.csv",row.names = FALSE)
