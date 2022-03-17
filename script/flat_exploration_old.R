sds_fullconfRT <- with(cost_df_fb,aggregate(cost,by=list(Vs=Vs,sub=sub,traindiffcond=traindiffcond),sd))
sds_fullconfRT <- cast(sds_fullconfRT,traindiffcond+sub~Vs)
for (i in 1:Nsub_2) {
    tempmean <- as.numeric(means_fullconfRT[Nsub_2*(2-1)+i,])
    tempmean <- tempmean[complete.cases(tempmean)]
    plot(as.numeric(tempmean[2:501]),main=paste(subs_2[i],"mean"))
    abline(v=which.min(tempmean))
}

test <- subset(Data2,sub==45&traindiffcond=="easy")

total <- c(highacc,c(37,15,34,33,24))
ind <- which(subs_2 %in% total)
v_train_easy[ind]
ind_ok <- which(!(subs_2 %in% total))
v_train_easy[ind_ok]

for (i in ind_ok) {
  tempmean <- as.numeric(means_fullconfRT[Nsub_2*(2-1)+i,])
  tempmean <- tempmean[complete.cases(tempmean)]
  plot(as.numeric(tempmean[2:501]),main=paste(subs_2[i],v_train_easy[i]))
  abline(v=which.min(tempmean))
}
for (i in ind) {
  tempmean <- as.numeric(means_fullconfRT[Nsub_2*(2-1)+i,])
  tempmean <- tempmean[complete.cases(tempmean)]
  plot(as.numeric(tempmean[2:501]),main=paste(subs_2[i],v_train_easy[i]))
  abline(v=which.min(tempmean[2:501]))
}

truc <- subset(Data2_train,sub==subs_2[10]&traindiffcond=="easy")
acc_train <- with(subset(Data2_train,traindiffcond=="easy"),aggregate(cor,by=list(sub),mean))
names(acc_train) <- c("sub","acc")
highacc <- acc_train[acc_train$acc>=.95]

bound_train[ind,2]
v_train_easy[ind]
v_train_easy[ind]/bound_train[ind,2]
v_train_easy[ind_ok]/bound_train[ind_ok,2]

maxs <- with(subset(cost_df_fb,traindiffcond=="easy"),aggregate(cost,by=list(sub=sub,Vs=Vs),max))
for (i in 1:Nsub_2) {
  tempdat <- subset(maxs,sub==subs_2[i])
  plot(tempdat$x,main=subs_2[i])
}
