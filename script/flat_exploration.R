for (i in 1:Nsub_2) {

    tempmean <- as.numeric(means_fullconfRT[Nsub_2*(2-1)+i,])
    tempmean <- tempmean[complete.cases(tempmean)]
    smoothed <- lowess(tempmean[2:501],f=.05)
    plot(tempmean[2:501],main=paste(subs_2[i],"easy","mean"),
         xlab="Vs",ylab="Mean over repetitions",xaxt='n')
    lines(smoothed,col="green",lwd=2)
    axis(1,at=seq(0,500,100),labels = seq(0,.5,.1))
    abline(v=which.min(tempmean[2:501]),col="red")
    abline(v=which.min(smoothed$y),col="green")
}

v_train_easy <- v_train[,2]
bound_train_easy <- bound_train[,2]

weird <- c(1:6,9,15,18,27,31,33,34,40,41)
flat <- c(1:6,9,13,14,18,20,23,27,31,40,41)
highacc <- subset(acc_train,acc>=.95)

weird_ind <- which(subs_2 %in% weird)
flat_ind <- which(subs_2 %in% flat)

acc_weird <- subset(acc_train,sub %in% weird)
v_weird <- v_train_easy[weird_ind]
bound_weird <- bound_train_easy[weird_ind]

acc_flat <- subset(acc_train,sub %in% flat)
v_flat <- v_train_easy[flat_ind]
bound_flat <- bound_train_easy[flat_ind]
