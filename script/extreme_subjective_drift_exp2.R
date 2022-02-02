## Experiment 2 ====
plot_drift <- with(df2_test,aggregate(Vs,by=list(sub=sub,condition=condition),mean));
plot_drift <- cast(plot_drift,sub~condition)
plot_drift <- plot_drift[complete.cases(plot_drift),]
subs_ok <- plot_drift$sub
plot_drift <- plot_drift[,c(4,2,3)] #Reorder columns to have hard -> easy
plot(colMeans(plot_drift),frame=F,type='n',cex.lab=2,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylim=c(min(plot_drift),max(plot_drift)),
     ylab="",xlab="Training Difficulty",xaxt='n');
axis(1,1:Ncond,c("Hard","Average","Easy"),cex.axis=1.75)
mtext("Subjective drift",side = 2, line = 2.5, cex = 2)
for(i in 1:Nsub_2) lines(1:Ncond,plot_drift[i,1:Ncond],type='b',lty=2,col="grey",pch=19)
points(colMeans(plot_drift),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_drift),colSds(plot_drift,na.rm=T)/sqrt(Nsub_2),lwd=3,length=.05)


Simul_test <- subset(Simuls2,sub %in% subs_ok)

Simul_test$cj <- Simul_test$cj_cont
#Aggregate conf for Simul_test
CJlow <- with(subset(Simul_test,condition=="hard"),aggregate(cj,by=list(sub,coh),mean));names(CJlow) <- c('sub','coh','cj')
CJlow <- cast(CJlow,sub~coh,)
CJmed <- with(subset(Simul_test,condition=="average"),aggregate(cj,by=list(sub,coh),mean));names(CJmed) <- c('sub','coh','cj')
CJmed <- cast(CJmed,sub~coh)
CJhigh <- with(subset(Simul_test,condition=="easy"),aggregate(cj,by=list(sub,coh),mean));names(CJhigh) <- c('sub','coh','cj')
CJhigh <- cast(CJhigh,sub~coh)


#aggregate cj for model
x <- CJlow[,c(2:4)];xmed <- CJmed[,c(2:4)];xhigh <- CJhigh[,c(2:4)]
n <- length(x)

x <- x[,c("hard","average","easy")];
xmed <- xmed[,c("hard","average","easy")];
xhigh <- xhigh[,c("hard","average","easy")]

stripchart(x, ylim=c(.45,.95), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main="Model prediction",cex.axis=1.25)
mtext("Confidence",2,at=.7,line=2.5,cex=1.75);axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);mtext("Trial difficulty",1,3,at=1,cex=1.75)
means <- sapply(x, mean);n<- length(x)
legend(0,.95,legend=c("Hard","Average","Easy"),title = "Training difficulty condition",pch=rep(16,3),bty = "n",inset=.1, cex = 1.25,col=c("red","orange","blue"))
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="red",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(Nsub_2),lwd=lwdgr,length=.05,col="red")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="orange",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(Nsub_2),lwd=lwdgr,length=.05,col="orange")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="blue",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(Nsub_2),lwd=lwdgr,length=.05,col="blue")
