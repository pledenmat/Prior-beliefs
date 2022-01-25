error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
cexkl <- 2;cexgr <- 2;lwdgr <- 3; lwddat <- 2

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

par(mfrow=c(1,1))
stripchart(x, ylim=c(0.6,1), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main="",cex.axis=1.25)
mtext("Accuracy",2,at=.8,line=3,cex=1.75);axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);mtext("Trial difficulty",1,3,at=1,cex=1.75)
means <- sapply(x, mean);n<- length(x)
# lines(0:(n-1),colMeans(x_sim,na.rm=T),type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=16,col=rgb(1,0,0,.5))
polygon(c(0:(n-1),(n-1):0),c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N)),(colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N))[3:1]),
        border=F,col=rgb(1,0,0,.2))
# lines(0:(n-1),colMeans(xmed_sim,na.rm=T),type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=16,col=rgb(1,.5,0,.5))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xmed_sim,na.rm=T) + (colSds(as.matrix(xmed_sim),na.rm=T)/sqrt(N)),(colMeans(xmed_sim,na.rm=T) - colSds(as.matrix(xmed_sim),na.rm=T)/sqrt(N))[3:1]),
        border=F,col=rgb(1,.5,0,.2))
# lines(0:(n-1),colMeans(xhigh_sim,na.rm=T),type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=16,col=rgb(0,0,1,.5))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xhigh_sim,na.rm=T) + (colSds(as.matrix(xhigh_sim),na.rm=T)/sqrt(N)),(colMeans(xhigh_sim,na.rm=T) - colSds(as.matrix(xhigh_sim),na.rm=T)/sqrt(N))[3:1]),
        border=F,col=rgb(0,0,1,.2))
legend(.05,1,legend=c("Negative","Average","Positive"),title = "Fake Feedback condition",pch=rep(16,3),bty = "n",inset=.1, cex = 1.25,col=c("red","orange","blue"))
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="red",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="red")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="orange",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="orange")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="blue",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="blue")

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

stripchart(x, ylim=c(.6,1), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main="",cex.axis=1.25)
mtext("Accuracy",2,at=.8,line=3,cex=1.75);axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);mtext("Trial difficulty",1,3,at=1,cex=1.75)
means <- sapply(x, mean);n<- length(x)
# lines(0:(n-1),colMeans(x_sim,na.rm=T),type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=16,col=rgb(1,0,0,.5))
polygon(c(0:(n-1),(n-1):0),c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N)),(colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N))[3:1]),
        border=F,col=rgb(1,0,0,.2))
# lines(0:(n-1),colMeans(xmed_sim,na.rm=T),type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=16,col=rgb(1,.5,0,.5))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xmed_sim,na.rm=T) + (colSds(as.matrix(xmed_sim),na.rm=T)/sqrt(N)),(colMeans(xmed_sim,na.rm=T) - colSds(as.matrix(xmed_sim),na.rm=T)/sqrt(N))[3:1]),
        border=F,col=rgb(1,.5,0,.2))
# lines(0:(n-1),colMeans(xhigh_sim,na.rm=T),type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=16,col=rgb(0,0,1,.5))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xhigh_sim,na.rm=T) + (colSds(as.matrix(xhigh_sim),na.rm=T)/sqrt(N)),(colMeans(xhigh_sim,na.rm=T) - colSds(as.matrix(xhigh_sim),na.rm=T)/sqrt(N))[3:1]),
        border=F,col=rgb(0,0,1,.2))
legend(.05,1,legend=c("Hard","Average","Easy"),title = "Training difficulty condition",pch=rep(16,3),bty = "n",inset=.1, cex = 1.5,col=c("red","orange","blue"))
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="red",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="red")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="orange",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="orange")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="blue",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="blue")