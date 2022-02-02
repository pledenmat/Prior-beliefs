error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
cexkl <- 2;cexgr <- 2;lwdgr <- 3; lwddat <- 2

## Experiment 1
#Aggregate conf for data
rtlow <- with(subset(Simuls,condition=="lowSC"),aggregate(rt,by=list(sub,coh),mean));names(rtlow) <- c('sub','coh','rt')
rtlow <- cast(rtlow,sub~coh,)
rtmed <- with(subset(Simuls,condition=="mediumSC"),aggregate(rt,by=list(sub,coh),mean));names(rtmed) <- c('sub','coh','rt')
rtmed <- cast(rtmed,sub~coh)
rthigh <- with(subset(Simuls,condition=="highSC"),aggregate(rt,by=list(sub,coh),mean));names(rthigh) <- c('sub','coh','rt')
rthigh <- cast(rthigh,sub~coh)


#aggregate rt for model
x <- rtlow[,c(2:4)];xmed <- rtmed[,c(2:4)];xhigh <- rthigh[,c(2:4)]
n <- length(x)

x <- x[,c("hard","average","easy")];
xmed <- xmed[,c("hard","average","easy")];
xhigh <- xhigh[,c("hard","average","easy")]

par(mfrow=c(1,1))
stripchart(x,ylim=c(.6,1.1), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main="",cex.axis=1.25)
mtext("RT (s)",2,at=.85,line=2.5,cex=1.75);axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);mtext("Trial difficulty",1,3,at=1,cex=1.75)
means <- sapply(x, mean);n<- length(x)
legend(0,1.1,legend=c("Negative","Average","Positive"),title = "Fake Feedback condition",pch=rep(16,3),bty = "n",inset=.1, cex = 1.25,col=c("red","orange","blue"))
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="red",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N1),lwd=lwdgr,length=.05,col="red")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="orange",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N1),lwd=lwdgr,length=.05,col="orange")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="blue",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N1),lwd=lwdgr,length=.05,col="blue")


## Experiment 2
#Aggregate conf for Simuls2
rtlow <- with(subset(Simuls2,condition=="hard"),aggregate(rt,by=list(sub,coh),mean));names(rtlow) <- c('sub','coh','rt')
rtlow <- cast(rtlow,sub~coh,)
rtmed <- with(subset(Simuls2,condition=="average"),aggregate(rt,by=list(sub,coh),mean));names(rtmed) <- c('sub','coh','rt')
rtmed <- cast(rtmed,sub~coh)
rthigh <- with(subset(Simuls2,condition=="easy"),aggregate(rt,by=list(sub,coh),mean));names(rthigh) <- c('sub','coh','rt')
rthigh <- cast(rthigh,sub~coh)


#aggregate rt for model
x <- rtlow[,c(2:4)];xmed <- rtmed[,c(2:4)];xhigh <- rthigh[,c(2:4)]
n <- length(x)

x <- x[,c("hard","average","easy")];
xmed <- xmed[,c("hard","average","easy")];
xhigh <- xhigh[,c("hard","average","easy")]

stripchart(x, ylim=c(.6,1.1), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main="",cex.axis=1.25)
mtext("RT (s)",2,at=.85,line=2.5,cex=1.75);axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);mtext("Trial difficulty",1,3,at=1,cex=1.75)
means <- sapply(x, mean);n<- length(x)
legend(0,.9,legend=c("Hard","Average","Easy"),title = "Training difficulty condition",pch=rep(16,3),bty = "n",inset=.1, cex = 1.25,col=c("red","orange","blue"))
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="red",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(Nsub_2),lwd=lwdgr,length=.05,col="red")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="orange",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(Nsub_2),lwd=lwdgr,length=.05,col="orange")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="blue",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(Nsub_2),lwd=lwdgr,length=.05,col="blue")
