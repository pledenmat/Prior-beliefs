error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
cexkl <- 2;cexgr <- 2;lwdgr <- 3; lwddat <- 2

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

par(mfrow=c(1,1))
stripchart(x,ylim=c(.45,.95), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main="Model prediction",cex.axis=1.25)
mtext("Confidence",2,at=.7,line=2.5,cex=1.75);axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);mtext("Trial difficulty",1,3,at=1,cex=1.75)
means <- sapply(x, mean);n<- length(x)
legend(0,.95,legend=c("Negative","Average","Positive"),title = "Fake Feedback condition",pch=rep(16,3),bty = "n",inset=.1, cex = 1.25,col=c("red","orange","blue"))
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="red",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="red")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="orange",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="orange")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="blue",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="blue")


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

stripchart(x, ylim=c(.45,.95), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main="Model prediction",cex.axis=1.25)
mtext("Confidence",2,at=.7,line=2.5,cex=1.75);axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);mtext("Trial difficulty",1,3,at=1,cex=1.75)
means <- sapply(x, mean);n<- length(x)
legend(0,.95,legend=c("Hard","Average","Easy"),title = "Training difficulty condition",pch=rep(16,3),bty = "n",inset=.1, cex = 1.25,col=c("red","orange","blue"))
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="red",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="red")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="orange",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="orange")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cexkl,col="blue",lwd=lwddat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="blue")
