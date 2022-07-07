curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)
library(reshape) #cast
library(timeSeries) #colSds
Data1 <- read.csv('data_exp1.csv')
Data2 <- read.csv('data_exp2.csv')
Simuls <- read.csv("model_prediction_exp1.csv")
Simuls2 <- read.csv("model_prediction_exp2.csv")
load("param_ddm_test_exp1.Rdata")
load("param_ddm_test_exp2.Rdata")
load("param_train_exp1.Rdata")
load("param_train_exp2.Rdata")

subs1 <- sort(unique(Data1$sub)); N1 <- length(subs1) 
cond_1 <- sort(unique(Data1$selfconf)); Ncond <- length(cond_1)
coh <- sort(unique(Data1$coh));Ndiff <- length(coh)

subs_2 <- sort(unique(Data2$sub)); Nsub_2 <- length(subs_2)
cond_2 <- sort(unique(Data2$traindiffcond)); Ncond_2 <- length(cond_2)
coh <- sort(unique(Data2$coh));Ndiff <- length(coh)
# Functions and parameters ----
## Transparent colors, Mark Gardener 2015, www.dataanalytics.org.uk
transp <- function(color, percent = 50, name = NULL) {
  #   color = color name
  #   percent = % transparency
  #   name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  transp <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                max = 255,
                alpha = (100 - percent) * 255 / 100,
                names = name)
  
  ## Save the color
  invisible(transp)
}

error.bar <- function(x, y, upper, lower=upper, length=0,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
cexkl <- 1.5;cexgr <- 2;lwdgr <- 3; lwddat <- 2
cex_lab <- 3; cex_legend <- 3; cex_title <- 2.5 
windowsFonts(A = windowsFont("Calibri")) 
par(family="A",font.main = 2, cex.main = cex_title)

### FIGURE 3
# Plot Layout -------------------------------------------------------------
jpeg(
  filename="results.jpeg",
  width=13,
  height=18,
  units="in",
  res=500)
# layout(matrix(c(1,3,7,9,11,12,1,5,7,9,11,12,2,4,8,10,11,13,2,6,8,10,11,13),ncol=4),heights = c(.4,1,1.5,1.5,.2,1.5))
layout(matrix(c(1,3,7,9,10,1,5,7,9,10,2,4,8,9,11,2,6,8,9,11),ncol=4),heights = c(.4,1,2,.2,2))

#' Add legend for both experiments on top
par(mar=c(0,0,0,0))
plot.new()
legend("top",legend=c("Negative","Average","Positive"),
       title = "Exp. 1: Fake feedback",pch=rep(16,3),bty = "n",inset=0,
       cex = cex_legend,col=c("brown3","cyan4","darkgoldenrod3"), horiz = T)

plot.new()
legend("top",legend=c("Difficult","Medium","Easy"),
       title = "Exp. 2: Training difficulty",pch=rep(16,3),bty = "n",inset=0, 
       cex = cex_legend,col=c("brown3","cyan4","darkgoldenrod3"), horiz = T)
par(mar=c(5,5,2,2)+0.1)
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
           frame=F,xaxt='n',yaxt='n',xlab="Trial difficulty",ylab="Accuracy",
           cex.lab=cex_lab/2)
mtext("A.", at = -.55, line = 1, cex = cex_title, font = 2)
axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);
axis(2,seq(0.6,1,.1),cex.axis=1.5)
# mtext("Accuracy",2,at=.8,line=2.5,cex=cex_lab);
# mtext("Trial difficulty",1,3,at=1,cex=cex_lab)
for(i in seq(.6,1,.1)) abline(h=i,col="lightgrey",lty = "dashed")
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N1)),
          (colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(205,51,51,51,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xmed_sim,na.rm=T) + (colSds(as.matrix(xmed_sim))/sqrt(N1)),
          (colMeans(xmed_sim,na.rm=T) - colSds(as.matrix(xmed_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(0,139,139,51,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xhigh_sim,na.rm=T) + (colSds(as.matrix(xhigh_sim))/sqrt(N1)),
          (colMeans(xhigh_sim,na.rm=T) - colSds(as.matrix(xhigh_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(205,149,12,51,maxColorValue = 255))
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
           frame=F,xaxt='n',main=NULL,yaxt='n',xlab="Trial difficulty",ylab="Accuracy",cex.lab=cex_lab/2)
mtext("B.", at = -.55, line = 1, cex = cex_title, font = 2)
axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);
axis(2,seq(0.6,1,.1),cex.axis=1.5)
# mtext("Accuracy",2,at=.8,line=2.5,cex=cex_lab);
# mtext("Trial difficulty",1,3,at=1,cex=cex_lab)
for(i in seq(.6,1,.1)) abline(h=i,col="lightgrey",lty = "dashed")
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N1)),
          (colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(205,51,51,51,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xmed_sim,na.rm=T) + (colSds(as.matrix(xmed_sim))/sqrt(N1)),
          (colMeans(xmed_sim,na.rm=T) - colSds(as.matrix(xmed_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(0,139,139,51,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xhigh_sim,na.rm=T) + (colSds(as.matrix(xhigh_sim))/sqrt(N1)),
          (colMeans(xhigh_sim,na.rm=T) - colSds(as.matrix(xhigh_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(205,149,12,51,maxColorValue = 255))
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
           frame=F,xaxt='n',main=NULL,yaxt='n',xlab="Trial difficulty",ylab="RT (s)",cex.lab=cex_lab/2)
axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);
axis(2,seq(0.6,1.1,.1),cex.axis=1.5)
# mtext("RT (s)",2,at=.85,line=2.5,cex=cex_lab);
# mtext("Trial difficulty",1,3,at=1,cex=cex_lab)
for(i in seq(.6,1.1,.1)) abline(h=i,col="lightgrey",lty = "dashed")

polygon(c(0:(n-1),(n-1):0),
        c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N1)),
          (colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(205,51,51,51,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xmed_sim,na.rm=T) + (colSds(as.matrix(xmed_sim))/sqrt(N1)),
          (colMeans(xmed_sim,na.rm=T) - colSds(as.matrix(xmed_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(0,139,139,51,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xhigh_sim,na.rm=T) + (colSds(as.matrix(xhigh_sim))/sqrt(N1)),
          (colMeans(xhigh_sim,na.rm=T) - colSds(as.matrix(xhigh_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(205,149,12,51,maxColorValue = 255))
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
           frame=F,xaxt='n',main=NULL,yaxt='n',xlab="Trial difficulty",ylab="RT (s)",cex.lab=cex_lab/2)
axis(1,at=0:(n-1),labels=names(x), cex.axis=1.5);
axis(2,seq(0.6,1.1,.1),cex.axis=1.5)
# mtext("RT (s)",2,at=.85,line=2.5,cex=cex_lab);
# mtext("Trial difficulty",1,3,at=1,cex=cex_lab)
for(i in seq(.6,1.1,.1)) abline(h=i,col="lightgrey",lty = "dashed")
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N1)),
          (colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(205,51,51,51,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xmed_sim,na.rm=T) + (colSds(as.matrix(xmed_sim))/sqrt(N1)),
          (colMeans(xmed_sim,na.rm=T) - colSds(as.matrix(xmed_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(0,139,139,51,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xhigh_sim,na.rm=T) + (colSds(as.matrix(xhigh_sim))/sqrt(N1)),
          (colMeans(xhigh_sim,na.rm=T) - colSds(as.matrix(xhigh_sim))/sqrt(N1))[3:1]),
        border=F,col=rgb(205,149,12,51,maxColorValue = 255))
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
##' Experiment 1
CJ_SC_diff_data <- with(Data1,aggregate(cj,by=list(sub=sub,selfconf=selfconf, coh=coh),mean));
CJ_SC_diff_data <- cast(CJ_SC_diff_data,sub~selfconf+coh)
average_CJ_SC_diff_data <- with(Data1,aggregate(cj,by=list(selfconf=selfconf,coh=coh),mean));
average_CJ_SC_diff_data <- cast(average_CJ_SC_diff_data,selfconf~coh)

# use family to adjust the font and cex. to adjust font size
CJ_SC_diff_plot = plot(as.numeric(average_CJ_SC_diff_data[1,]),type='n',frame=F,
                       main=NULL,
                       ylab="Confidence",
                       xlab="",
                       xaxt='n',
                       xlim=c(1,3.3),ylim=c(3,6),
                       cex.axis = cex_lab-1, 
                       cex.lab = cex_lab,
                       family="A")
mtext("C.", at=.7, line = 1, cex = cex_title, font = 2)
axis(1,at=1.1:3.1,labels=c("hard","average","easy"),cex.axis=cex_lab-1,family="A")
abline(h = seq(3,6,0.5), col = "lightgrey", lty = "dashed")

# High SC
for(i in 1:N1) points(jitter(1:3,0.1),CJ_SC_diff_data[i,c(4,2,3)],lty=i,type='p',pch=21,col='white',bg=transp('darkgoldenrod2'))
lines(1:3,average_CJ_SC_diff_data[1,c(4,2,3)],lty=2,type='b',pch=21,
      col='darkgoldenrod3',bg='darkgoldenrod2',lwd=lwddat,cex=cexkl)
# Medium SC
for(i in 1:N1) points(jitter(1.1:3.1,0.1),CJ_SC_diff_data[i,c(10,8,9)],lty=i,type='p',pch=24,col='white',bg=transp('cyan3'))
lines(1.1:3.1,average_CJ_SC_diff_data[3,c(4,2,3)],lty=2,type='b',pch=24,
      col='cyan4',bg='cyan3',lwd=lwddat,cex=cexkl)
# Low SC
for(i in 1:N1) points(jitter(1.2:3.2,0.1),CJ_SC_diff_data[i,c(7,5,6)],lty=i,type='p',pch=22,col='white',bg=transp('brown2'))
lines(1.2:3.2,average_CJ_SC_diff_data[2,c(4,2,3)],lty=2,type='b',pch=22,
      col='brown3',bg="brown2",lwd=lwddat,cex=cexkl)

# plot error bars
error.bar(1:3,colMeans(CJ_SC_diff_data[,c(4,2,3)]),colSds(as.matrix(CJ_SC_diff_data[,c(4,2,3)])/sqrt(N1)),
          length=0,lwd=lwdgr, col='darkgoldenrod3')
error.bar(1.1:3.1,colMeans(CJ_SC_diff_data[,c(10,8,9)]),colSds(as.matrix(CJ_SC_diff_data[,c(10,8,9)])/sqrt(N1)),
          length=0,lwd=lwdgr, col='cyan4')
error.bar(1.2:3.2,colMeans(CJ_SC_diff_data[,c(7,5,6)]),colSds(as.matrix(CJ_SC_diff_data[,c(7,5,6)])/sqrt(N1)),
          length=0,lwd=lwdgr, col='brown3')


## Experiment 2
CJ_SC_diff_data <- with(Data2,aggregate(cj,by=list(sub=sub,traindiffcond=traindiffcond, coh=coh),mean));
CJ_SC_diff_data <- cast(CJ_SC_diff_data,sub~traindiffcond+coh)
average_CJ_SC_diff_data <- with(Data2,aggregate(cj,by=list(traindiffcond=traindiffcond,coh=coh),mean));
average_CJ_SC_diff_data <- cast(average_CJ_SC_diff_data,traindiffcond~coh)

# use family to adjust the font and cex. to adjust font size
CJ_SC_diff_plot = plot(as.numeric(average_CJ_SC_diff_data[1,]),type='n',frame=F,
                       main=NULL,
                       ylab="Confidence",
                       xlab="",
                       xaxt='n',
                       xlim=c(1,3.3),ylim=c(3,6),
                       cex.axis = cex_lab-1, 
                       cex.lab = cex_lab,
                       family="A")
mtext("D.", at = .7, line = 1, cex = cex_title, font = 2)
axis(1,at=1.1:3.1,labels=c("hard","average","easy"),cex.axis=cex_lab-1,family="A")
abline(h = seq(3,6,0.5), col = "lightgrey", lty = "dashed")

# High SC
for(i in 1:Nsub_2) points(jitter(1:3,0.1),CJ_SC_diff_data[i,c(7,5,6)],lty=i,type='p',pch=21,col='white',bg=transp('darkgoldenrod2'))
lines(1:3,average_CJ_SC_diff_data[2,c(4,2,3)],lty=2,type='b',pch=21,
      col='darkgoldenrod3',bg='darkgoldenrod2',lwd=lwddat,cex=cexkl)
# Medium SC
for(i in 1:Nsub_2) points(jitter(1.1:3.1,0.1),CJ_SC_diff_data[i,c(4,2,3)],lty=i,type='p',pch=24,col='white',bg=transp('cyan3'))
lines(1.1:3.1,average_CJ_SC_diff_data[1,c(4,2,3)],lty=2,type='b',pch=24,
      col='cyan4',bg='cyan3',lwd=lwddat,cex=cexkl)
# Low SC
for(i in 1:Nsub_2) points(jitter(1.2:3.2,0.1),CJ_SC_diff_data[i,c(10,8,9)],lty=i,type='p',pch=22,col='white',bg=transp('brown2'))
lines(1.2:3.2,average_CJ_SC_diff_data[3,c(4,2,3)],lty=2,type='b',pch=22,
      col='brown3',bg="brown2",lwd=lwddat,cex=cexkl)

# plot error bars
error.bar(1:3,colMeans(CJ_SC_diff_data[,c(7,5,6)]),colSds(as.matrix(CJ_SC_diff_data[,c(7,5,6)])/sqrt(Nsub_2)),
          length=0,lwd=lwdgr, col='darkgoldenrod3')
error.bar(1.1:3.1,colMeans(CJ_SC_diff_data[,c(4,2,3)]),colSds(as.matrix(CJ_SC_diff_data[,c(4,2,3)])/sqrt(Nsub_2)),
          length=0,lwd=lwdgr, col='cyan4')
error.bar(1.2:3.2,colMeans(CJ_SC_diff_data[,c(10,8,9)]),colSds(as.matrix(CJ_SC_diff_data[,c(10,8,9)])/sqrt(Nsub_2)),
          length=0,lwd=lwdgr, col='brown3')


# Plot confidence prediction ----------------------------------------------
par(mar=c(0,0,0,0))
plot.new()
text(.5,.75, labels="E. Model Predictions",cex = cex_legend+.5,font=2)
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
           main=NULL, yaxt = 'n',family="A",ylab = "Confidence",cex.lab=cex_lab )
# mtext("Confidence",2,at=.75,line=2.5,cex=cex_lab);
# mtext("Trial difficulty",1,3,at=1,cex=cex_lab)
# mtext("E.", at=.7, line = 1, cex = cex_title, font = 2)
axis(1,at=0:(n-1),labels=names(x), cex.axis=cex_lab-1);
axis(2, seq(.5,1,.1), cex.axis=cex_lab-1)
means <- sapply(x, mean);n<- length(x)
for(i in seq(.5,1,length.out = 6)) abline(h=i,col="lightgrey",lty = "dashed")
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
           main=NULL, yaxt = 'n',family="A",ylab="Confidence",cex.lab=cex_lab)
# mtext("Confidence",2,at=.75,line=2.5,cex=cex_lab);
# mtext("Trial difficulty",1,3,at=1,cex=cex_lab)
# mtext("F.", at=.7, line = 1, cex = cex_title, font = 2)
axis(1,at=0:(n-1),labels=names(x), cex.axis=cex_lab-1);
axis(2, seq(.5,1,.1), cex.axis=cex_lab-1)
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

#'save
dev.off()
### FIGURE 4
# EXP 1 Estimated parameters plot ----------------------------------------------
# Layout ------------------------------------------------------------------
jpeg(
  filename="fit_results1.jpeg",
  width=13,
  height=11,
  units="in",
  res=500)
# layout(matrix(c(1,2,1,3,1,4),ncol=3),heights = c(1.5,1))
layout(matrix(c(1,2,3,1,5,4),ncol=2),heights = c(.2,1,1))
par(mar=c(0,0,0,0))
plot.new()
text(.5,.75, labels="Exp. 1: Fake feedback",cex = cex_legend+.5,font=2)
par(mar=c(5,5,4,0)+0.1)
# Plot Subjective drift -----------------------------------------------------
##Exp1
plot_drift <- with(df,aggregate(Vs,by=list(sub=sub,condition=condition),mean))
plot_drift <- cast(plot_drift,sub~condition)
plot_drift <- plot_drift[,c(3,4,2)] #Reorder columns to have hard -> easy
plot(colMeans(plot_drift),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylab='Subjective drift rate',ylim=c(0,.16),
     xlab="Feedback condition",xaxt='n', yaxt='n')
mtext("A.", at = .55, line = 3, cex = cex_title, font = 2)
segments(y0 = seq(0,.16,.04),y1 = seq(0,.16,.04),x0 = 0, x1 = Ncond, col = "lightgrey", lty = "dotted")
axis(1,1:Ncond,c("Negative","Average","Positive"),cex.axis=1.75)
axis(2,at=seq(0,.16,.04),cex.axis=1.75)
for(i in 1:Nsub_2){
  x <- jitter(1:Ncond,.2)
  lines(x,plot_drift[i,1:Ncond],lty=2,col=transp('grey'))
  points(x,plot_drift[i,1:Ncond],col="white", bg = transp('grey'),pch=21, cex = 2)
} 
points(colMeans(plot_drift),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_drift),colSds(plot_drift,na.rm=T)/sqrt(N1),lwd=3)

# Plot DDM parameters test phase EXP1 ------------------------------------------
par(mar=c(5,5,2,0)+0.1)
##Non-decision time
plot_ter <- with(param_ddm_test_exp1,aggregate(ter,by=list(sub=sub,condition=condition),mean))
plot_ter <- cast(plot_ter,sub~condition)
plot_ter <- plot_ter[,c(3,4,2)] #Reorder columns to have easy -> hard
plot(colMeans(plot_ter),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylim=c(.15,.8),ylab="Non-decision time",
     xlab="Feedback condition",xaxt='n',main="", cex.main = 2);
mtext("C.", at = .55, line = 1, cex = cex_title, font = 2)
axis(1,1:Ncond,c("Negative","Average","Positive"),cex.axis=1.75)
segments(y0 = seq(.2,.8,.1),y1 = seq(.2,.8,.1),x0 = 0, x1 = Ncond, col = "lightgrey", lty = "dotted")
for(i in 1:Nsub_2){
  x <- jitter(1:Ncond,.2)
  lines(x,plot_ter[i,1:Ncond],lty=2,col=transp('grey'))
  points(x,plot_ter[i,1:Ncond],col="white", bg = transp('grey'),pch=21, cex = 2)
} 
points(colMeans(plot_ter),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_ter),colSds(plot_ter,na.rm=T)/sqrt(N1),lwd=3,length=0)

##Bound
plot_bound <- with(param_ddm_test_exp1,aggregate(bound,by=list(sub=sub,condition=condition),mean))
plot_bound <- cast(plot_bound,sub~condition)
plot_bound <- plot_bound[,c(3,4,2)] #Reorder columns to have easy -> hard
plot(colMeans(plot_bound),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylim=c(0,.2),ylab="Bound",
     xlab="Feedback condition",xaxt='n',main="", cex.main = 2);
mtext("D.", at = .55, line = 1, cex = cex_title, font = 2)
segments(y0 = seq(0,.2,.05),y1 = seq(0,.2,.05),x0 = 0, x1 = Ncond, col = "lightgrey", lty = "dotted")
axis(1,1:Ncond,c("Negative","Average","Positive"),cex.axis=1.75)
for(i in 1:Nsub_2){
  x <- jitter(1:Ncond,.2)
  lines(x,plot_bound[i,1:Ncond],lty=2,col=transp('grey'))
  points(x,plot_bound[i,1:Ncond],col="white", bg = transp('grey'),pch=21, cex = 2)
} 
points(colMeans(plot_bound),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_bound),colSds(plot_bound,na.rm=T)/sqrt(N1),lwd=3,length=0)

##Drift interaction
plot_drift_minus <- with(subset(param_ddm_test_exp1,difflevel=="hard"),
                         aggregate(drift,by=list(sub=sub,condition=condition),mean))
plot_drift_minus <- cast(plot_drift_minus,sub~condition)
plot_drift_minus <- plot_drift_minus[,c(3,4,2)] #Reorder columns to have easy -> hard
plot(colMeans(plot_drift_minus),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,xlim=c(.8,Ncond+.2),
     ylim=c(min(plot_drift_minus),.27),ylab="Drift rate",xlab="Feedback condition",xaxt='n');
axis(1,1:Ncond,c("Negative","Average","Positive"),cex.axis=1.75)
mtext("B.", at = .55, line = 1, cex = cex_title, font = 2)
segments(y0 = seq(0,.25,.05),y1 = seq(0,.25,.05),x0 = 0, x1 = Ncond, col = "lightgrey", lty = "dotted")
points(colMeans(plot_drift_minus),type='b',lwd=5,col="darkolivegreen",lty="dashed")
error.bar(1:Ncond,colMeans(plot_drift_minus),
          colSds(plot_drift_minus,na.rm=T)/sqrt(N1),lwd=3,length=0,col="darkolivegreen")

plot_drift_control <- with(subset(param_ddm_test_exp1,difflevel=="average"),
                           aggregate(drift,by=list(sub=sub,condition=condition),mean))
plot_drift_control <- cast(plot_drift_control,sub~condition)
plot_drift_control <- plot_drift_control[,c(3,4,2)] #Reorder columns to have easy -> hard
points(colMeans(plot_drift_control),type='b',lwd=5,col="darkolivegreen3",lty="dotdash")
error.bar(1:Ncond,colMeans(plot_drift_control),
          colSds(plot_drift_control,na.rm=T)/sqrt(N1),lwd=3,length=0,col="darkolivegreen3")

plot_drift_plus <- with(subset(param_ddm_test_exp1,difflevel=="easy"),
                        aggregate(drift,by=list(sub=sub,condition=condition),mean))
plot_drift_plus <- cast(plot_drift_plus,sub~condition)
plot_drift_plus <- plot_drift_plus[,c(3,4,2)] #Reorder columns to have easy -> hard
points(colMeans(plot_drift_plus),type='b',lwd=5,col="darkolivegreen1",)
error.bar(1:Ncond,colMeans(plot_drift_plus),
          colSds(plot_drift_plus,na.rm=T)/sqrt(N1),lwd=3,length=0,col="darkolivegreen1")
par(xpd=T)
legend("top",border=F,legend=c("Hard","Average","Easy"),lwd=3, horiz = T, inset = c(0,-.065),
       col=c("darkolivegreen","darkolivegreen3","darkolivegreen1"),bty="n",cex=2,
       title = "Trial difficulty", lty = c("dashed","dotdash","solid"))
par(xpd=F)

#'save
dev.off()
# EXP 2 Estimated parameters plot -----------------------------------------
# Layout ------------------------------------------------------------------
jpeg(
  filename="fit_results2.jpeg",
  width=13,
  height=11,
  units="in",
  res=500)
layout(matrix(c(1,2,3,1,5,4),ncol=2),heights = c(.2,1,1))
par(mar=c(0,0,0,0))
plot.new()
text(.5,.75, labels="Exp. 2: Training difficulty",cex = cex_legend+.5,font=2)
par(mar=c(5,5,4,0)+0.1)
# Plot Subjective drift -----------------------------------------------------
##Exp1
plot_drift <- with(df2,aggregate(Vs,by=list(sub=sub,condition=condition),mean))
plot_drift <- cast(plot_drift,sub~condition)
plot_drift <- plot_drift[,c(4,2,3)] #Reorder columns to have hard -> easy
plot(colMeans(plot_drift),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylab='Subjective drift rate',ylim=c(0,.5),
     xlab="Training condition",xaxt='n', yaxt='n')
mtext("E.", at = .55, line = 3, cex = cex_title, font = 2)
segments(y0 = seq(0,.5,.1),y1 = seq(0,.5,.1),x0 = 0, x1 = Ncond, col = "lightgrey", lty = "dotted")
axis(1,1:Ncond,c("Difficult","Medium","Easy"),cex.axis=1.75)
axis(2,at=seq(0,.5,.1),cex.axis=1.75)
for(i in 1:Nsub_2){
  x <- jitter(1:Ncond,.2)
  lines(x,plot_drift[i,1:Ncond],lty=2,col=transp('grey'))
  points(x,plot_drift[i,1:Ncond],col="white", bg = transp('grey'),pch=21, cex = 2)
} 
points(colMeans(plot_drift),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_drift),colSds(plot_drift,na.rm=T)/sqrt(Nsub_2),lwd=3)

# Plot DDM parameters test phase EXP2 ------------------------------------------
par(mar=c(5,5,2,0)+0.1)
##Non-decision time
plot_ter <- with(param_ddm_test_exp2,aggregate(ter,by=list(sub=sub,condition=condition),mean))
plot_ter <- cast(plot_ter,sub~condition)
plot_ter <- plot_ter[,c(4,2,3)] #Reorder columns to have easy -> hard
plot(colMeans(plot_ter),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylim=c(.15,.8),ylab="Non-decision time",
     xlab="Training condition",xaxt='n',main="", cex.main = 2);
mtext("G.", at = .55, line = 1, cex = cex_title, font = 2)
axis(1,1:Ncond,c("Difficult","Medium","Easy"),cex.axis=1.75)
segments(y0 = seq(.2,.8,.1),y1 = seq(.2,.8,.1),x0 = 0, x1 = Ncond, col = "lightgrey", lty = "dotted")
for(i in 1:Nsub_2){
  x <- jitter(1:Ncond,.2)
  lines(x,plot_ter[i,1:Ncond],lty=2,col=transp('grey'))
  points(x,plot_ter[i,1:Ncond],col="white", bg = transp('grey'),pch=21, cex = 2)
} 
points(colMeans(plot_ter),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_ter),colSds(plot_ter,na.rm=T)/sqrt(Nsub_2),lwd=3,length=0)

##Bound
plot_bound <- with(param_ddm_test_exp2,aggregate(bound,by=list(sub=sub,condition=condition),mean))
plot_bound <- cast(plot_bound,sub~condition)
plot_bound <- plot_bound[,c(4,2,3)] #Reorder columns to have easy -> hard
plot(colMeans(plot_bound),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,
     xlim=c(.8,Ncond+.2),ylim=c(.04,.12),ylab="Bound",
     xlab="Training condition",xaxt='n',main="", cex.main = 2);
mtext("H.", at = .55, line = 1, cex = cex_title, font = 2)
segments(y0 = seq(.04,.12,.02),y1 = seq(.04,.12,.02),x0 = 0, x1 = Ncond, col = "lightgrey", lty = "dotted")
axis(1,1:Ncond,c("Difficult","Medium","Easy"),cex.axis=1.75)
for(i in 1:Nsub_2){
  x <- jitter(1:Ncond,.2)
  lines(x,plot_bound[i,1:Ncond],lty=2,col=transp('grey'))
  points(x,plot_bound[i,1:Ncond],col="white", bg = transp('grey'),pch=21, cex = 2)
} 
points(colMeans(plot_bound),type='b',lwd=5)
error.bar(1:Ncond,colMeans(plot_bound),colSds(plot_bound,na.rm=T)/sqrt(Nsub_2),lwd=3,length=0)

##Drift interaction
plot_drift_minus <- with(subset(param_ddm_test_exp2,difflevel=="hard"),
                         aggregate(drift,by=list(sub=sub,condition=condition),mean))
plot_drift_minus <- cast(plot_drift_minus,sub~condition)
plot_drift_minus <- plot_drift_minus[,c(4,2,3)] #Reorder columns to have easy -> hard
plot(colMeans(plot_drift_minus),frame=F,type='n',cex.lab=2.5,cex.axis=1.75,xlim=c(.8,Ncond+.2),
     ylim=c(min(plot_drift_minus),.27),ylab="Drift rate",xlab="Training condition",xaxt='n',yaxt='n');
axis(1,1:Ncond,c("Difficult","Medium","Easy"),cex.axis=1.75)
axis(2,seq(0,.25,.05),cex.axis=1.75)
mtext("F.", at = .55, line = 1, cex = cex_title, font = 2)
segments(y0 = seq(0,.25,.05),y1 = seq(0,.25,.05),x0 = 0, x1 = Ncond, col = "lightgrey", lty = "dotted")
points(colMeans(plot_drift_minus),type='b',lwd=5,col="darkolivegreen",lty="dashed")
error.bar(1:Ncond,colMeans(plot_drift_minus),
          colSds(plot_drift_minus,na.rm=T)/sqrt(Nsub_2),lwd=3,length=0,col="darkolivegreen")

plot_drift_control <- with(subset(param_ddm_test_exp2,difflevel=="average"),
                           aggregate(drift,by=list(sub=sub,condition=condition),mean))
plot_drift_control <- cast(plot_drift_control,sub~condition)
plot_drift_control <- plot_drift_control[,c(4,2,3)] #Reorder columns to have easy -> hard
points(colMeans(plot_drift_control),type='b',lwd=5,col="darkolivegreen3",lty="dotdash")
error.bar(1:Ncond,colMeans(plot_drift_control),
          colSds(plot_drift_control,na.rm=T)/sqrt(Nsub_2),lwd=3,length=0,col="darkolivegreen3")

plot_drift_plus <- with(subset(param_ddm_test_exp2,difflevel=="easy"),
                        aggregate(drift,by=list(sub=sub,condition=condition),mean))
plot_drift_plus <- cast(plot_drift_plus,sub~condition)
plot_drift_plus <- plot_drift_plus[,c(4,2,3)] #Reorder columns to have easy -> hard
points(colMeans(plot_drift_plus),type='b',lwd=5,col="darkolivegreen1",)
error.bar(1:Ncond,colMeans(plot_drift_plus),
          colSds(plot_drift_plus,na.rm=T)/sqrt(Nsub_2),lwd=3,length=0,col="darkolivegreen1")
par(xpd=T)
legend("top",border=F,legend=c("Hard","Average","Easy"),lwd=3, horiz = T, inset = c(0,-.07),
       col=c("darkolivegreen","darkolivegreen3","darkolivegreen1"),bty="n",cex=2,
       title = "Trial difficulty", lty = c("dashed","dotdash","solid"))
par(xpd=F)

#'save
dev.off()
