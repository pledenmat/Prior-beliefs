library(myPackage)

curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir) ## change our current working directory

create_modelHM2 <- function(alpha,beta,type="upper",plot=T,cexax=1.5,cexlab=3,cexmain=2.5,
                           dt = .001, ev_bound = .5, ev_window = .01, upperRT = 5){
  
  ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
  timesteps <- upperRT/dt
  
  hm <- matrix(NA,nrow=timesteps,ncol=length(ev_mapping))
  hm_low <- matrix(NA,nrow=timesteps,ncol=length(ev_mapping))
  hm_up <- matrix(NA,nrow=timesteps,ncol=length(ev_mapping))
  i=1;j=1
  for(e in ev_mapping){
    for(t in seq(0,5,length.out=timesteps)){
      hm[i,j] <- 1/(1 + exp(1/sqrt(t+.00001)*(-alpha*e*custom_sign(e) - beta) ) )
      hm_up[i,j] <- 1/(1+exp((1/sqrt(t+.00001))*(-alpha * e - beta) ))
      hm_low[i,j] <- 1/(1+exp((1/sqrt(t+.00001))*(alpha * e - beta) ))
      i<-i+1
    }
    i<-1;j<-j+1
  }
  #Custom colormap
  
  colors2spaceThrough <- c(rgb(112,48,160,maxColorValue = 255), rgb(146,208,80,maxColorValue = 255),"yellow")
  # colMap <- heat.colors(length(hm))
  colMap <- colorRampPalette(colors2spaceThrough)(length(hm))
  if (plot==T) {
    if (type=="upper") {
      fields::image.plot(1:dim(hm_up)[1],1:dim(hm_up)[2],hm_up,zlim=c(0,1),
                         col=colMap,ylab="",xlab='',legend.shrink=.25,
                         main=paste('alpha = ',round(alpha,3),'; beta =',round(beta,3)),
                         axes=F,cex.main=cexmain,
                         axis.args=list(at=seq(0,1,.5),labels=seq(0,1,.5),cex.axis=cexax))
      mtext("Evidence",2,at=dim(hm_up)[2]/2,line=2,cex=cexlab);
      mtext("Time (s)",1,at=dim(hm_up)[1]/2,line=2,cex=cexlab)
      # axis(1,at=c(1,timesteps),labels=c(0,upperRT),cex.axis=cexax);
      # axis(2,at=c(1,length(ev_mapping)),labels=c(-ev_bound,ev_bound),cex.axis=cexax)
      
    }else if(type=="lower"){
      fields::image.plot(1:dim(hm_low)[1],1:dim(hm_low)[2],hm_low,zlim=c(0,1),col=colMap,ylab="Evidence",xlab='time (s)',main=paste('alpha = ',round(alpha,3),'; beta = ',round(beta,3)),col=colMap,axes=F)
      mtext("Evidence",2,at=dim(hm_low)[2]/2,line=2,cex=cexlab);
      mtext("Time (s)",1,at=dim(hm_low)[1]/2,line=2,cex=cexlab)
      # axis(1,at=c(1,timesteps),labels=c(0,upperRT),cex.axis=cexax,cex.lab=2)
      # axis(2,at=c(1,length(ev_mapping)),labels=c(-ev_bound,ev_bound))
      
    }else if(type=="symmetrical"){ #
      fields::image.plot(1:dim(hm)[1],1:dim(hm)[2],hm,zlim=c(0,1),col=colMap,ylab="Evidence",xlab='time (s)',main=paste('alpha = ',round(alpha,3),'; beta = ',round(beta,3)),col=colMap,axes=F)
      mtext("Evidence",2,at=dim(hm)[2]/2,line=2,cex=cexlab);
      mtext("Time (s)",1,at=dim(hm)[1]/2,line=2,cex=cexlab)
      # axis(1,at=c(1,timesteps),labels=c(0,upperRT))
      # axis(2,at=c(1,length(ev_mapping)),labels=c(-ev_bound,ev_bound))
    }
  }else{
    if (type=="upper") {
      return(hm_up)
    }else if(type=="lower"){
      return(hm_low)
    }else if(type=="symmetrical"){
      return(hm)
    }
  }
}

jpeg(
  filename="hm_low.jpeg",
  width=9,
  height=7,
  units="in",
  res=500)
create_modelHM2(5,0,dt=.0025,ev_window=.0125)
dev.off()

jpeg(
  filename="hm_med.jpeg",
  width=9,
  height=7,
  units="in",
  res=500)
create_modelHM2(15,0,dt=.0025,ev_window=.0125)
dev.off()

jpeg(
  filename="hm_high.jpeg",
  width=9,
  height=7,
  units="in",
  res=500)
create_modelHM2(40,0,dt=.0025,ev_window=.0125)
dev.off()
