library(fields)
library(mnormt)

fast_hm <- function(mu, sigma = .1, dt = .001, ev_bound = .5, ev_window = .01, upperRT = 5){
  ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
  time_vec <- seq(0,upperRT,dt)
  if (length(mu)==1) {
    hm <- sapply(time_vec,function(x){
      if (x==0) {
        res <- rep(NA,length(ev_mapping))
        res[ceiling(length(res)/2)] <- .5
      }else{
        dist_pos <- dnorm(ev_mapping,mean=mu,sd=sigma*sqrt(x))
        dist_neg <- dnorm(ev_mapping,mean=-mu,sd=sigma*sqrt(x))
        res <- dist_pos / (dist_pos + dist_neg)
      }
      c(res)
    })
  }else{
    Nv <- length(mu)
    sigma_mat <- diag(sigma,Nv)
    hm <- sapply(time_vec,function(x){
      if (x==0) {
        res <- rep(NA,length(ev_mapping))
        res[ceiling(length(res)/2)] <- .5
      }else{
        dist_pos <- dmnorm(matrix(rep(ev_mapping,Nv),ncol = Nv),mean=mu,varcov=sigma_mat*sqrt(x))
        dist_neg <- dmnorm(matrix(rep(ev_mapping,Nv),ncol = Nv),mean=-mu,varcov=sigma_mat*sqrt(x))
        res <- dist_pos / (dist_pos + dist_neg)
      }
      c(res)
    })
  }
  return(t(hm))
  
}

## Heat map resolution
dt <- .001; ev_bound <- .5; ev_window <- dt*10; upperRT <- 5
ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
timesteps <- upperRT/dt
v <- .1
Nv <- 100
mu <- rnorm(Nv,v,.02)

hm <- fast_hm(v)

colMap <- viridis(length(hm))

image.plot(1:dim(hm)[1],1:dim(hm)[2],hm,zlim=c(0,1),
          col=colMap,ylab="",xlab='',legend.shrink=.25,
          axes=F)
mtext("Evidence",2,at=dim(hm)[2]/2,line=2);
mtext("Time (s)",1,at=dim(hm)[1]/2,line=2)
axis(1,at=c(1,timesteps),labels=c(0,upperRT));
axis(2,at=c(1,length(ev_mapping)),labels=c(-ev_bound,ev_bound))

hm <- fast_hm(mu)
image.plot(1:dim(hm)[1],1:dim(hm)[2],hm,zlim=c(0,1),
           col=colMap,ylab="",xlab='',legend.shrink=.25,
           axes=F)
mtext("Evidence",2,at=dim(hm)[2]/2,line=2);
mtext("Time (s)",1,at=dim(hm)[1]/2,line=2)
axis(1,at=c(1,timesteps),labels=c(0,upperRT));
axis(2,at=c(1,length(ev_mapping)),labels=c(-ev_bound,ev_bound))


hm <- sapply(mu,function(i){
  temp_hm <- sapply(
    time_vec,function(x){
      if (x==0) {
        res <- rep(NA,length(ev_mapping))
        res[ceiling(length(res)/2)] <- .5
      }else{
        dist_pos <- dnorm(ev_mapping,mean=i,sd=sigma*sqrt(x))
        dist_neg <- dnorm(ev_mapping,mean=-i,sd=sigma*sqrt(x))
        res <- dist_pos / (dist_pos + dist_neg)
      }
      c(res)
    })
})