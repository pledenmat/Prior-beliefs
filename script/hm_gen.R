source("rw_createHM_driftsign.r")

v_min <- .001
v_max <- .5
step <- .001

drifts <- seq(v_min,v_max,step)

# Simulation parameters
dt <- .001; nsim <- 100000; ev_bound <- .5; ev_window <- dt*10; upperRT <- 5
timesteps <- upperRT/dt; ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)

for (v in drifts) {
  print(v)
  if (!(file.exists(paste0("heatmaps/hm_",v,".Rdata")))) {
    mu <- c(-v,v)
    output <- RW_createHM_driftsign(mu, dt=dt, nsim=nsim, ev_bound=ev_bound, ev_window=ev_window, upperRT=upperRT)
    save(output,file=paste0("heatmaps/hm_",v,".Rdata"))
  }
}
v <- drifts[1]
# load(paste0("heatmaps/hm_",v,".Rdata"))
# library(fields)
# 
# ref_output <- output
# 
# ## Fill with linear regression
# for (i in 1:dim(output$upper)[2]) {
#   t <- seq(1:timesteps)
#   fit <- lm(output$upper[,i]~ t)
#   fit.c <- coef(fit)
#   pred <- fit.c[1] + fit.c[2]*t
#   pred[pred>1] <- 1;pred[pred<0] <- 0
#   output$upper[,i][is.na(output$upper[,i])] <- pred[is.na(output$upper[,i])]
# }


## Fill with median for each evidence level
for (d in drifts) {
  if (!(file.exists(paste0("heatmaps/hm_",d,"_filled.Rdata")))) {
    print(d)
    load(paste0("heatmaps/hm_",d,".Rdata"))
    for (i in 1:dim(output$upper)[2]) {
      output$upper[is.na(output$upper[,i]),i] <- median(output$upper[,i],na.rm=T)
      output$lower[is.na(output$lower[,i]),i] <- median(output$lower[,i],na.rm=T)
    }
    save(output,file=paste0("heatmaps/hm_",d,"_filled.Rdata"))
  }
}

hm_up <- output$upper

image.plot(1:dim(hm_up)[1],1:dim(hm_up)[2],hm_up,zlim=c(0,1),ylab="Evidence",xlab='time (s)',main=paste("sub",i),axes=F)
axis(1,at=c(1,timesteps),labels=c(0,upperRT));axis(2,at=c(1,length(ev_mapping)),labels=c(-ev_bound,ev_bound))

for (d in drifts){
  load(paste0("heatmaps/hm_",d,"_filled.Rdata"))
  hm_up <- output$upper
  
  image.plot(1:dim(hm_up)[1],1:dim(hm_up)[2],hm_up,zlim=c(0,1),ylab="Evidence",xlab='time (s)',main=paste("drift =",d),axes=F)
  axis(1,at=c(1,timesteps),labels=c(0,upperRT));axis(2,at=c(1,length(ev_mapping)),labels=c(-ev_bound,ev_bound))
}  
