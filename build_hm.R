library(fields)
library(mnormt)
library(distr)

#' Fast generation of p(cor) heatmap
#'
#' @param mu Drift rate used to generate the heatmap (only correctly works with 1 drift rate for now)
#' @param sigma 
#' @param dt 
#' @param ev_bound 
#' @param ev_window 
#' @param upperRT 
#'
#' @return
#' @export
#'
#' @examples
build_hm <- function(mu, sigma = .1, dt = .001, ev_bound = .5, ev_window = .01, upperRT = 5){
  ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
  time_vec <- seq(0,upperRT,dt)
  if (length(mu)==1) {
    hm <- sapply(time_vec,function(time){
      if (time==0) {
        res <- rep(NA,length(ev_mapping))
        res[ceiling(length(res)/2)] <- .5
      }else{
        dist_pos <- dnorm(ev_mapping,mean=mu*time,sd=sigma*sqrt(time))
        dist_neg <- dnorm(ev_mapping,mean=-mu*time,sd=sigma*sqrt(time))
        res <- dist_pos / (dist_pos + dist_neg)
      }
      c(res)
    })
  }else{
    hm <- sapply(time_vec,function(time){
      if (time==0) {
        res <- rep(NA,length(ev_mapping))
        res[ceiling(length(res)/2)] <- .5
      }else{
        dist_pos <- lapply(mu, function(x){Norm(x*time,sigma*sqrt(time))})
        dist_neg <- lapply(-mu, function(x){Norm(x*time,sigma*sqrt(time))})
        dist_pos <- UnivarMixingDistribution(Dlist = dist_pos)
        dist_neg <- UnivarMixingDistribution(Dlist = dist_neg)
        dpos <- d(dist_pos)
        dneg <- d(dist_neg)
        res <- dpos(ev_mapping)/(dpos(ev_mapping)+dneg(ev_mapping))
        }
      c(res)
    })
  }
  return(t(hm))
}

# test <- create_mixture_gaussian(list(mu = 1,sigma=.1,weight=.5),list(mu=0,sigma=.1,weight=.5))
# test2 <- UnivarMixingDistribution(Norm(mean=1, sd=.1), 
#                          Norm(mean=0, sd=.1),
#                          mixCoeff=c(.5,.5))
# test3 <- UnivarMixingDistribution(Dlist=c(Norm(mean=1, sd=.1), 
#                                   Norm(mean=0, sd=.1)),
#                                   mixCoeff=c(.5,.5))
# truc <- UnivarDistrList(Norm(mean=1, sd=.1), 
#                              Norm(mean=0, sd=.1))
# truc <- UnivarDistrList(truc,Norm(mean=2,sd=.1))
# test4 <- UnivarMixingDistribution(truc)
# sigma <- .1
# mu <- c(.05,.1,.15)
# for (i in 1:length(mu)) {
#   if (i==1) {
#     dist_pos <- UnivarMixingDistribution(Norm(mu[i],sigma))
#     dist_neg <- UnivarMixingDistribution(Norm(-mu[i],sigma))
#   } else if (i==length(mu)) {
#     dist_pos <- UnivarMixingDistribution(dist_pos,Norm(mu[i],sigma),mixCoeff = c((length(mu)-1)/length(mu),1/length(mu)))
#     dist_neg <- UnivarMixingDistribution(dist_neg,Norm(-mu[i],sigma),mixCoeff = c((length(mu)-1)/length(mu),1/length(mu)))
#   } else{
#     dist_pos <- UnivarMixingDistribution(dist_pos,Norm(mu[i],sigma))
#     dist_neg <- UnivarMixingDistribution(dist_neg,Norm(-mu[i],sigma))
#   }
# }
# res <- dpos(ev_mapping)/(dpos+dneg)
# plot(test5)
# 
# 
# gauss_dist_pos <- c()
# gauss_dist_neg <- c()
# time <- 5
# 
# plot(res)
# par(mfrow=c(1,1))
