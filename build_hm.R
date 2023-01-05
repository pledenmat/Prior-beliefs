library(fields)
library(mnormt)

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
    hm <- sapply(time_vec,function(x){
      if (x==0) {
        res <- rep(NA,length(ev_mapping))
        res[ceiling(length(res)/2)] <- .5
      }else{
        dist_pos <- dnorm(ev_mapping,mean=mu*x,sd=sigma*sqrt(x))
        dist_neg <- dnorm(ev_mapping,mean=-mu*x,sd=sigma*sqrt(x))
        res <- dist_pos / (dist_pos + dist_neg)
      }
      c(res)
    })
  }else{
    print("Case with length(mu) > 1 is not implemented yet")
  }
  return(t(hm))
  
}
