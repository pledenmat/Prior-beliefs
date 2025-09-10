#' quantile_optim_DDM_Vs_bias
#'
#' @param params Parameters of the model to be fitted + fixed parameters. Order : bound, ter, z, ntrials, sigma, dt, vratio, v_s, bias, drift(s)
#' @param observations Data to be fitted
#' @param returnFit if 0, then returns predicted data. Else, fit the parameters
#' @param confRT_name Name of the column in observations containing the confidence reaction time
#' @param ev_bound Min and Max of the evidence space of the heatmap (goes from -ev_bound to ev_bound)
#' @param ev_window Evidence resolution of the heatmap
#' @param upperRT Maximum time (in seconds) computed by the heatmap
#' @param conf_min Scale of confidence reports (default to 0 = from sure error to sure correct)
#' @param binning Bin confidence to 6 point scale if true, continuous confidence otherwise
#'
#' @return predicted data if returnFit == 0, else returns the prediction error
#'
#' @details This is the function used to simultaneously fit reaction time, accuracy and confidence judgments via quantile optimization.
#'
#' @examples
#' @export
library(MALDIquant)
library(Rcpp)
library(myPackage)
sourceCpp("DDM_with_confidence_slow_fullconfRT.cpp")
source('build_hm.R')
source('fastmerge.R')

quantile_optim_DDM_Vs_bias <- function(params, observations, returnFit,confRT_name = "RTconf",condition_name = "selfconf",
                                       ev_bound = .5, ev_window = .01, upperRT = 5, conf_min = 0, binning = T){ # Heatmap parameters
  #First, generate predictions:
  drift <- params[10:length(params)]
  params <- params[1:9]
  # v_s = subjective drift/prior belief parameter
  names(params) <- c('a','ter','z','ntrials','sigma','dt','vratio','v_s','bias')
  
  coherences <- sort(unique(observations$coh)) #/!\ CHANGE ACCORDING TO DATASET
  
  # Generate trials from DDM parameters
  trial = data.frame(matrix(NA,nrow=0,ncol=7))
  names(trial) <- c('rt','resp','cor','evidence2','rt2','cj','drift')
  for (d in 1:length(drift)) {
    predictions <- data.frame(DDM_with_confidence_slow_fullconfRT(
      v=drift[d],a=params['a'],ter=params['ter'],z=params['z'],
      ntrials=params['ntrials']*dim(observations)[1]/2/length(drift),s=params['sigma'],
      dt=params['dt'],t2distribution=rep(observations[observations$coh==coherences[d],confRT_name],times=params['ntrials']),
      postdriftmod=params['vratio']))
    predictionsneg <- data.frame(DDM_with_confidence_slow_fullconfRT(
      v=-drift[d],a=params['a'],ter=params['ter'],z=params['z'],
      ntrials=params['ntrials']*dim(observations)[1]/2/length(drift),s=params['sigma'],
      dt=params['dt'],t2distribution=rep(observations[observations$coh==coherences[d],confRT_name],times=params['ntrials']),
      postdriftmod=params['vratio']))
    names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
    names(predictionsneg) <- c('rt','resp','cor','evidence2','rt2','cj')
    predictions <- fastmerge(predictions,predictionsneg)
    predictions['drift'] <- drift[d]
    trial <- fastmerge(trial,predictions)
  }
  predictions <- trial
  
  # Generate the heatmap from v_s
  hm_up <- fast_hm(params["v_s"],sigma = params["sigma"], dt = params["dt"], 
                   ev_bound = ev_bound, ev_window = ev_window, upperRT = upperRT)
  hm_low <- 1-hm_up
  hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
  
  # match data to the heatmap
  timesteps <- upperRT/params["dt"]
  ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
  
  predictions$closest_evdnc2 <- match.closest(predictions$evidence2,ev_mapping)
  predictions$temprt2 <- predictions$rt2;
  predictions$temprt2[predictions$temprt2>upperRT] <- upperRT #heatmap doesn't go higher
  predictions$temprt2 <- predictions$temprt2*timesteps/upperRT #scale with the heatmap
  
  # Compute confidence
  predictions$cj_raw <- NA
  predictions[predictions$resp==1,]$cj_raw <- hmvec_up[(predictions[predictions$resp==1,]$closest_evdnc2-1)*timesteps+round(predictions[predictions$resp==1,]$temprt2)]
  predictions[predictions$resp==-1,]$cj_raw <- hmvec_low[(predictions[predictions$resp==-1,]$closest_evdnc2-1)*timesteps+round(predictions[predictions$resp==-1,]$temprt2)]
  
  params <- c(params,drift)
  if (any(is.na(predictions$cj_raw))) {
    save(params,file=paste0("parameters_hm_",unique(observations$sub),unique(observations[,condition_name]),".Rdata"))
    save(predictions,file=paste0("predictions_hm_",unique(observations$sub),unique(observations[,condition_name]),".Rdata"))
    stop(paste(unique(observations$sub),unique(observations[,condition_name]),"hm"))
  }
  
  #add a bias, but keep values bounded between 0 and 1
  predictions$cj_raw <- predictions$cj_raw + params['bias']
  predictions$cj_raw[predictions$cj_raw>1] <- 1
  predictions$cj_raw[predictions$cj_raw<0] <- 0
  
  if (any(is.na(predictions$cj_raw))) {
    save(params,file=paste0("parameters_bias_",unique(observations$sub),unique(observations[,condition_name]),".Rdata"))
    save(predictions,file=paste0("predictions_bias_",unique(observations$sub),unique(observations[,condition_name]),".Rdata"))
    stop(paste(unique(observations$sub),unique(observations[,condition_name]),"bias"))
  }
  
  # Adjust the cut function to match the scale of your confidence report
  if (binning) {
    predictions$cj <- as.numeric(cut(predictions$cj_raw,breaks=seq(0,1,length.out = 7),include.lowest = TRUE))  
  }else{
    predictions$cj <- predictions$cj_raw
    # Scale confidence values to the range of the reports
    predictions[predictions$cj<conf_min,"cj"] <- conf_min
  }
  
  if (any(is.na(predictions$cj))) {
    save(params,file=paste0("parameters_binning_",unique(observations$sub),unique(observations[,condition_name]),".Rdata"))
    save(predictions,file=paste0("predictions_binning_",unique(observations$sub),unique(observations[,condition_name]),".Rdata"))
    stop(paste(unique(observations$sub),unique(observations[,condition_name]),"binning"))
  }
  
  #if we're only simulating data, return the predictions
  if(returnFit==0){
    return(predictions)
    
    #If we are fitting the model, now compare these predictions to the observations
  }else{
    # again, separate the predections according to the response
    c_predicted <- predictions[predictions$cor == 1,]
    e_predicted <- predictions[predictions$cor == 0,]
    
    # First, separate the data in correct and error trials
    c_observed <- observations[observations$cor == 1,]
    e_observed <- observations[observations$cor == 0,]
    
    obs_props <- NULL; pred_props <- NULL;obs_props_cj <- NULL; pred_props_cj <- NULL
    for (d in 1:length(drift)) {
      # Now, get the quantile RTs on the "observed data" for correct and error distributions separately (for quantiles .1, .3, .5, .7, .9)
      c_quantiles <- quantile(c_observed[c_observed$coh == coherences[d],]$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
      e_quantiles <- quantile(e_observed[e_observed$coh == coherences[d],]$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
      if (any(is.na(e_quantiles))) {
        e_quantiles <- rep(0,5)
      }
      if (any(is.na(c_quantiles))) {
        c_quantiles <- rep(0,5)
      }
      # to combine correct and incorrect we scale the expected interquantile probability by the proportion of correct and incorect respectively
      prop_obs_c <- dim(c_observed[c_observed$coh == coherences[d],])[1] / dim(observations)[1]
      prop_obs_e <- dim(e_observed[e_observed$coh == coherences[d],])[1] / dim(observations)[1]
      
      c_obs_proportion = prop_obs_c * c(.1, .2, .2, .2, .2, .1)
      e_obs_proportion = prop_obs_e * c(.1, .2, .2, .2, .2, .1)
      obs_props <- c(obs_props,c_obs_proportion,e_obs_proportion)
      
      c_predicted_rt <- sort(c_predicted[c_predicted$drift == drift[d],]$rt)
      e_predicted_rt <- sort(e_predicted[e_predicted$drift == drift[d],]$rt)
      # now, get the proportion of responses that fall between the observed quantiles when applied to the predicted data
      c_pred_proportion <- c(
        sum(c_predicted_rt <= c_quantiles[1]),
        sum(c_predicted_rt <= c_quantiles[2]) - sum(c_predicted_rt <= c_quantiles[1]),
        sum(c_predicted_rt <= c_quantiles[3]) - sum(c_predicted_rt <= c_quantiles[2]),
        sum(c_predicted_rt <= c_quantiles[4]) - sum(c_predicted_rt <= c_quantiles[3]),
        sum(c_predicted_rt <= c_quantiles[5]) - sum(c_predicted_rt <= c_quantiles[4]),
        sum(c_predicted_rt > c_quantiles[5])
      ) / dim(predictions)[1]
      
      e_pred_proportion <- c(
        sum(e_predicted_rt <= e_quantiles[1]),
        sum(e_predicted_rt <= e_quantiles[2]) - sum(e_predicted_rt <= e_quantiles[1]),
        sum(e_predicted_rt <= e_quantiles[3]) - sum(e_predicted_rt <= e_quantiles[2]),
        sum(e_predicted_rt <= e_quantiles[4]) - sum(e_predicted_rt <= e_quantiles[3]),
        sum(e_predicted_rt <= e_quantiles[5]) - sum(e_predicted_rt <= e_quantiles[4]),
        sum(e_predicted_rt > e_quantiles[5])
      ) / dim(predictions)[1]
      pred_props <- c(pred_props,c_pred_proportion,e_pred_proportion)
      
      
      # Now, do the same for confidence
      c_predicted_cj <- c_predicted[c_predicted$drift == drift[d],]$cj
      e_predicted_cj <- e_predicted[e_predicted$drift == drift[d],]$cj  
      
      if (binning) {
        c_obs_proportion_cj <- data.frame(var1=1:6,Freq=0) #Change according to the scale of your reports
        e_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
        
        # Change cj to your affect column to fit affect
        c_props_cj <- as.data.frame(table(c_observed[c_observed$coh == coherences[d],]$cj)/dim(observations)[1])
        e_props_cj <- as.data.frame(table(e_observed[e_observed$coh == coherences[d],]$cj)/dim(observations)[1])
        
        c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] <- c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] + c_props_cj$Freq
        e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] <- e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] + e_props_cj$Freq
        obs_props_cj <- c(obs_props_cj,c_obs_proportion_cj$Freq,e_obs_proportion_cj$Freq)
        
        c_pred_proportion_cj <- c( #Change according to the scale of your reports
          sum(c_predicted_cj == 1),
          sum(c_predicted_cj == 2),
          sum(c_predicted_cj == 3),
          sum(c_predicted_cj == 4),
          sum(c_predicted_cj == 5),
          sum(c_predicted_cj == 6)
        ) / dim(predictions)[1]
        
        e_pred_proportion_cj <- c( #Change according to the scale of your reports
          sum(e_predicted_cj == 1),
          sum(e_predicted_cj == 2),
          sum(e_predicted_cj == 3),
          sum(e_predicted_cj == 4),
          sum(e_predicted_cj == 5),
          sum(e_predicted_cj == 6)
        ) / dim(predictions)[1]
        
      }else{
        c_quantiles_cj <- quantile(c_observed[c_observed$coh == coherences[d],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
        e_quantiles_cj <- quantile(e_observed[e_observed$coh == coherences[d],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
        if (any(is.na(e_quantiles_cj))) {
          e_quantiles_cj <- rep(0,5)
        }
        if (any(is.na(c_quantiles_cj))) {
          c_quantiles_cj <- rep(0,5)
        }  
        c_pred_proportion_cj <- c(
          sum(c_predicted_cj <= c_quantiles_cj[1]),
          sum(c_predicted_cj <= c_quantiles_cj[2]) - sum(c_predicted_cj <= c_quantiles_cj[1]),
          sum(c_predicted_cj <= c_quantiles_cj[3]) - sum(c_predicted_cj <= c_quantiles_cj[2]),
          sum(c_predicted_cj <= c_quantiles_cj[4]) - sum(c_predicted_cj <= c_quantiles_cj[3]),
          sum(c_predicted_cj <= c_quantiles_cj[5]) - sum(c_predicted_cj <= c_quantiles_cj[4]),
          sum(c_predicted_cj > c_quantiles_cj[5])
        ) / dim(predictions)[1]
        
        e_pred_proportion_cj <- c(
          sum(e_predicted_cj <= e_quantiles_cj[1]),
          sum(e_predicted_cj <= e_quantiles_cj[2]) - sum(e_predicted_cj <= e_quantiles_cj[1]),
          sum(e_predicted_cj <= e_quantiles_cj[3]) - sum(e_predicted_cj <= e_quantiles_cj[2]),
          sum(e_predicted_cj <= e_quantiles_cj[4]) - sum(e_predicted_cj <= e_quantiles_cj[3]),
          sum(e_predicted_cj <= e_quantiles_cj[5]) - sum(e_predicted_cj <= e_quantiles_cj[4]),
          sum(e_predicted_cj > e_quantiles_cj[5])
        ) / dim(predictions)[1]
      }
      
      
      pred_props_cj <- c(pred_props_cj,c_pred_proportion_cj,e_pred_proportion_cj)
    }    
    
    # Combine the quantiles for rts and cj
    if (binning) {obs_props <- c(obs_props,obs_props_cj)
    }else{obs_props <- c(obs_props,obs_props)}
    
    pred_props <- c(pred_props,pred_props_cj)
    # calculate chi square
    
    # calculate chi square
    chiSquare = sum( (obs_props - pred_props) ^ 2)
    return(chiSquare)
  }
}
quantile_optim_DDM_Vs_biasfixed <- function(params, observations, returnFit,ddm_params, confRT_name = "RTconf",condition_name = "selfconf",
                                            ev_bound = .5, ev_window = .01, upperRT = 5, conf_min = 0, binning = T){ # Heatmap parameters
  #First, generate predictions:
  v_ratio <- params[6:8]
  # v_s = subjective drift/prior belief parameter
  v_s <- params[9:11]
  params <- params[1:5]
  names(params) <- c('z','ntrials','sigma','dt','bias')
  
  condition = sort(unique(observations[,condition_name]))
  coherences <- sort(unique(observations$coh)) #/!\ CHANGE ACCORDING TO DATASET
  
  # Generate trials from DDM parameters
  for (cond in 1:length(condition)) {
    trial = data.frame(matrix(NA,nrow=0,ncol=9))
    names(trial) <- c('rt','resp','cor','evidence2','rt2','cj','drift','condition','difflevel')
    for (d in 1:length(coherences)) {
      index_ddm <- ddm_params$condition==condition[cond] & ddm_params$difflevel == coherences[d]
      predictions <- data.frame(DDM_with_confidence_slow_fullconfRT(
        v=ddm_params[index_ddm,"drift"],a=ddm_params[index_ddm,"bound"],ter=ddm_params[index_ddm,"ter"],z=params['z'],
        ntrials=params['ntrials']*dim(observations)[1]/2/length(coherences)/length(condition),s=params['sigma'],
        dt=params['dt'],t2distribution=rep(observations[observations$coh==coherences[d]&observations[,condition_name]==condition[cond],confRT_name],times=params['ntrials']),
        postdriftmod=v_ratio[cond]))
      predictionsneg <- data.frame(DDM_with_confidence_slow_fullconfRT(
        v=-ddm_params[index_ddm,"drift"],a=ddm_params[index_ddm,"bound"],ter=ddm_params[index_ddm,"ter"],z=params['z'],
        ntrials=params['ntrials']*dim(observations)[1]/2/length(coherences)/length(condition),s=params['sigma'],
        dt=params['dt'],t2distribution=rep(observations[observations$coh==coherences[d]&observations[,condition_name]==condition[cond],confRT_name],times=params['ntrials']),
        postdriftmod=v_ratio[cond]))
      names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
      names(predictionsneg) <- c('rt','resp','cor','evidence2','rt2','cj')
      predictions <- fastmerge(predictions,predictionsneg)
      predictions['drift'] <- ddm_params[index_ddm,"drift"]
      predictions$condition <- condition[cond]
      predictions$difflevel <- coherences[d]
      trial <- fastmerge(trial,predictions)
    }
    predictions <- trial
    
    # Generate the heatmap from v_s
    hm_up <- fast_hm(v_s[cond],sigma = params["sigma"], dt = params["dt"], 
                     ev_bound = ev_bound, ev_window = ev_window, upperRT = upperRT)
    hm_low <- 1-hm_up
    hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
    
    # match data to the heatmap
    timesteps <- upperRT/params["dt"]
    ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
    
    
    predictions$closest_evdnc2 <- match.closest(predictions$evidence2,ev_mapping)
    predictions$temprt2 <- predictions$rt2;
    predictions$temprt2[predictions$temprt2>upperRT] <- upperRT #heatmap doesn't go higher
    predictions$temprt2 <- predictions$temprt2*timesteps/upperRT #scale with the heatmap
    
    # Compute confidence
    predictions$cj_raw <- NA
    predictions[predictions$resp==1,]$cj_raw <- hmvec_up[(predictions[predictions$resp==1,]$closest_evdnc2-1)*timesteps+round(predictions[predictions$resp==1,]$temprt2)]
    predictions[predictions$resp==-1,]$cj_raw <- hmvec_low[(predictions[predictions$resp==-1,]$closest_evdnc2-1)*timesteps+round(predictions[predictions$resp==-1,]$temprt2)]
    
    #add a bias, but keep values bounded between 0 and 1
    predictions$cj_raw <- predictions$cj_raw + params['bias']
    predictions$cj_raw[predictions$cj_raw>1] <- 1
    predictions$cj_raw[predictions$cj_raw<0] <- 0
    
    # Adjust the cut function to match the scale of your confidence report
    if (binning) {
      predictions$cj <- as.numeric(cut(predictions$cj_raw,breaks=seq(0,1,length.out = 7),include.lowest = TRUE))  
    }else{
      predictions$cj <- predictions$cj_raw
      # Scale confidence values to the range of the reports
      predictions[predictions$cj<conf_min,"cj"] <- conf_min
    }
    if (cond==1) {
      predictions_tot <- predictions
    }else{
      predictions_tot <- fastmerge(predictions_tot,predictions)
    }
  }
  predictions <- predictions_tot  
  
  
  #if we're only simulating data, return the predictions
  if(returnFit==0){
    return(predictions)
    
    #If we are fitting the model, now compare these predictions to the observations
  }else{
    
    # again, separate the predections according to the response
    c_predicted <- predictions[predictions$cor == 1,]
    e_predicted <- predictions[predictions$cor == 0,]
    
    # First, separate the data in correct and error trials
    c_observed <- observations[observations$cor == 1,]
    e_observed <- observations[observations$cor == 0,]
    
    obs_props <- NULL; pred_props <- NULL;obs_props_cj <- NULL; pred_props_cj <- NULL
    for (cond in 1:length(condition)) {
      for (d in 1:length(coherences)) {
        # Now, do the same for confidence
        c_predicted_cj <- c_predicted[c_predicted$condition == condition[cond] & c_predicted$difflevel == coherences[d],]$cj
        e_predicted_cj <- e_predicted[e_predicted$condition == condition[cond] & e_predicted$difflevel == coherences[d],]$cj  
        
        if (binning) {
          c_obs_proportion_cj <- data.frame(var1=1:6,Freq=0) #Change according to the scale of your reports
          e_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
          
          # Change cj to your affect column to fit affect
          c_props_cj <- as.data.frame(table(c_observed[c_observed[,condition_name]==condition[cond] & c_observed$coh == coherences[d],]$cj)/dim(observations)[1])
          e_props_cj <- as.data.frame(table(e_observed[e_observed[,condition_name]==condition[cond] & e_observed$coh == coherences[d],]$cj)/dim(observations)[1])
          
          c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] <- c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] + c_props_cj$Freq
          e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] <- e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] + e_props_cj$Freq
          obs_props_cj <- c(obs_props_cj,c_obs_proportion_cj$Freq,e_obs_proportion_cj$Freq)
          
          c_pred_proportion_cj <- c( #Change according to the scale of your reports
            sum(c_predicted_cj == 1),
            sum(c_predicted_cj == 2),
            sum(c_predicted_cj == 3),
            sum(c_predicted_cj == 4),
            sum(c_predicted_cj == 5),
            sum(c_predicted_cj == 6)
          ) / dim(predictions)[1]
          
          e_pred_proportion_cj <- c( #Change according to the scale of your reports
            sum(e_predicted_cj == 1),
            sum(e_predicted_cj == 2),
            sum(e_predicted_cj == 3),
            sum(e_predicted_cj == 4),
            sum(e_predicted_cj == 5),
            sum(e_predicted_cj == 6)
          ) / dim(predictions)[1]
          
        }else{
          c_quantiles_cj <- quantile(c_observed[c_observed[,condition_name]==condition[cond] & c_observed$coh == coherences[d],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
          e_quantiles_cj <- quantile(e_observed[e_observed[,condition_name]==condition[cond] & e_observed$coh == coherences[d],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
          if (any(is.na(e_quantiles_cj))) {
            e_quantiles_cj <- rep(0,5)
          }
          if (any(is.na(c_quantiles_cj))) {
            c_quantiles_cj <- rep(0,5)
          }  
          c_pred_proportion_cj <- c(
            sum(c_predicted_cj <= c_quantiles_cj[1]),
            sum(c_predicted_cj <= c_quantiles_cj[2]) - sum(c_predicted_cj <= c_quantiles_cj[1]),
            sum(c_predicted_cj <= c_quantiles_cj[3]) - sum(c_predicted_cj <= c_quantiles_cj[2]),
            sum(c_predicted_cj <= c_quantiles_cj[4]) - sum(c_predicted_cj <= c_quantiles_cj[3]),
            sum(c_predicted_cj <= c_quantiles_cj[5]) - sum(c_predicted_cj <= c_quantiles_cj[4]),
            sum(c_predicted_cj > c_quantiles_cj[5])
          ) / dim(predictions)[1]
          
          e_pred_proportion_cj <- c(
            sum(e_predicted_cj <= e_quantiles_cj[1]),
            sum(e_predicted_cj <= e_quantiles_cj[2]) - sum(e_predicted_cj <= e_quantiles_cj[1]),
            sum(e_predicted_cj <= e_quantiles_cj[3]) - sum(e_predicted_cj <= e_quantiles_cj[2]),
            sum(e_predicted_cj <= e_quantiles_cj[4]) - sum(e_predicted_cj <= e_quantiles_cj[3]),
            sum(e_predicted_cj <= e_quantiles_cj[5]) - sum(e_predicted_cj <= e_quantiles_cj[4]),
            sum(e_predicted_cj > e_quantiles_cj[5])
          ) / dim(predictions)[1]
        }
        
        
        pred_props_cj <- c(pred_props_cj,c_pred_proportion_cj,e_pred_proportion_cj)
      }    
      
    }
    
    # Combine the quantiles for rts and cj
    if (binning) {obs_props <- c(obs_props,obs_props_cj)
    }else{obs_props <- c(obs_props,obs_props)}
    
    pred_props <- c(pred_props,pred_props_cj)
    # calculate chi square
    
    # calculate chi square
    chiSquare = sum( (obs_props - pred_props) ^ 2)
    return(chiSquare)
  }
}
quantile_optim_DDM_Vs_bias_noddm <- function(params, observations, returnFit,ddm_params, confRT_name = "RTconf",condition_name = "selfconf",
                                             ev_bound = .5, ev_window = .01, upperRT = 5, conf_min = 0, binning = T,ddm="condition"){ # Heatmap parameters
  #First, generate predictions:
  bias <- params[5:7]
  v_ratio <- params[8:10]
  # v_s = subjective drift/prior belief parameter
  v_s <- params[11:13]
  params <- params[1:4]
  names(params) <- c('z','ntrials','sigma','dt')
  
  condition = sort(unique(observations[,condition_name]))
  coherences <- sort(unique(observations$coh)) #/!\ CHANGE ACCORDING TO DATASET
  
  # Generate trials from DDM parameters
  for (cond in 1:length(condition)) {
    trial = data.frame(matrix(NA,nrow=0,ncol=9))
    names(trial) <- c('rt','resp','cor','evidence2','rt2','cj','drift','condition','difflevel')
    for (d in 1:length(coherences)) {
      if (ddm=="condition") {
        index_ddm <- ddm_params$condition==condition[cond] & ddm_params$difflevel == coherences[d]
      } else {
        index_ddm <- ddm_params$difflevel == coherences[d]
      }
      predictions <- data.frame(DDM_with_confidence_slow_fullconfRT(
        v=ddm_params[index_ddm,"drift"],a=ddm_params[index_ddm,"bound"],ter=ddm_params[index_ddm,"ter"],z=params['z'],
        ntrials=params['ntrials']*dim(observations)[1]/2/length(coherences)/length(condition),s=params['sigma'],
        dt=params['dt'],t2distribution=rep(observations[observations$coh==coherences[d]&observations[,condition_name]==condition[cond],confRT_name],times=params['ntrials']),
        postdriftmod=v_ratio[cond]))
      predictionsneg <- data.frame(DDM_with_confidence_slow_fullconfRT(
        v=-ddm_params[index_ddm,"drift"],a=ddm_params[index_ddm,"bound"],ter=ddm_params[index_ddm,"ter"],z=params['z'],
        ntrials=params['ntrials']*dim(observations)[1]/2/length(coherences)/length(condition),s=params['sigma'],
        dt=params['dt'],t2distribution=rep(observations[observations$coh==coherences[d]&observations[,condition_name]==condition[cond],confRT_name],times=params['ntrials']),
        postdriftmod=v_ratio[cond]))
      names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
      names(predictionsneg) <- c('rt','resp','cor','evidence2','rt2','cj')
      predictions <- fastmerge(predictions,predictionsneg)
      predictions['drift'] <- ddm_params[index_ddm,"drift"]
      predictions$condition <- condition[cond]
      predictions$difflevel <- coherences[d]
      trial <- fastmerge(trial,predictions)
    }
    predictions <- trial
    
    # Generate the heatmap from v_s
    hm_up <- fast_hm(v_s[cond],sigma = params["sigma"], dt = params["dt"], 
                     ev_bound = ev_bound, ev_window = ev_window, upperRT = upperRT)
    hm_low <- 1-hm_up
    hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
    
    # match data to the heatmap
    timesteps <- upperRT/params["dt"]
    ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
    
    
    predictions$closest_evdnc2 <- match.closest(predictions$evidence2,ev_mapping)
    predictions$temprt2 <- predictions$rt2;
    predictions$temprt2[predictions$temprt2>upperRT] <- upperRT #heatmap doesn't go higher
    predictions$temprt2 <- predictions$temprt2*timesteps/upperRT #scale with the heatmap
    
    # Compute confidence
    predictions$cj_raw <- NA
    predictions[predictions$resp==1,]$cj_raw <- hmvec_up[(predictions[predictions$resp==1,]$closest_evdnc2-1)*timesteps+round(predictions[predictions$resp==1,]$temprt2)]
    predictions[predictions$resp==-1,]$cj_raw <- hmvec_low[(predictions[predictions$resp==-1,]$closest_evdnc2-1)*timesteps+round(predictions[predictions$resp==-1,]$temprt2)]
    
    #add a bias, but keep values bounded between 0 and 1
    predictions$cj_raw <- predictions$cj_raw + bias[cond]
    predictions$cj_raw[predictions$cj_raw>1] <- 1
    predictions$cj_raw[predictions$cj_raw<0] <- 0
    
    # Adjust the cut function to match the scale of your confidence report
    if (binning) {
      predictions$cj <- as.numeric(cut(predictions$cj_raw,breaks=seq(0,1,length.out = 7),include.lowest = TRUE))  
    }else{
      predictions$cj <- predictions$cj_raw
      # Scale confidence values to the range of the reports
      predictions[predictions$cj<conf_min,"cj"] <- conf_min
    }
    if (cond==1) {
      predictions_tot <- predictions
    }else{
      predictions_tot <- fastmerge(predictions_tot,predictions)
    }
  }
  predictions <- predictions_tot  
  
  
  #if we're only simulating data, return the predictions
  if(returnFit==0){
    return(predictions)
    
    #If we are fitting the model, now compare these predictions to the observations
  }else{
    
    # again, separate the predections according to the response
    c_predicted <- predictions[predictions$cor == 1,]
    e_predicted <- predictions[predictions$cor == 0,]
    
    # First, separate the data in correct and error trials
    c_observed <- observations[observations$cor == 1,]
    e_observed <- observations[observations$cor == 0,]
    
    obs_props <- NULL; pred_props <- NULL;obs_props_cj <- NULL; pred_props_cj <- NULL
    for (cond in 1:length(condition)) {
      for (d in 1:length(coherences)) {
        # Now, do the same for confidence
        c_predicted_cj <- c_predicted[c_predicted$condition == condition[cond] & c_predicted$difflevel == coherences[d],]$cj
        e_predicted_cj <- e_predicted[e_predicted$condition == condition[cond] & e_predicted$difflevel == coherences[d],]$cj  
        
        if (binning) {
          c_obs_proportion_cj <- data.frame(var1=1:6,Freq=0) #Change according to the scale of your reports
          e_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
          
          # Change cj to your affect column to fit affect
          c_props_cj <- as.data.frame(table(c_observed[c_observed[,condition_name]==condition[cond] & c_observed$coh == coherences[d],]$cj)/dim(observations)[1])
          e_props_cj <- as.data.frame(table(e_observed[e_observed[,condition_name]==condition[cond] & e_observed$coh == coherences[d],]$cj)/dim(observations)[1])
          
          c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] <- c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] + c_props_cj$Freq
          e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] <- e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] + e_props_cj$Freq
          obs_props_cj <- c(obs_props_cj,c_obs_proportion_cj$Freq,e_obs_proportion_cj$Freq)
          
          c_pred_proportion_cj <- c( #Change according to the scale of your reports
            sum(c_predicted_cj == 1),
            sum(c_predicted_cj == 2),
            sum(c_predicted_cj == 3),
            sum(c_predicted_cj == 4),
            sum(c_predicted_cj == 5),
            sum(c_predicted_cj == 6)
          ) / dim(predictions)[1]
          
          e_pred_proportion_cj <- c( #Change according to the scale of your reports
            sum(e_predicted_cj == 1),
            sum(e_predicted_cj == 2),
            sum(e_predicted_cj == 3),
            sum(e_predicted_cj == 4),
            sum(e_predicted_cj == 5),
            sum(e_predicted_cj == 6)
          ) / dim(predictions)[1]
          
        }else{
          c_quantiles_cj <- quantile(c_observed[c_observed[,condition_name]==condition[cond] & c_observed$coh == coherences[d],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
          e_quantiles_cj <- quantile(e_observed[e_observed[,condition_name]==condition[cond] & e_observed$coh == coherences[d],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
          if (any(is.na(e_quantiles_cj))) {
            e_quantiles_cj <- rep(0,5)
          }
          if (any(is.na(c_quantiles_cj))) {
            c_quantiles_cj <- rep(0,5)
          }  
          c_pred_proportion_cj <- c(
            sum(c_predicted_cj <= c_quantiles_cj[1]),
            sum(c_predicted_cj <= c_quantiles_cj[2]) - sum(c_predicted_cj <= c_quantiles_cj[1]),
            sum(c_predicted_cj <= c_quantiles_cj[3]) - sum(c_predicted_cj <= c_quantiles_cj[2]),
            sum(c_predicted_cj <= c_quantiles_cj[4]) - sum(c_predicted_cj <= c_quantiles_cj[3]),
            sum(c_predicted_cj <= c_quantiles_cj[5]) - sum(c_predicted_cj <= c_quantiles_cj[4]),
            sum(c_predicted_cj > c_quantiles_cj[5])
          ) / dim(predictions)[1]
          
          e_pred_proportion_cj <- c(
            sum(e_predicted_cj <= e_quantiles_cj[1]),
            sum(e_predicted_cj <= e_quantiles_cj[2]) - sum(e_predicted_cj <= e_quantiles_cj[1]),
            sum(e_predicted_cj <= e_quantiles_cj[3]) - sum(e_predicted_cj <= e_quantiles_cj[2]),
            sum(e_predicted_cj <= e_quantiles_cj[4]) - sum(e_predicted_cj <= e_quantiles_cj[3]),
            sum(e_predicted_cj <= e_quantiles_cj[5]) - sum(e_predicted_cj <= e_quantiles_cj[4]),
            sum(e_predicted_cj > e_quantiles_cj[5])
          ) / dim(predictions)[1]
        }
        
        
        pred_props_cj <- c(pred_props_cj,c_pred_proportion_cj,e_pred_proportion_cj)
      }    
      
    }
    
    # Combine the quantiles for rts and cj
    if (binning) {obs_props <- c(obs_props,obs_props_cj)
    }else{obs_props <- c(obs_props,obs_props)}
    
    pred_props <- c(pred_props,pred_props_cj)
    # calculate chi square
    
    # calculate chi square
    chiSquare = sum( (obs_props - pred_props) ^ 2)
    return(chiSquare)
  }
}
quantile_optim_DDM_Vs_biasfixed_ddmfixed <- function(params, observations, returnFit,ddm_params, confRT_name = "RTconf",condition_name = "selfconf",
                                                     ev_bound = .5, ev_window = .01, upperRT = 5, conf_min = 0, binning = T){ # Heatmap parameters
  #First, generate predictions:
  v_ratio <- params[6:8]
  # v_s = subjective drift/prior belief parameter
  v_s <- params[9:11]
  params <- params[1:5]
  names(params) <- c('z','ntrials','sigma','dt','bias')
  
  condition = sort(unique(observations[,condition_name])) # Sort used to get the same order for all participants
  coherences <- sort(unique(observations$coh)) #/!\ CHANGE ACCORDING TO DATASET
  
  # Generate trials from DDM parameters
  for (cond in 1:length(condition)) {
    trial = data.frame(matrix(NA,nrow=0,ncol=9))
    names(trial) <- c('rt','resp','cor','evidence2','rt2','cj','drift','condition','difflevel')
    for (d in 1:length(coherences)) {
      index_ddm <- ddm_params$difflevel == coherences[d]
      predictions <- data.frame(DDM_with_confidence_slow_fullconfRT(
        v=ddm_params[index_ddm,"drift"],a=ddm_params[index_ddm,"bound"],ter=ddm_params[index_ddm,"ter"],z=params['z'],
        ntrials=params['ntrials']*dim(observations)[1]/2/length(coherences)/length(condition),s=params['sigma'],
        dt=params['dt'],t2distribution=rep(observations[observations$coh==coherences[d]&observations[,condition_name]==condition[cond],confRT_name],times=params['ntrials']),
        postdriftmod=v_ratio[cond]))
      predictionsneg <- data.frame(DDM_with_confidence_slow_fullconfRT(
        v=-ddm_params[index_ddm,"drift"],a=ddm_params[index_ddm,"bound"],ter=ddm_params[index_ddm,"ter"],z=params['z'],
        ntrials=params['ntrials']*dim(observations)[1]/2/length(coherences)/length(condition),s=params['sigma'],
        dt=params['dt'],t2distribution=rep(observations[observations$coh==coherences[d]&observations[,condition_name]==condition[cond],confRT_name],times=params['ntrials']),
        postdriftmod=v_ratio[cond]))
      names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
      names(predictionsneg) <- c('rt','resp','cor','evidence2','rt2','cj')
      predictions <- fastmerge(predictions,predictionsneg)
      predictions['drift'] <- ddm_params[index_ddm,"drift"]
      predictions$condition <- condition[cond]
      predictions$difflevel <- coherences[d]
      trial <- fastmerge(trial,predictions)
    }
    predictions <- trial
    
    # Generate the heatmap from v_s
    hm_up <- fast_hm(v_s[cond],sigma = params["sigma"], dt = params["dt"], 
                     ev_bound = ev_bound, ev_window = ev_window, upperRT = upperRT)
    hm_low <- 1-hm_up
    hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
    
    # match data to the heatmap
    timesteps <- upperRT/params["dt"]
    ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
    
    
    predictions$closest_evdnc2 <- match.closest(predictions$evidence2,ev_mapping)
    predictions$temprt2 <- predictions$rt2;
    predictions$temprt2[predictions$temprt2>upperRT] <- upperRT #heatmap doesn't go higher
    predictions$temprt2 <- predictions$temprt2*timesteps/upperRT #scale with the heatmap
    
    # Compute confidence
    predictions$cj_raw <- NA
    predictions[predictions$resp==1,]$cj_raw <- hmvec_up[(predictions[predictions$resp==1,]$closest_evdnc2-1)*timesteps+round(predictions[predictions$resp==1,]$temprt2)]
    predictions[predictions$resp==-1,]$cj_raw <- hmvec_low[(predictions[predictions$resp==-1,]$closest_evdnc2-1)*timesteps+round(predictions[predictions$resp==-1,]$temprt2)]
    
    #add a bias, but keep values bounded between 0 and 1
    predictions$cj_raw <- predictions$cj_raw + params['bias']
    predictions$cj_raw[predictions$cj_raw>1] <- 1
    predictions$cj_raw[predictions$cj_raw<0] <- 0
    
    # Adjust the cut function to match the scale of your confidence report
    if (binning) {
      predictions$cj <- as.numeric(cut(predictions$cj_raw,breaks=seq(0,1,length.out = 7),include.lowest = TRUE))  
    }else{
      predictions$cj <- predictions$cj_raw
      # Scale confidence values to the range of the reports
      predictions[predictions$cj<conf_min,"cj"] <- conf_min
    }
    if (cond==1) {
      predictions_tot <- predictions
    }else{
      predictions_tot <- fastmerge(predictions_tot,predictions)
    }
  }
  predictions <- predictions_tot  
  
  
  #if we're only simulating data, return the predictions
  if(returnFit==0){
    return(predictions)
    
    #If we are fitting the model, now compare these predictions to the observations
  }else{
    
    # again, separate the predections according to the response
    c_predicted <- predictions[predictions$cor == 1,]
    e_predicted <- predictions[predictions$cor == 0,]
    
    # First, separate the data in correct and error trials
    c_observed <- observations[observations$cor == 1,]
    e_observed <- observations[observations$cor == 0,]
    
    obs_props <- NULL; pred_props <- NULL;obs_props_cj <- NULL; pred_props_cj <- NULL
    for (cond in 1:length(condition)) {
      for (d in 1:length(coherences)) {
        # Now, do the same for confidence
        c_predicted_cj <- c_predicted[c_predicted$condition == condition[cond] & c_predicted$difflevel == coherences[d],]$cj
        e_predicted_cj <- e_predicted[e_predicted$condition == condition[cond] & e_predicted$difflevel == coherences[d],]$cj  
        
        if (binning) {
          c_obs_proportion_cj <- data.frame(var1=1:6,Freq=0) #Change according to the scale of your reports
          e_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
          
          # Change cj to your affect column to fit affect
          c_props_cj <- as.data.frame(table(c_observed[c_observed[,condition_name]==condition[cond] & c_observed$coh == coherences[d],]$cj)/dim(observations)[1])
          e_props_cj <- as.data.frame(table(e_observed[e_observed[,condition_name]==condition[cond] & e_observed$coh == coherences[d],]$cj)/dim(observations)[1])
          
          c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] <- c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] + c_props_cj$Freq
          e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] <- e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] + e_props_cj$Freq
          obs_props_cj <- c(obs_props_cj,c_obs_proportion_cj$Freq,e_obs_proportion_cj$Freq)
          
          c_pred_proportion_cj <- c( #Change according to the scale of your reports
            sum(c_predicted_cj == 1),
            sum(c_predicted_cj == 2),
            sum(c_predicted_cj == 3),
            sum(c_predicted_cj == 4),
            sum(c_predicted_cj == 5),
            sum(c_predicted_cj == 6)
          ) / dim(predictions)[1]
          
          e_pred_proportion_cj <- c( #Change according to the scale of your reports
            sum(e_predicted_cj == 1),
            sum(e_predicted_cj == 2),
            sum(e_predicted_cj == 3),
            sum(e_predicted_cj == 4),
            sum(e_predicted_cj == 5),
            sum(e_predicted_cj == 6)
          ) / dim(predictions)[1]
          
        }else{
          c_quantiles_cj <- quantile(c_observed[c_observed[,condition_name]==condition[cond] & c_observed$coh == coherences[d],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
          e_quantiles_cj <- quantile(e_observed[e_observed[,condition_name]==condition[cond] & e_observed$coh == coherences[d],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
          if (any(is.na(e_quantiles_cj))) {
            e_quantiles_cj <- rep(0,5)
          }
          if (any(is.na(c_quantiles_cj))) {
            c_quantiles_cj <- rep(0,5)
          }  
          c_pred_proportion_cj <- c(
            sum(c_predicted_cj <= c_quantiles_cj[1]),
            sum(c_predicted_cj <= c_quantiles_cj[2]) - sum(c_predicted_cj <= c_quantiles_cj[1]),
            sum(c_predicted_cj <= c_quantiles_cj[3]) - sum(c_predicted_cj <= c_quantiles_cj[2]),
            sum(c_predicted_cj <= c_quantiles_cj[4]) - sum(c_predicted_cj <= c_quantiles_cj[3]),
            sum(c_predicted_cj <= c_quantiles_cj[5]) - sum(c_predicted_cj <= c_quantiles_cj[4]),
            sum(c_predicted_cj > c_quantiles_cj[5])
          ) / dim(predictions)[1]
          
          e_pred_proportion_cj <- c(
            sum(e_predicted_cj <= e_quantiles_cj[1]),
            sum(e_predicted_cj <= e_quantiles_cj[2]) - sum(e_predicted_cj <= e_quantiles_cj[1]),
            sum(e_predicted_cj <= e_quantiles_cj[3]) - sum(e_predicted_cj <= e_quantiles_cj[2]),
            sum(e_predicted_cj <= e_quantiles_cj[4]) - sum(e_predicted_cj <= e_quantiles_cj[3]),
            sum(e_predicted_cj <= e_quantiles_cj[5]) - sum(e_predicted_cj <= e_quantiles_cj[4]),
            sum(e_predicted_cj > e_quantiles_cj[5])
          ) / dim(predictions)[1]
        }
        
        
        pred_props_cj <- c(pred_props_cj,c_pred_proportion_cj,e_pred_proportion_cj)
      }    
      
    }
    
    # Combine the quantiles for rts and cj
    if (binning) {obs_props <- c(obs_props,obs_props_cj)
    }else{obs_props <- c(obs_props,obs_props)}
    
    pred_props <- c(pred_props,pred_props_cj)
    # calculate chi square
    
    # calculate chi square
    chiSquare = sum( (obs_props - pred_props) ^ 2)
    return(chiSquare)
  }
}
quantile_optim_DDM_Vsfixed_bias <- function(params, observations, returnFit,ddm_params, confRT_name = "RTconf",condition_name = "selfconf",
                                            ev_bound = .5, ev_window = .01, upperRT = 5, conf_min = 0, binning = T,ddm="condition"){ # Heatmap parameters
  #First, generate predictions:
  v_ratio <- params[8:10]
  # v_s = subjective drift/prior belief parameter
  v_s <- params[11]
  bias <- params[5:7]
  params <- params[1:4]
  names(params) <- c('z','ntrials','sigma','dt')
  
  condition = sort(unique(observations[,condition_name]))
  coherences <- sort(unique(observations$coh)) #/!\ CHANGE ACCORDING TO DATASET
  
  # Generate trials from DDM parameters
  for (cond in 1:length(condition)) {
    trial = data.frame(matrix(NA,nrow=0,ncol=9))
    names(trial) <- c('rt','resp','cor','evidence2','rt2','cj','drift','condition','difflevel')
    for (d in 1:length(coherences)) {
      if (ddm=="condition") {
        index_ddm <- ddm_params$condition==condition[cond] & ddm_params$difflevel == coherences[d]
      } else {
        index_ddm <- ddm_params$difflevel == coherences[d]
      }
      predictions <- data.frame(DDM_with_confidence_slow_fullconfRT(
        v=ddm_params[index_ddm,"drift"],a=ddm_params[index_ddm,"bound"],ter=ddm_params[index_ddm,"ter"],z=params['z'],
        ntrials=params['ntrials']*dim(observations)[1]/2/length(coherences)/length(condition),s=params['sigma'],
        dt=params['dt'],t2distribution=rep(observations[observations$coh==coherences[d]&observations[,condition_name]==condition[cond],confRT_name],times=params['ntrials']),
        postdriftmod=v_ratio[cond]))
      predictionsneg <- data.frame(DDM_with_confidence_slow_fullconfRT(
        v=-ddm_params[index_ddm,"drift"],a=ddm_params[index_ddm,"bound"],ter=ddm_params[index_ddm,"ter"],z=params['z'],
        ntrials=params['ntrials']*dim(observations)[1]/2/length(coherences)/length(condition),s=params['sigma'],
        dt=params['dt'],t2distribution=rep(observations[observations$coh==coherences[d]&observations[,condition_name]==condition[cond],confRT_name],times=params['ntrials']),
        postdriftmod=v_ratio[cond]))
      names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
      names(predictionsneg) <- c('rt','resp','cor','evidence2','rt2','cj')
      predictions <- fastmerge(predictions,predictionsneg)
      predictions['drift'] <- ddm_params[index_ddm,"drift"]
      predictions$condition <- condition[cond]
      predictions$difflevel <- coherences[d]
      trial <- fastmerge(trial,predictions)
    }
    predictions <- trial
    
    # Generate the heatmap from v_s
    hm_up <- fast_hm(v_s,sigma = params["sigma"], dt = params["dt"], 
                     ev_bound = ev_bound, ev_window = ev_window, upperRT = upperRT)
    hm_low <- 1-hm_up
    hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
    
    # match data to the heatmap
    timesteps <- upperRT/params["dt"]
    ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
    
    
    predictions$closest_evdnc2 <- match.closest(predictions$evidence2,ev_mapping)
    predictions$temprt2 <- predictions$rt2;
    predictions$temprt2[predictions$temprt2>upperRT] <- upperRT #heatmap doesn't go higher
    predictions$temprt2 <- predictions$temprt2*timesteps/upperRT #scale with the heatmap
    
    # Compute confidence
    predictions$cj_raw <- NA
    predictions[predictions$resp==1,]$cj_raw <- hmvec_up[(predictions[predictions$resp==1,]$closest_evdnc2-1)*timesteps+round(predictions[predictions$resp==1,]$temprt2)]
    predictions[predictions$resp==-1,]$cj_raw <- hmvec_low[(predictions[predictions$resp==-1,]$closest_evdnc2-1)*timesteps+round(predictions[predictions$resp==-1,]$temprt2)]
    
    #add a bias, but keep values bounded between 0 and 1
    predictions$cj_raw <- predictions$cj_raw + bias[cond]
    predictions$cj_raw[predictions$cj_raw>1] <- 1
    predictions$cj_raw[predictions$cj_raw<0] <- 0
    
    # Adjust the cut function to match the scale of your confidence report
    if (binning) {
      predictions$cj <- as.numeric(cut(predictions$cj_raw,breaks=seq(0,1,length.out = 7),include.lowest = TRUE))  
    }else{
      predictions$cj <- predictions$cj_raw
      # Scale confidence values to the range of the reports
      predictions[predictions$cj<conf_min,"cj"] <- conf_min
    }
    if (cond==1) {
      predictions_tot <- predictions
    }else{
      predictions_tot <- fastmerge(predictions_tot,predictions)
    }
  }
  predictions <- predictions_tot  
  
  
  #if we're only simulating data, return the predictions
  if(returnFit==0){
    return(predictions)
    
    #If we are fitting the model, now compare these predictions to the observations
  }else{
    
    # again, separate the predections according to the response
    c_predicted <- predictions[predictions$cor == 1,]
    e_predicted <- predictions[predictions$cor == 0,]
    
    # First, separate the data in correct and error trials
    c_observed <- observations[observations$cor == 1,]
    e_observed <- observations[observations$cor == 0,]
    
    obs_props <- NULL; pred_props <- NULL;obs_props_cj <- NULL; pred_props_cj <- NULL
    for (cond in 1:length(condition)) {
      for (d in 1:length(coherences)) {
        # Now, do the same for confidence
        c_predicted_cj <- c_predicted[c_predicted$condition == condition[cond] & c_predicted$difflevel == coherences[d],]$cj
        e_predicted_cj <- e_predicted[e_predicted$condition == condition[cond] & e_predicted$difflevel == coherences[d],]$cj  
        
        if (binning) {
          c_obs_proportion_cj <- data.frame(var1=1:6,Freq=0) #Change according to the scale of your reports
          e_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
          
          # Change cj to your affect column to fit affect
          c_props_cj <- as.data.frame(table(c_observed[c_observed[,condition_name]==condition[cond] & c_observed$coh == coherences[d],]$cj)/dim(observations)[1])
          e_props_cj <- as.data.frame(table(e_observed[e_observed[,condition_name]==condition[cond] & e_observed$coh == coherences[d],]$cj)/dim(observations)[1])
          
          c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] <- c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] + c_props_cj$Freq
          e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] <- e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] + e_props_cj$Freq
          obs_props_cj <- c(obs_props_cj,c_obs_proportion_cj$Freq,e_obs_proportion_cj$Freq)
          
          c_pred_proportion_cj <- c( #Change according to the scale of your reports
            sum(c_predicted_cj == 1),
            sum(c_predicted_cj == 2),
            sum(c_predicted_cj == 3),
            sum(c_predicted_cj == 4),
            sum(c_predicted_cj == 5),
            sum(c_predicted_cj == 6)
          ) / dim(predictions)[1]
          
          e_pred_proportion_cj <- c( #Change according to the scale of your reports
            sum(e_predicted_cj == 1),
            sum(e_predicted_cj == 2),
            sum(e_predicted_cj == 3),
            sum(e_predicted_cj == 4),
            sum(e_predicted_cj == 5),
            sum(e_predicted_cj == 6)
          ) / dim(predictions)[1]
          
        }else{
          c_quantiles_cj <- quantile(c_observed[c_observed[,condition_name]==condition[cond] & c_observed$coh == coherences[d],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
          e_quantiles_cj <- quantile(e_observed[e_observed[,condition_name]==condition[cond] & e_observed$coh == coherences[d],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
          if (any(is.na(e_quantiles_cj))) {
            e_quantiles_cj <- rep(0,5)
          }
          if (any(is.na(c_quantiles_cj))) {
            c_quantiles_cj <- rep(0,5)
          }  
          c_pred_proportion_cj <- c(
            sum(c_predicted_cj <= c_quantiles_cj[1]),
            sum(c_predicted_cj <= c_quantiles_cj[2]) - sum(c_predicted_cj <= c_quantiles_cj[1]),
            sum(c_predicted_cj <= c_quantiles_cj[3]) - sum(c_predicted_cj <= c_quantiles_cj[2]),
            sum(c_predicted_cj <= c_quantiles_cj[4]) - sum(c_predicted_cj <= c_quantiles_cj[3]),
            sum(c_predicted_cj <= c_quantiles_cj[5]) - sum(c_predicted_cj <= c_quantiles_cj[4]),
            sum(c_predicted_cj > c_quantiles_cj[5])
          ) / dim(predictions)[1]
          
          e_pred_proportion_cj <- c(
            sum(e_predicted_cj <= e_quantiles_cj[1]),
            sum(e_predicted_cj <= e_quantiles_cj[2]) - sum(e_predicted_cj <= e_quantiles_cj[1]),
            sum(e_predicted_cj <= e_quantiles_cj[3]) - sum(e_predicted_cj <= e_quantiles_cj[2]),
            sum(e_predicted_cj <= e_quantiles_cj[4]) - sum(e_predicted_cj <= e_quantiles_cj[3]),
            sum(e_predicted_cj <= e_quantiles_cj[5]) - sum(e_predicted_cj <= e_quantiles_cj[4]),
            sum(e_predicted_cj > e_quantiles_cj[5])
          ) / dim(predictions)[1]
        }
        
        
        pred_props_cj <- c(pred_props_cj,c_pred_proportion_cj,e_pred_proportion_cj)
      }    
      
    }
    
    # Combine the quantiles for rts and cj
    if (binning) {obs_props <- c(obs_props,obs_props_cj)
    }else{obs_props <- c(obs_props,obs_props)}
    
    pred_props <- c(pred_props,pred_props_cj)
    # calculate chi square
    
    # calculate chi square
    chiSquare = sum( (obs_props - pred_props) ^ 2)
    return(chiSquare)
  }
}
