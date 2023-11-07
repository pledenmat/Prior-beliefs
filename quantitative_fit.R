curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)
source("1_preprocessing.R")
source("quantile_optim_DDM_Vs_bias.R")
library(lmerTest)
library(reshape)
library(timeSeries) #colSds
library(emmeans)
library(car)

# Plot parameters & functions ---------------------------------------------------------
error.bar <- function(x, y, upper, lower=upper, length=0,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
transp <- function(color, percent = 70, name = NULL) {
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

cexkl <- 1.5;cexgr <- 2;lwdgr <- 3; lwddat <- 2
cex_lab <- 3; cex_legend <- 2.5; cex_title <- 2.25 
cex_axis <- 1.5

## Data load ====
Data_exp1$fbcond <- Data_exp1$selfconf
Data_exp1$trialdifflevel <- Data_exp1$coh
subs1 <- sort(unique(Data_exp1$sub)); Nsub_1 <- length(subs1) 
cond_1 <- sort(unique(Data_exp1$fbcond)); Ncond <- length(cond_1)

Data_exp2$trialdifflevel <- Data_exp2$coh
subs2 <- sort(unique(Data_exp2$sub)); Nsub_2 <- length(subs2) 
cond_2 <- sort(unique(Data_exp2$traindiffcond)); Ncond <- length(cond_2)

trialdifflevel <- sort(unique(Data_exp1$trialdifflevel));Ndiff <- length(trialdifflevel)
# Fits load Exp 1 ---------------------------------------------------------------
setwd("fits")

#' Parameters null model
bound <- matrix(NA,Nsub_1,Ncond);ter <- matrix(NA,Nsub_1,Ncond)
#Adjust the number of drift parameters to the model loaded
v1 <- matrix(NA,Nsub_1,Ncond) 
v2 <- matrix(NA,Nsub_1,Ncond) 
v3 <- matrix(NA,Nsub_1,Ncond) 
resid <- matrix(NA,Nsub_1,Ncond)

#' 2 Step models
ddm_params <- data.frame(drift=NA,bound=NA,ter=NA,sub=rep(subs1,Ndiff*Ncond),
                         condition=rep(cond_1,each=Nsub_1,length.out=Nsub_1*Ncond*Ndiff),
                         difflevel=rep(trialdifflevel,each=Nsub_1*Ncond))
bound <- matrix(NA,Nsub_1,Ncond);ter <- matrix(NA,Nsub_1,Ncond)
v1 <- matrix(NA,Nsub_1,Ncond) 
v2 <- matrix(NA,Nsub_1,Ncond) 
v3 <- matrix(NA,Nsub_1,Ncond) 
resid <- matrix(NA,Nsub_1,Ncond)


#' Parameters bias fixed + v-ratio model
v_s <- matrix(NA,Nsub_1,Ncond) 
bias <- matrix(NA,Nsub_1,Ncond) 
resid <- matrix(NA,Nsub_1,Ncond)

v_s_bias <- matrix(NA,Nsub_1,Ncond) 
bias_bias <- matrix(NA,Nsub_1,Ncond) 
resid_bias <- matrix(NA,Nsub_1,Ncond)

v_s_bias_vssub <- matrix(NA,Nsub_1,Ncond) 
bias_bias_vssub <- matrix(NA,Nsub_1,Ncond) 
resid_bias_vssub <- matrix(NA,Nsub_1,Ncond)

v_s_full <- matrix(NA,Nsub_1,Ncond) 
bias_full <- matrix(NA,Nsub_1,Ncond) 
resid_full <- matrix(NA,Nsub_1,Ncond)

for(i in 1:Nsub_1){
  print(paste('Running participant',i,'from',Nsub_1))
  file_ddm <- paste0('Exp1/ddm/testfit',subs1[i],'.Rdata')
  file_conf <- paste0('Exp1/v_s_biassub/testfit',subs1[i],'.Rdata')
  file_conf_bias <- paste0('Exp1/bias/testfit',subs1[i],'.Rdata')
  file_conf_bias_vssub <- paste0('Exp1/v_ssub_bias/testfit',subs1[i],'.Rdata')
  file_conf_full <- paste0('Exp1/v_s_bias/testfit',subs1[i],'.Rdata')
  file_ok <- all(file.exists(file_conf,file_ddm,file_conf_bias,file_conf_bias_vssub,file_conf_full))
  if(file_ok|TRUE){
    load(file_ddm)
    bound[i,] <- results_ddm$optim$bestmem[1]
    ter[i,] <- results_ddm$optim$bestmem[2]
    v1[i,] <- results_ddm$optim$bestmem[7]
    v2[i,] <- results_ddm$optim$bestmem[8]
    v3[i,] <- results_ddm$optim$bestmem[9]
    resid[i,] <- results_ddm$optim$bestval
    ddm_params[ddm_params$sub==subs1[i],"drift"] <- results_ddm$optim$bestmem[7:9]
    ddm_params[ddm_params$sub==subs1[i],"bound"] <- results_ddm$optim$bestmem[1] 
    ddm_params[ddm_params$sub==subs1[i],"ter"] <- results_ddm$optim$bestmem[2]
    
    load(file_conf)
    v_s[i,] <- results$optim$bestmem[9:11]
    bias[i,] <- results$optim$bestmem[5]
    resid[i,] <- results$optim$bestval
    
    load(file_conf_bias)
    v_s_bias[i,] <- results$optim$bestmem[11]
    bias_bias[i,] <- results$optim$bestmem[5:7]
    resid_bias[i,] <- results$optim$bestval
    
    load(file_conf_bias_vssub)
    v_s_bias_vssub[i,] <- results$optim$bestmem[11]
    bias_bias_vssub[i,] <- results$optim$bestmem[5:7]
    resid_bias_vssub[i,] <- results$optim$bestval
    
    load(file_conf_full)
    v_s_full[i,] <- results$optim$bestmem[11:13]
    bias_full[i,] <- results$optim$bestmem[5:7]
    resid_full[i,] <- results$optim$bestval
  }
  else{ 
    print(paste("Unable to retrieve fit for",subs1[Nsub_1]))
  }
}
setwd("..")

params <- data.frame(drift = c(v1,v2,v3),bound=rep(bound,Ndiff),
                     ter=rep(ter,Ndiff), sub=rep(subs1,Ndiff*Ncond), 
                     v_s = rep(v_s,Ndiff), bias = rep(bias,Ndiff), 
                     condition=rep(cond_1,each=Nsub_1,length.out=Nsub_1*Ncond*Ndiff),
                     difflevel=rep(trialdifflevel,each=Nsub_1*Ncond),
                     exp=1,resid=rep(resid,Ndiff),model="v_s")
params$Npar <- 9

params_bias <- data.frame(drift = c(v1,v2,v3),bound=rep(bound,Ndiff),
                          ter=rep(ter,Ndiff), sub=rep(subs1,Ndiff*Ncond), 
                          v_s = rep(v_s_bias,Ndiff), bias = rep(bias_bias,Ndiff), 
                          condition=rep(cond_1,each=Nsub_1,length.out=Nsub_1*Ncond*Ndiff),
                          difflevel=rep(trialdifflevel,each=Nsub_1*Ncond),
                          exp=1,resid=rep(resid_bias,Ndiff),model="bias")
params_bias$Npar <- 8

params_bias_vssub <- data.frame(drift = c(v1,v2,v3),bound=rep(bound,Ndiff),
                                ter=rep(ter,Ndiff), sub=rep(subs1,Ndiff*Ncond), 
                                v_s = rep(v_s_bias_vssub,Ndiff), bias = rep(bias_bias_vssub,Ndiff), 
                                condition=rep(cond_1,each=Nsub_1,length.out=Nsub_1*Ncond*Ndiff),
                                difflevel=rep(trialdifflevel,each=Nsub_1*Ncond),
                                exp=1,resid=rep(resid_bias_vssub,Ndiff),model="bias_vssub")
params_bias_vssub$Npar <- 9

params_full <- data.frame(drift = c(v1,v2,v3),bound=rep(bound,Ndiff),
                          ter=rep(ter,Ndiff), sub=rep(subs1,Ndiff*Ncond),
                          v_s = rep(v_s_full,Ndiff), bias = rep(bias_full,Ndiff),
                          condition=rep(cond_1,each=Nsub_1,length.out=Nsub_1*Ncond*Ndiff),
                          difflevel=rep(trialdifflevel,each=Nsub_1*Ncond),
                          exp=1,resid=rep(resid_full,Ndiff),model="full")
params_full$Npar <- 11

params <- rbind(params,params_bias,params_bias_vssub,params_full)
# Generate model simulations ----------------------------------------------
rm(Simuls1)
for (s in 1:Nsub_1) {
  print(paste("Generating participant",s,"of",Nsub_1))
  tempddm <- subset(ddm_params,sub==subs1[s])
  tempDat <- subset(Data_exp1,sub==subs1[s])
  temp <- quantile_optim_DDM_Vs_biasfixed(c(0,1,.1,.001,bias[s],rep(1,Ncond),v_s[s,]),
                                          tempDat,returnFit=0,ddm_params=tempddm,
                                          condition_name = "selfconf")
  if(!exists('Simuls1')){ Simuls1 <- cbind(temp,subs1[s])
  }else{ Simuls1 <- rbind(Simuls1,cbind(temp,subs1[s]))
  }
}

names(Simuls1) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','condition','trialdifflevel','closest_evdnc2',"temp_rt2",'cj_raw','sub')
# Fits load Exp 2 ---------------------------------------------------------------
setwd("fits")

#' Parameters null model
bound <- matrix(NA,Nsub_2,Ncond);ter <- matrix(NA,Nsub_2,Ncond)
#Adjust the number of drift parameters to the model loaded
v1 <- matrix(NA,Nsub_2,Ncond) 
v2 <- matrix(NA,Nsub_2,Ncond) 
v3 <- matrix(NA,Nsub_2,Ncond) 
resid <- matrix(NA,Nsub_2,Ncond)

#' 2 Step models
ddm_params <- data.frame(drift=NA,bound=NA,ter=NA,sub=rep(subs2,Ndiff*Ncond),
                         condition=rep(cond_2,each=Nsub_2,length.out=Nsub_2*Ncond*Ndiff),
                         difflevel=rep(trialdifflevel,each=Nsub_2*Ncond))
bound <- matrix(NA,Nsub_2,Ncond);ter <- matrix(NA,Nsub_2,Ncond)
v1 <- matrix(NA,Nsub_2,Ncond) 
v2 <- matrix(NA,Nsub_2,Ncond) 
v3 <- matrix(NA,Nsub_2,Ncond) 
resid <- matrix(NA,Nsub_2,Ncond)


#' Parameters bias fixed + v-ratio model
v_s <- matrix(NA,Nsub_2,Ncond) 
bias <- matrix(NA,Nsub_2,Ncond) 
resid <- matrix(NA,Nsub_2,Ncond)

v_s_bias <- matrix(NA,Nsub_2,Ncond) 
bias_bias <- matrix(NA,Nsub_2,Ncond) 
resid_bias <- matrix(NA,Nsub_2,Ncond)

v_s_bias_vssub <- matrix(NA,Nsub_2,Ncond) 
bias_bias_vssub <- matrix(NA,Nsub_2,Ncond) 
resid_bias_vssub <- matrix(NA,Nsub_2,Ncond)

v_s_full <- matrix(NA,Nsub_2,Ncond) 
bias_full <- matrix(NA,Nsub_2,Ncond) 
resid_full <- matrix(NA,Nsub_2,Ncond)

for(i in 1:Nsub_2){
  print(paste('Running participant',i,'from',Nsub_2))
  file_ddm <- paste0('Exp2/ddm/testfit',subs2[i],'.Rdata')
  file_conf <- paste0('Exp2/v_s_biassub/testfit',subs2[i],'.Rdata')
  file_conf_bias <- paste0('Exp2/bias/testfit',subs2[i],'.Rdata')
  file_conf_bias_vssub <- paste0('Exp2/v_ssub_bias/testfit',subs2[i],'.Rdata')
  file_conf_full <- paste0('Exp2/v_s_bias/testfit',subs2[i],'.Rdata')
  file_ok <- all(file.exists(file_conf,file_ddm,file_conf_bias,file_conf_bias_vssub,file_conf_full))
  if(file_ok|T){
    load(file_ddm)
    bound[i,] <- results_ddm$optim$bestmem[1]
    ter[i,] <- results_ddm$optim$bestmem[2]
    v1[i,] <- results_ddm$optim$bestmem[7]
    v2[i,] <- results_ddm$optim$bestmem[8]
    v3[i,] <- results_ddm$optim$bestmem[9]
    resid[i,] <- results_ddm$optim$bestval
    ddm_params[ddm_params$sub==subs2[i],"drift"] <- results_ddm$optim$bestmem[7:9]
    ddm_params[ddm_params$sub==subs2[i],"bound"] <- results_ddm$optim$bestmem[1] 
    ddm_params[ddm_params$sub==subs2[i],"ter"] <- results_ddm$optim$bestmem[2]
    
    load(file_conf)
    v_s[i,] <- results$optim$bestmem[9:11]
    bias[i,] <- results$optim$bestmem[5]
    resid[i,] <- results$optim$bestval
    
    load(file_conf_bias)
    v_s_bias[i,] <- results$optim$bestmem[11]
    bias_bias[i,] <- results$optim$bestmem[5:7]
    resid_bias[i,] <- results$optim$bestval
    
    load(file_conf_bias_vssub)
    v_s_bias_vssub[i,] <- results$optim$bestmem[11]
    bias_bias_vssub[i,] <- results$optim$bestmem[5:7]
    resid_bias_vssub[i,] <- results$optim$bestval
    
    load(file_conf_full)
    v_s_full[i,] <- results$optim$bestmem[11:13]
    bias_full[i,] <- results$optim$bestmem[5:7]
    resid_full[i,] <- results$optim$bestval
  }
  else{ 
    print(paste("Unable to retrieve fit for",subs1[Nsub_2]))
  }
}
setwd("..")

params2 <- data.frame(drift = c(v1,v2,v3),bound=rep(bound,Ndiff),
                      ter=rep(ter,Ndiff), sub=rep(subs2,Ndiff*Ncond), 
                      v_s = rep(v_s,Ndiff), bias = rep(bias,Ndiff), 
                      condition=rep(cond_2,each=Nsub_2,length.out=Nsub_2*Ncond*Ndiff),
                      difflevel=rep(trialdifflevel,each=Nsub_2*Ncond),
                      exp=2,resid=rep(resid,Ndiff),model="v_s")
params2$Npar <- 9
params2_bias <- data.frame(drift = c(v1,v2,v3),bound=rep(bound,Ndiff),
                           ter=rep(ter,Ndiff), sub=rep(subs2,Ndiff*Ncond), 
                           v_s = rep(v_s_bias,Ndiff), bias = rep(bias_bias,Ndiff), 
                           condition=rep(cond_2,each=Nsub_2,length.out=Nsub_2*Ncond*Ndiff),
                           difflevel=rep(trialdifflevel,each=Nsub_2*Ncond),
                           exp=2,resid=rep(resid_bias,Ndiff),model="bias")
params2_bias$Npar <- 8

params2_bias_vssub <- data.frame(drift = c(v1,v2,v3),bound=rep(bound,Ndiff),
                                 ter=rep(ter,Ndiff), sub=rep(subs2,Ndiff*Ncond), 
                                 v_s = rep(v_s_bias_vssub,Ndiff), bias = rep(bias_bias_vssub,Ndiff), 
                                 condition=rep(cond_2,each=Nsub_2,length.out=Nsub_2*Ncond*Ndiff),
                                 difflevel=rep(trialdifflevel,each=Nsub_2*Ncond),
                                 exp=2,resid=rep(resid_bias_vssub,Ndiff),model="bias_vssub")
params2_bias_vssub$Npar <- 8

params2_full <- data.frame(drift = c(v1,v2,v3),bound=rep(bound,Ndiff),
                           ter=rep(ter,Ndiff), sub=rep(subs2,Ndiff*Ncond),
                           v_s = rep(v_s_full,Ndiff), bias = rep(bias_full,Ndiff),
                           condition=rep(cond_2,each=Nsub_2,length.out=Nsub_2*Ncond*Ndiff),
                           difflevel=rep(trialdifflevel,each=Nsub_2*Ncond),
                           exp=2,resid=rep(resid_full,Ndiff),model="full")
params2_full$Npar <- 11

params2 <- rbind(params2,params2_bias,params2_bias_vssub,params2_full)
# Generate model simulations ----------------------------------------------
rm(Simuls2)
for (s in 1:Nsub_2) {
  print(paste("Generating participant",s,"of",Nsub_2))
  tempddm <- subset(ddm_params,sub==subs2[s])
  tempDat <- subset(Data_exp2,sub==subs2[s])
  temp <- quantile_optim_DDM_Vs_biasfixed(c(0,1,.1,.001,mean(bias[s,]),rep(1,Ncond),v_s[s,]),
                                          tempDat,returnFit=0,ddm_params=tempddm,
                                          condition_name = "traindiffcond")
  if(!exists('Simuls2')){ Simuls2 <- cbind(temp,subs2[s])
  }else{ Simuls2 <- rbind(Simuls2,cbind(temp,subs2[s]))
  }
}

names(Simuls2) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','condition','trialdifflevel','closest_evdnc2',"temp_rt2",'cj_raw','sub')

# Model comparison --------------------------------------------------------
bic_custom <- function(Residuals,k,n){
  return(log(n)*k+n*log(Residuals/n))
}

param <- rbind(params,params2)
param$Ndata_points <- 12*Ncond # 6 CJ quantiles for correct/error
param$bic <- bic_custom(param$resid,param$Npar,param$Ndata_points)

mean_bic <- with(param,aggregate(bic,by=list(model=model,exp=exp),mean))
mean_resid <- with(param,aggregate(resid,by=list(model=model,exp=exp),mean))

mean_bic$delta <- -99
mean_bic[mean_bic$exp==1,"delta"] <- 
  mean_bic[mean_bic$exp==1,]$x - 
  min(mean_bic[mean_bic$exp==1,]$x)
mean_bic[mean_bic$exp==2,"delta"] <- 
  mean_bic[mean_bic$exp==2,]$x -
  min(mean_bic[mean_bic$exp==2,]$x)

# Stat tests --------------------------------------------------------------
Simuls1$condition <- as.factor(Simuls1$condition)
Simuls1$trialdifflevel <- as.factor(Simuls1$trialdifflevel)
#' Accuracy
m.int <- glmer(data = Simuls1, cor ~ condition*trialdifflevel + (1|sub),family = "binomial")
m.cond <- glmer(data = Simuls1, cor ~ condition*trialdifflevel + (condition|sub),
                family = "binomial")
m.diff <- glmer(data = Simuls1, cor ~ condition*trialdifflevel + (trialdifflevel|sub),
                family = "binomial")
anova(m.int,m.cond)
plot(resid(m.cond),Simuls1$cor) #Linearity
leveneTest(residuals(m.cond) ~ Simuls1$condition*Simuls1$trialdifflevel) #Homogeneity of variance
qqmath(m.cond) #Normality
Anova(m.cond)

#' RT
m.int <- lmer(data = Simuls1, rt ~ condition*trialdifflevel + (1|sub),REML = F)
m.cond <- lmer(data = Simuls1, rt ~ condition*trialdifflevel + (condition|sub),
               REML=F)
m.diff <- lmer(data = Simuls1, rt ~ condition*trialdifflevel + (trialdifflevel|sub),
               REML=F,control=lmerControl(optimizer = 'bobyqa'))
anova(m.int,m.cond)
plot(resid(m.cond),Simuls1$rt) #Linearity
leveneTest(residuals(m.cond) ~ Simuls1$condition*Simuls1$trialdifflevel) #Homogeneity of variance
qqmath(m.cond) #Normality
anova(m.cond)

# Also do the analysis with log-transformed RT for robustness of results
m.int <- lmer(data = Simuls1, log(rt) ~ condition*trialdifflevel + (1|sub),REML = F)
m.cond <- lmer(data = Simuls1, log(rt) ~ condition*trialdifflevel + (condition|sub),
               REML=F)
m.diff <- lmer(data = Simuls1, log(rt) ~ condition*trialdifflevel + (trialdifflevel|sub),
               REML=F,control=lmerControl(optimizer = 'bobyqa'))
anova(m.int,m.cond)
plot(resid(m.cond),log(Simuls1$rt)) #Linearity
leveneTest(residuals(m.cond) ~ Simuls1$condition*Simuls1$trialdifflevel) #Homogeneity of variance
qqmath(m.cond) #Normality looks better
anova(m.cond)

#' Confidence
m.int <- lmer(data = Simuls1, cj ~ condition*trialdifflevel + (1|sub),REML = F)
m.cond <- lmer(data = Simuls1, cj ~ condition*trialdifflevel + (condition|sub),
               REML=F,control=lmerControl(optimizer = 'bobyqa'))
m.diff <- lmer(data = Simuls1, cj ~ condition*trialdifflevel + (trialdifflevel|sub),
               REML=F,control=lmerControl(optimizer = 'bobyqa'))
anova(m.int,m.cond)
anova(m.int,m.diff)
anova(m.diff,m.cond)
m.both <- lmer(data = Simuls1, cj ~ condition*trialdifflevel + (condition+trialdifflevel|sub),
               REML=F,control=lmerControl(optimizer = 'bobyqa'))
anova(m.cond,m.both)
m.full <- lmer(data = Simuls1, cj ~ condition*trialdifflevel + (condition*trialdifflevel|sub),
               REML=F,control=lmerControl(optimizer = 'bobyqa'))
m.full2 <- lmer(data = Simuls1, cj ~ condition*trialdifflevel + (condition*trialdifflevel|sub),
                REML=F,control=lmerControl(optimizer = 'bobyqa',optCtrl = list(maxfun=15000)))
plot(resid(m.both),Simuls1$cj) #Linearity
leveneTest(residuals(m.both) ~ Simuls1$condition*Simuls1$trialdifflevel) #Homogeneity of variance
qqmath(m.both) #Normality
anova(m.both)

# Output ------------------------------------------------------------------
setwd("Data/aggregated")
write.csv(Simuls1,file = "quantitative_prediction_exp1.csv",row.names = F)
write.csv(Simuls2,file = "quantitative_prediction_exp2.csv",row.names = F)
setwd("..")
setwd("..")
# BF ----------------------------------------------------------------------

agg_cjData_exp1 <- with(Simuls1,aggregate(cj,by=list(sub=sub,trialdifflevel = trialdifflevel, condition=condition),mean));
agg_cjData_exp1$sub <- as.factor(agg_cjData_exp1$sub)
agg_cjData_exp1$condition <- as.factor(agg_cjData_exp1$condition)
agg_cjData_exp1$trialdifflevel <- as.factor(agg_cjData_exp1$trialdifflevel)
head(agg_cjData_exp1)
bf = BayesFactor::anovaBF(x ~ trialdifflevel*condition + sub, data = agg_cjData_exp1, whichRandom="sub")
bf[3]/bf[2] #main effect 1
bf[3]/bf[1] # main effect 2
bf[4]/bf[3] #interaction

agg_cjData_exp1 <- with(Simuls2,aggregate(cj,by=list(sub=sub,trialdifflevel = trialdifflevel, condition=condition),mean));
agg_cjData_exp1$sub <- as.factor(agg_cjData_exp1$sub)
agg_cjData_exp1$condition <- as.factor(agg_cjData_exp1$condition)
agg_cjData_exp1$trialdifflevel <- as.factor(agg_cjData_exp1$trialdifflevel)
head(agg_cjData_exp1)
bf = BayesFactor::anovaBF(x ~ trialdifflevel*condition + sub, data = agg_cjData_exp1, whichRandom="sub")
bf[3]/bf[2] #main effect 1
bf[3]/bf[1] # main effect 2
bf[3]/bf[4] #interaction (BF01)


# RT and accuracy scatter plots --------------------------------------------------------
jpeg(filename = "individual_fit_scatterplot.jpg",width=12,height=19.5,units='cm',res=600)
layout(matrix(c(1,3,5,7,2,4,6,8),ncol=2),heights = c(.1,.4,.4,.4))
par(mar=c(0,0,0,0))
plot.new()
title("Experiment 1",line=-2)
plot.new()
title("Experiment 2",line=-2)
par(mar=c(5,4,0,2)+.1)
## RT
# Experiment 1
sim_rt <- with(Simuls1,aggregate(rt,by=list(sub=sub),median))
dat_rt <- with(Data_exp1,aggregate(rt,by=list(sub=sub),median))
cor.test(sim_rt$x,dat_rt$x)

plot(sim_rt$x~dat_rt$x, bty='n',ylab="median RT - model fit", xlab = "median RT - behavior",
     xlim=c(0,1.4),ylim=c(0,1.4))
abline(a=0,b=1)

# Experiment 2
sim_rt <- with(Simuls2,aggregate(rt,by=list(sub=sub),median))
dat_rt <- with(Data_exp2,aggregate(rt,by=list(sub=sub),median))
cor.test(sim_rt$x,dat_rt$x)

plot(sim_rt$x~dat_rt$x, bty='n',ylab="median RT - model fit", xlab = "median RT - behavior",
     xlim=c(0,1.4),ylim=c(0,1.4))
abline(a=0,b=1)

# Accuracy
sim_cor <- with(Simuls1,aggregate(cor,by=list(sub=sub),mean))
dat_cor <- with(Data_exp1,aggregate(cor,by=list(sub=sub),mean))
cor.test(sim_cor$x,dat_cor$x)

plot(sim_cor$x~dat_cor$x, bty='n',ylab="Accuracy - model fit", xlab = "Accuracy - behavior",
     xlim=c(0.5,1),ylim=c(.5,1))
abline(a=0,b=1)

sim_cor <- with(Simuls2,aggregate(cor,by=list(sub=sub),mean))
dat_cor <- with(Data_exp2,aggregate(cor,by=list(sub=sub),mean))
cor.test(sim_cor$x,dat_cor$x)

plot(sim_cor$x~dat_cor$x, bty='n',ylab="Accuracy - model fit", xlab = "Accuracy - behavior",
     xlim=c(0.5,1),ylim=c(.5,1))
abline(a=0,b=1)

sim_cj <- with(Simuls2,aggregate(cj,by=list(sub=sub),mean))
dat_cj <- with(Data_exp2,aggregate(cj,by=list(sub=sub),mean))
cor.test(sim_cj$x,dat_cj$x)

## Confidence
plot(sim_cj$x~dat_cj$x, bty='n',ylab="Confidence - model fit", xlab = "Confidence - behavior",
     xlim=c(3,6),ylim=c(3,6))
abline(a=0,b=1)

sim_cj <- with(Simuls1,aggregate(cj,by=list(sub=sub),mean))
dat_cj <- with(Data_exp1,aggregate(cj,by=list(sub=sub),mean))
cor.test(sim_cj$x,dat_cj$x)

plot(sim_cj$x~dat_cj$x, bty='n',ylab="Confidence - model fit", xlab = "Confidence - behavior",
     xlim=c(3,6),ylim=c(3,6))
abline(a=0,b=1)

dev.off()
par(mar=c(5,4,4,2)+.1)


# Plot CJ & RT distribution - Experiment 1 ----------------------------------------------------
Simuls1[Simuls1$rt>5,"rt"] <- 5
Simuls2[Simuls2$rt>5,"rt"] <- 5
breaks <- seq(0,5,.25)
setwd("plots")
setwd("indiv")
setwd("exp1")
for (s in 1:Nsub_1) {
  jpeg(filename = paste0("dist_rt_conf_sub_exp1_",subs1[s],".jpg"),width=16,height=10,units="cm",res = 600)
  layout(matrix(c(1,2,1,3),ncol=2),heights = c(.1,.9))
  par(mar=c(0,0,0,0))
  plot.new()
  title(paste("Participant",subs1[s]),line=-1)
  par(mar=c(5,4,0,2))
  
  # Confidence
  c_data <- Data_exp1[Data_exp1$sub == subs1[s] &Data_exp1$cor == 1,]
  e_data <- Data_exp1[Data_exp1$sub == subs1[s] &Data_exp1$cor == 0,]
  c_simul <- Simuls1[Simuls1$sub == subs1[s] &Simuls1$cor==1,]
  e_simul <- Simuls1[Simuls1$sub == subs1[s] &Simuls1$cor==0,]
  ylim <- c(0, max(max(table(c_simul$cj)), max(table(c_data$cj))))
  tempC <- hist(c_data$cj, breaks = 0:6, xlim = c(0,6), prob = F,ylim=ylim, 
                col = rgb(0,1,0,.25), border = "white", 
                ylab = "Frequency", xlab = "Confidence",
                cex.lab = 1, cex.main = 1, cex.axis = 1, main = "")
  tempE <- hist(e_data$cj,breaks=0:6,prob=F,add=T,col=rgb(1,0,0,.25),border='white')
  Cors <- hist(c_simul$cj,breaks=0:6,plot=F)
  Errs <- hist(e_simul$cj,breaks=0:6,plot=F)
  lines(Cors$mids,Cors$counts,type='l',col='green',lwd=3)
  lines(Errs$mids,Errs$counts,type='l',col='red',lwd=3)
  legend("topleft",fill=c("white","white","green","red"),border=F,cex=.8,
         legend=c("Simulated corrects","Simulated errors","Empirical corrects","Empirical errors"),
         col=rep(c("Green","Red"),2),bty='n',lwd=c(1,1,-1,-1))
  
  
  # RT
  c_data <- Data_exp1[Data_exp1$sub == subs1[s] &Data_exp1$cor == 1,]
  e_data <- Data_exp1[Data_exp1$sub == subs1[s] &Data_exp1$cor == 0,]
  c_simul <- Simuls1[Simuls1$sub == subs1[s] &Simuls1$cor==1,]
  e_simul <- Simuls1[Simuls1$sub == subs1[s] &Simuls1$cor==0,]
  ylim <- c(0, max(max(table(cut(c_simul$rt,breaks))), max(table(cut(c_data$rt,breaks=breaks))),
                   max(table(cut(e_simul$rt,breaks))), max(table(cut(e_data$rt,breaks)))))
  tempC <- hist(c_data$rt, breaks = breaks, xlim = c(0,5), prob = F,ylim=ylim,
                col = rgb(0,1,0,.25), border = "white",
                ylab = "Frequency", xlab = "Reaction time (s)",
                cex.lab = 1, cex.main = 1, cex.axis = 1, main = "")
  tempE <- hist(e_data$rt,breaks=breaks,prob=F,add=T,col=rgb(1,0,0,.25),border='white')
  Cors <- hist(c_simul$rt,breaks=breaks,plot=F)
  Errs <- hist(e_simul$rt,breaks=breaks,plot=F)
  lines(Cors$mids,Cors$counts,type='l',col='green',lwd=3)
  lines(Errs$mids,Errs$counts,type='l',col='red',lwd=3)
  
  dev.off()
}
setwd("..")
setwd("exp2")
for (s in 1:Nsub_1) {
  jpeg(filename = paste0("dist_rt_conf_sub_exp2_",subs2[s],".jpg"),width=16,height=10,units="cm",res = 600)
  layout(matrix(c(1,2,1,3),ncol=2),heights = c(.1,.9))
  par(mar=c(0,0,0,0))
  plot.new()
  title(paste("Participant",subs2[s]),line=-1)
  par(mar=c(5,4,0,2))
  
  # Confidence
  c_data <- Data_exp2[Data_exp2$sub == subs2[s] &Data_exp2$cor == 1,]
  e_data <- Data_exp2[Data_exp2$sub == subs2[s] &Data_exp2$cor == 0,]
  c_simul <- Simuls2[Simuls2$sub == subs2[s] &Simuls2$cor==1,]
  e_simul <- Simuls2[Simuls2$sub == subs2[s] &Simuls2$cor==0,]
  ylim <- c(0, max(max(table(c_simul$cj)), max(table(c_data$cj))))
  tempC <- hist(c_data$cj, breaks = 0:6, xlim = c(0,6), prob = F,ylim=ylim, 
                col = rgb(0,1,0,.25), border = "white", 
                ylab = "Frequency", xlab = "Confidence",
                cex.lab = 1, cex.main = 1, cex.axis = 1, main = "")
  tempE <- hist(e_data$cj,breaks=0:6,prob=F,add=T,col=rgb(1,0,0,.25),border='white')
  Cors <- hist(c_simul$cj,breaks=0:6,plot=F)
  Errs <- hist(e_simul$cj,breaks=0:6,plot=F)
  lines(Cors$mids,Cors$counts,type='l',col='green',lwd=3)
  lines(Errs$mids,Errs$counts,type='l',col='red',lwd=3)
  legend("topleft",fill=c("white","white","green","red"),border=F,cex=.8,
         legend=c("Simulated corrects","Simulated errors","Empirical corrects","Empirical errors"),
         col=rep(c("Green","Red"),2),bty='n',lwd=c(1,1,-1,-1))
  
  
  # RT
  c_data <- Data_exp2[Data_exp2$sub == subs2[s] &Data_exp2$cor == 1,]
  e_data <- Data_exp2[Data_exp2$sub == subs2[s] &Data_exp2$cor == 0,]
  c_simul <- Simuls2[Simuls2$sub == subs2[s] &Simuls2$cor==1,]
  e_simul <- Simuls2[Simuls2$sub == subs2[s] &Simuls2$cor==0,]
  ylim <- c(0, max(max(table(cut(c_simul$rt,breaks))), max(table(cut(c_data$rt,breaks=breaks))),
                   max(table(cut(e_simul$rt,breaks))), max(table(cut(e_data$rt,breaks)))))
  tempC <- hist(c_data$rt, breaks = breaks, xlim = c(0,5), prob = F,ylim=ylim,
                col = rgb(0,1,0,.25), border = "white",
                ylab = "Frequency", xlab = "Reaction time (s)",
                cex.lab = 1, cex.main = 1, cex.axis = 1, main = "")
  tempE <- hist(e_data$rt,breaks=breaks,prob=F,add=T,col=rgb(1,0,0,.25),border='white')
  Cors <- hist(c_simul$rt,breaks=breaks,plot=F)
  Errs <- hist(e_simul$rt,breaks=breaks,plot=F)
  lines(Cors$mids,Cors$counts,type='l',col='green',lwd=3)
  lines(Errs$mids,Errs$counts,type='l',col='red',lwd=3)
  
  dev.off()
}
setwd('..')  
setwd('..')  
