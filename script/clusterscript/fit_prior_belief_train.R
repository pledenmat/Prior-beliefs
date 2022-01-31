rm(list=ls())
library(Rcpp) 
sourceCpp("DDM_with_confidence_slow.cpp")
source('fastmerge.r')
source('quantilefit_function_DDMonly_new.R')
library(DEoptim)
library(optparse)
option_list = list(
  make_option(c("-s","--subject"), type = "character",default=NULL,metavar="character")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
p = as.numeric(opt$subject)

# Data Load ---------------------------------------------------------------

Data <- read.table('data_exp1_training.csv',sep=',',header=T)
condLab <- unique(Data$condition); Ncond <- length(condLab)
N<-length(p);subs <- unique(Data$sub)
difficulty <- sort(unique(Data$coh));
Ntrials <- 10

# Model fit ---------------------------------------------------------------

for (c in 1:Ncond) {
  tempDat <- subset(Data,sub==subs[p]&condition==condlab[c])
  print(paste('Running participant',p,'from',N,'condition',c))
  #Load existing individual results if already exist
  file_name <- paste0('prior_belief/train/trainfit',condLab[c],subs[p],'.Rdata')
  if(file.exists(file_name)){
    load(file_name)
  }else{ #if not, fit the model
    optimal_params <- DEoptim(chi_square_optim, # function to optimize
                              lower = c(0, 0, 0, Ntrials, .1, .001,1,0,0,0), # a,ter,z,ntrials,sigma,dt,vratio,3 alpha,beta, 3 v
                              upper = c(.2, 2, 0, Ntrials, .1, .001,1,.5,.5,.5),
                              observations = tempDat,control=c(itermax=1000,steptol=100,reltol=.001,NP=50), returnFit = 1) # observed data is a parameter for the ks function we pass
    results <- summary(optimal_params)
    #save individual results
    save(results, file=file_name)
  }
  
  
}
