## EXP 2: Training difficulty manipulation ##-----------------------------------------------------------------------------------


# Load libraries
library(readxl)
library(multcomp)
library(reshape)
library(lme4)
library(lmerTest)
library(afex)
library(ggplot2)
library(effects)
library(car)
library(quickpsy)
library(BayesFactor)
library(ggmcmc)
library(ggthemes)
library(ggridges)

#################################----------------------------------------------------------------------------------------------------------------------------------
## DATA PRERAPATION AND CLEANING ----
#################################-----------------------------------------------------------------------------------------------------------------------------------

rm(list=ls())
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)
source("fastmerge.R")

# Import data
data_files <- list.files(path=c("Data/Data_2"),pattern = "*.csv", recursive = FALSE)

for(i in 1:length(data_files)){
  if(i == 1){ 
    AllData <- read.csv(paste0("Data/Data_2/",data_files[i]))
  }else{
    temp <- read.csv(paste0("Data/Data_2/",data_files[i]))
    AllData <- fastmerge(AllData,temp)
  }
}


## SPLIT DATA AND POSTCHECKS ---------------------------------------------------

Data <- subset(subset(AllData,task!=""), select=-c(post1,post2,post3,post4))
head(Data)

# drop unused levels
Data$task <- Data$task[drop=T]
Data$running <- Data$running[drop=T]
Data$traindiffcond <- Data$traindiffcond[drop=T]
Data$trialdifflevel <- Data$trialdifflevel[drop=T]
Data$resp <- Data$resp[drop=T]
Data$cresp <- Data$cresp[drop=T]
Data$conf_press <- Data$conf_press[drop=T]
Data$gender <- Data$gender[drop=T]
Data$handedness <- Data$handedness[drop=T]


train <- subset(Data, running == "training")
Data <- subset(Data, running == "main")

## ENCODE NUMBER OF PP
N = length(unique(Data$sub))

Data$traindiffcond <- factor(paste0(Data$traindiffcond,"training"))
unique(Data$traindiffcond)
contrasts(Data$traindiffcond)

x <- rep(c(rep("firsttask",216),rep("secondtask",216),rep("thirdtask",216)),N)
Data$whichtask <- x



## DESCRIPTIVES ----

# demographics
demographics <- subset(Data, task=="dotcolour"&block==5&withinblocktrial==1)
table(demographics$gender)
#
ages <- subset(demographics,age!=0)
mean(ages$age)
sd(ages$age)
range(ages$age)


## REMOVE OUTLIERS ----

#######
# RTS #
#######

mean(Data$rt)
table(Data$rt)
sd(Data$rt)

# get rid of obvious noise: everything below 100 ms and above 4000 ms
Data <- subset(Data, rt < 4000 & rt > 100)
mean(Data$rt)
sd(Data$rt)

#######
# ACC #
#######

Data['response'] <- 0
Data$response[Data$resp == "['n']"] <- 1

## Diagnostic plot per participant and task + chance performance testing
N <- length(unique(Data$sub)); subs <- unique(Data$sub); exclusion <- c()
tasks <- unique(Data$task)
par(mfrow=c(2,2))
for(i in 1:N){
  for (t in tasks) {
    tempDat <- subset(Data,sub==subs[i]&task==t)
    acc_block <- with(tempDat,aggregate(cor,by=list(block=block),mean))
    bias_block <- with(tempDat,aggregate(response,by=list(block=block),mean))
    test <- binom.test(length(tempDat$cor[tempDat$cor==1]),n=length(tempDat$cor),alternative = "greater")
    print(paste("In", t, "sub",subs[i],"p =", round(test$p.value,3),"compared to chance"))
    if (test$p.value > .05) {
      exclusion <- c(exclusion,subs[i])
      plot(acc_block,ylab="Acc (.) and bias (x)",frame=F,ylim=c(0,1));abline(h=.5,lty=2,col="grey")
      points(bias_block,pch=4)
      plot(tempDat$rt/1000,frame=F,col=c("black"),main=paste('subject',i,"task :",t),ylab="RT",ylim=c(0,5))
      plot(tempDat$cj,frame=F,col=c("black"),ylim=c(1,6),ylab="conf")
      plot(tempDat$RTconf,frame=F,col=c("black"),ylab="RT_conf")
    }
  }
}

Data <- subset(Data,!(sub %in% exclusion)) #Sub 12, 32 and 50 removed 

length(unique(Data$sub)) # 47 pps left



#################################----------------------------------------------------------------------------------------------------------------------------------
## DATA ANALYSIS ----
#################################-----------------------------------------------------------------------------------------------------------------------------------

## DATA PREPARATION ---------------------

# Load libraries
library(readxl)
library(multcomp)
library(reshape)
library(lme4)
library(lmerTest)
library(afex)
library(ggplot2)
library(effects)
library(car)
library(quickpsy)
library(BayesFactor)
library(ggmcmc)
library(ggthemes)
library(ggridges)
library(matrixStats)


maindata <- subset(Data, running == "main")

# explicitly create factors
maindata$block <- factor(maindata$block)
maindata$traindiffcond <- factor(maindata$traindiffcond, levels = c("easytraining","averagetraining","hardtraining"))
maindata$trialdifflevel <- factor(maindata$trialdifflevel)
maindata$cor <- factor(maindata$cor)


# Set coding schemes
options(contrasts = c("contr.sum","contr.poly"))
contrasts(maindata$traindiffcond)
contrasts(maindata$trialdifflevel)


##--------------------
## MAIN ANALYSIS: CJ -
##--------------------


nullmod <- lmer(cj ~ 1 + traindiffcond*trialdifflevel + (1 | sub), data = maindata)

### adding random slopes 
fit1 <- lmer(cj ~ 1 + traindiffcond*trialdifflevel + (1 + traindiffcond | sub), data = maindata,
             control = lmerControl(optimizer = "bobyqa"))
anova(nullmod,fit1) 
#
fit2 <- lmer(cj ~ 1 + traindiffcond*trialdifflevel + (1 + trialdifflevel | sub), data = maindata,
              control = lmerControl(optimizer = "bobyqa"))
anova(nullmod,fit2) 


summary(fit1)
Anova(fit1, type="III")

plot(allEffects(fit1))
plot(effect("traindiffcond", fit1))


##----------------
## FOLLOW-UP: COR-
##----------------

nullmod_cor <- glmer(cor ~ 1 + traindiffcond*trialdifflevel + (1 | sub), 
                     data = maindata, 
                     control = glmerControl(optimizer = "bobyqa"),
                     family = binomial)
#
# Adding random slopes
# 
fitcor1 <- glmer(cor ~ 1 + traindiffcond*trialdifflevel + (1 + traindiffcond | sub), 
                 data = maindata, 
                 control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)),
                 family = binomial)
anova(nullmod_cor,fitcor1)
#
fitcor2 <- glmer(cor ~ 1 + traindiffcond*trialdifflevel + (1 + trialdifflevel | sub), 
                 data = maindata, 
                 control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)),
                 family = binomial)
anova(nullmod_cor,fitcor2)


Anova(fitcor1,type="III")



##----------------
## FOLLOW-UP: RT-
##----------------

nullmod_rt <- lmer(rt ~ 1 + traindiffcond*trialdifflevel + (1 | sub), data = maindata)
#
# adding random slopes
#
fitrt1 <- lmer(rt ~ 1 + traindiffcond*trialdifflevel + (1 + traindiffcond | sub), data = maindata,
               control = lmerControl(optimizer = "bobyqa"))
anova(nullmod_rt,fitrt1)
#
fitrt2 <- lmer(rt ~ 1 + traindiffcond*trialdifflevel + (1 + trialdifflevel | sub), data = maindata,
               control = lmerControl(optimizer = "bobyqa"))
anova(nullmod_rt,fitrt2)
#
fitrt3 <- lmer(rt ~ 1 + traindiffcond*trialdifflevel + (1 + traindiffcond+trialdifflevel | sub), data = maindata,
               control = lmerControl(optimizer = "bobyqa"))
anova(fitrt1,fitrt3)
#
# adding both fails to converge
Anova(fitrt1, type="III")




##------------------------------
## FOLLOW UP: EFFECT OVER TIME -
##------------------------------

#---
# 1. Effect over blocks
#---

block_fit1 <- lmer(cj ~ 1 + traindiffcond*trialdifflevel*block + (1  | sub), data = maindata,
                   control = lmerControl(optimizer = "bobyqa"))
# Add random slopes
block_fit11 <- lmer(cj ~ 1 + traindiffcond*trialdifflevel*block + (1 +traindiffcond | sub), data = maindata,
                    control = lmerControl(optimizer = "bobyqa"))
anova(block_fit1,block_fit11)
#
block_fit12 <- lmer(cj ~ 1 + traindiffcond*trialdifflevel*block + (1 +trialdifflevel | sub), data = maindata,
                    control = lmerControl(optimizer = "bobyqa"))
anova(block_fit1,block_fit12)

block_fit13 <- lmer(cj ~ 1 + traindiffcond*trialdifflevel*block + (1 +trialdifflevel+traindiffcond | sub), data = maindata,
                    control = lmerControl(optimizer = "bobyqa"))
# fails to converge

Anova(block_fit11,type="III")



#---
# 2. Effect within each block 
#---

data_block1 <- subset(Data,block==5)
data_block2 <- subset(Data,block==6)
data_block3 <- subset(Data,block==7)

# Model selection
b1_fit12 <- lmer(cj ~ 1 + traindiffcond*trialdifflevel + (1 | sub), data = data_block1,
                 control = lmerControl(optimizer = "bobyqa"))
b1_fit13 <- lmer(cj ~ 1 + traindiffcond*trialdifflevel + (1 + traindiffcond| sub), data = data_block1,
                 control = lmerControl(optimizer = "bobyqa"))
b1_fit14 <- lmer(cj ~ 1 + traindiffcond*trialdifflevel + (1 + trialdifflevel| sub), data = data_block1,
                 control = lmerControl(optimizer = "bobyqa")) # overfit
anova(b1_fit12,b1_fit13)
#
b2_fit12 <- lmer(cj ~ 1 + traindiffcond*trialdifflevel + (1 | sub), data = data_block2,
                 control = lmerControl(optimizer = "bobyqa"))
b2_fit13 <- lmer(cj ~ 1 + traindiffcond*trialdifflevel + (1 + traindiffcond | sub), data = data_block2,
                 control = lmerControl(optimizer = "bobyqa"))
b2_fit14 <- lmer(cj ~ 1 + traindiffcond*trialdifflevel + (1 + trialdifflevel | sub), data = data_block2,
                 control = lmerControl(optimizer = "bobyqa")) # doesn't converge
anova(b2_fit12,b2_fit13)
#
b3_fit12 <- lmer(cj ~ 1 + traindiffcond*trialdifflevel + (1 | sub), data = data_block3,
                 control = lmerControl(optimizer = "bobyqa"))
b3_fit13 <- lmer(cj ~ 1 + traindiffcond*trialdifflevel + (1 + traindiffcond | sub), data = data_block3,
                 control = lmerControl(optimizer = "bobyqa"))
b3_fit14 <- lmer(cj ~ 1 + traindiffcond*trialdifflevel + (1 + trialdifflevel | sub), data = data_block3,
                 control = lmerControl(optimizer = "bobyqa")) # fits but higher AIC
b3_fit15 <- lmer(cj ~ 1 + traindiffcond*trialdifflevel + (1 + trialdifflevel + traindiffcond | sub), data = data_block3,
                 control = lmerControl(optimizer = "bobyqa")) 
anova(b3_fit12,b3_fit13)
anova(b3_fit12,b3_fit14) 
anova(b3_fit12,b3_fit15) # best model includes both slopes

# outputs from winning models: 
anova(b1_fit13)
anova(b2_fit13)
anova(b3_fit15)


