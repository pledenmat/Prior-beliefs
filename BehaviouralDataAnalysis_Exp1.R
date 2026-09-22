## EXP 1: Fake feedback manipulation ##-----------------------------------------------------------------------------------

# Load libraries
library(readxl)

#################################----------------------------------------------------------------------------------------------------------------------------------
## DATA PRERAPATION AND CLEANING ----
#################################-----------------------------------------------------------------------------------------------------------------------------------

rm(list=ls())
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)
source("fastmerge.R")

# Import data
data_files <- list.files(path=c("Data/Data_1"),pattern = "*.csv", recursive = FALSE)

for(i in 1:length(data_files)){
  if(i == 1){ 
    AllData <- read.csv(paste0("Data/Data_1/",data_files[i]))
  }else{
    temp <- read.csv(paste0("Data/Data_1/",data_files[i]))
    AllData <- fastmerge(AllData,temp)
  }
}


## SUBJECT NUMBERS
# Subject number 41938 should be 38 
# change for plotting and looping purposes
sort(unique(AllData$sub))
for(i in 1:nrow(AllData)){
  AllData$sub[AllData$sub==41938] <- 38
}
sort(unique(AllData$sub))


## SPLIT DATA AND POSTCHECKS ---------------------------------------------------

Data <- subset(subset(AllData,task!=""), select=-c(post1,post2,post3,post4,post5,post6))
postchecks <- read.csv(paste0("RecodedPostchecks_Exp1.csv")) #Recoded postcheck answers

head(Data)
head(postchecks)

names(Data)[names(Data) == 'selfconf'] <- 'fbcond'

# Drop training data
Data <- subset(Data, running == "main")
table(Data$sub, Data$fbcond)

# Number of subjects
N = length(unique(Data$sub))

# New variable: 
x <- rep(c(rep("firsttask",216),rep("secondtask",216),rep("thirdtask",216)),N)
Data$whichtask <- x


## DESCRIPTIVES ----

## demographics
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

Data <- subset(Data,!(sub %in% exclusion)) #Sub 10 and 49 removed
postchecks <- subset(postchecks,!(sub %in% exclusion)) #Sub 10 and 49 removed

length(unique(Data$sub)) # 48 pps left
length(unique(postchecks$sub))




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
library(ggmcmc)
library(ggthemes)
library(ggridges)
library(matrixStats)

# FIXED
postchecks$SubjectiveInfluenceFb <- factor(postchecks$SubjectiveInfluenceFb)
postchecks$FeedbackCredibility   <- factor(postchecks$FeedbackCredibility)
postchecks$ManipulationAwareness <- factor(postchecks$ManipulationAwareness)

posts_table <- table(postchecks$FeedbackCredibility,postchecks$ManipulationAwareness)
posts_table <- posts_table[,c(2,3,1)]


maindata <- subset(Data, running == "main")
head(maindata)

# explicitly create factors
maindata$block <- factor(maindata$block)
maindata$fbcond <- factor(maindata$fbcond, levels=c("lowSC","mediumSC","highSC"))
maindata$difflevel <- factor(maindata$difflevel, levels=c("easy","average","hard"))
maindata$cor <- factor(maindata$cor)

maindata$fbcond <- factor(maindata$fbcond)

# Set coding schemes
options(contrasts = c("contr.sum","contr.poly"))
contrasts(maindata$fbcond)
contrasts(maindata$difflevel)



##--------------------
## MAIN ANALYSIS: CJ -
##--------------------
# Is there an effect of FEEDBACK on CONFIDENCE JUDGEMENTS (controlling for trial difficulty)?

nullmod <- lmer(cj ~ 1 + fbcond*difflevel + (1 | sub), data = maindata)

# Add random slopes 
fit1 <- lmer(cj ~ 1 + fbcond*difflevel + (1 + fbcond | sub), data = maindata,
             control = lmerControl(optimizer = "bobyqa"))
anova(nullmod,fit1) 
#
fit2 <- lmer(cj ~ 1 + fbcond*difflevel + (1 + difflevel | sub), data = maindata,
              control = lmerControl(optimizer = "bobyqa"))
anova(nullmod,fit2) 


## Both random slopes significant
#  AIC: 89299 fbcond slope
#       90368 difflevel slope
# Add both 
fit3 <- lmer(cj ~ 1 + fbcond*difflevel + (1 +difflevel+fbcond | sub), data = maindata,
              control = lmerControl(optimizer = "bobyqa"))
anova(nullmod,fit3)
# Converges, AIC 89046
# Best fit: (1 + difflevel + fbcond | sub)

# FIX (reproducibility report, Sept 2026): every reference below this point was
# `fit11`, a model that is never defined anywhere in this file.  
# Comment suggests that fit3 is the winning model (matches random-effects structure)
# so fit11 -> fit3 throughout this section.
summary(fit3)
Anova(fit3,type="III") # Doesn't produce the F(2,47)=16.65 test (original anova call)
anova(fit3) # This does but wasn't in the OSF script, might want to double-check all stats

## plot the effects
plot(allEffects(fit3))

# FIX NEEDED (reproducibility report, Sept 2026): this hand-built contrast matrix C
# is 3x27, but fit3 only has 9 fixed-effect coefficients , so glht(fit3, linfct = C) 
# errors with "'ncol(linfct)' is not equal to 'length(coef(model))'". 
# 27 columns is exactly what a 3-way interaction with another 3-level factor would need, 
# which matches fitFC/fitMA/fitSI (the postcheck follow-up models further below, e.g.
# `fbcond*difflevel*FeedbackCredibility`) rather than fit3. This looks like another
# casualty of the same renaming/refactor the other bugs came from: this post-hoc
# block was very likely written against one of those 3-way models, not fit3.
# testInteractions() below doesn't have this problem, so it's left running on fit3.
C <- rep(0,27)
C <- rep(C,3)
C <- matrix(C,3)
C[1,1] <- 1
C[2,2] <- 1
C[3,3] <- 1
C
# summary(glht(fit3, linfct = C)) 
summary(glht(fit3))

library(phia)
testInteractions(fit3, fixed = c("difflevel","fbcond"))



##----------------
## FOLLOW-UP: COR-
##----------------
# Does FEEDBACK also influence objective performance such as ACCURACY?

nullmod_cor <- glmer(cor ~ 1 + fbcond*difflevel + (1 | sub), 
                     data = maindata, 
                     control = glmerControl(optimizer = "bobyqa"),
                     family = binomial)

# Add random slopes
fitcor1 <- glmer(cor ~ 1 + fbcond*difflevel + (1 + fbcond | sub), 
                 data = maindata, 
                 control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)),
                 family = binomial)
anova(nullmod_cor,fitcor1)
#
fitcor2 <- glmer(cor ~ 1 + fbcond*difflevel + (1 + difflevel | sub), 
                 data = maindata, 
                 control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)),
                 family = binomial)
anova(nullmod_cor,fitcor2)

# Both slopes significant but difflevel is overfitted  
# --> Final model: only fbcond slope 

Anova(fitcor1,type="III")



##---------------
## FOLLOW-UP: RT-
##---------------
# Does FEEDBACK also influence objective performance such as REACTION TIME?

nullmod_rt <- lmer(rt ~ 1 + fbcond*difflevel + (1 | sub), data = maindata)
#
# adding random slopes
#
fitrt1 <- lmer(rt ~ 1 + fbcond*difflevel + (1 + fbcond | sub), data = maindata,
               control = lmerControl(optimizer = "bobyqa"))
anova(nullmod_rt,fitrt1)
#
fitrt2 <- lmer(rt ~ 1 + fbcond*difflevel + (1 + difflevel | sub), data = maindata,
               control = lmerControl(optimizer = "bobyqa"))
anova(nullmod_rt,fitrt2)

# Both slopes are significant but difflevel is overfitted  
# --> Final model: only fbcond slope 

Anova(fitrt1,type="III")



##------------------------------
## FOLLOW UP: EFFECT OVER TIME -
##------------------------------

#---
# 1. Effect over blocks
#---

# Model selection
block_fit0 <- lmer(cj ~ 1 + fbcond*difflevel*block + (1 | sub), data = maindata,
                    control = lmerControl(optimizer = "bobyqa"))

# Add random slopes
block_fit1 <- lmer(cj ~ 1 + fbcond*difflevel*block + (1 +fbcond | sub), data = maindata,
                    control = lmerControl(optimizer = "bobyqa"))
anova(block_fit0,block_fit1)
#
block_fit2 <- lmer(cj ~ 1 + fbcond*difflevel*block + (1 +difflevel | sub), data = maindata,
                    control = lmerControl(optimizer = "bobyqa"))
anova(block_fit0,block_fit2)


block_fit3 <- lmer(cj ~ 1 + fbcond*difflevel*block + (1 +fbcond*difflevel | sub), data = maindata,
                      control = lmerControl(optimizer = "bobyqa"))
# doesn't converge
# --> Final model: only fbcond slope (lowest AIC)
Anova(block_fit1, type="III")


#---
# 2. Effect within each block 
#---

data_block1 <- subset(Data,block==5)
data_block2 <- subset(Data,block==6)
data_block3 <- subset(Data,block==7)

## Model selection
b1_fit1 <- lmer(cj ~ 1 + fbcond*difflevel + (1 + fbcond | sub), data = data_block1,
                 control = lmerControl(optimizer = "bobyqa"))
b1_fit2 <- lmer(cj ~ 1 + fbcond*difflevel + (1 + difflevel | sub), data = data_block1,
                 control = lmerControl(optimizer = "bobyqa"))
anova(b1_fit1,b1_fit2) 
#
b2_fit1 <- lmer(cj ~ 1 + fbcond*difflevel + (1 + fbcond | sub), data = data_block2,
                 control = lmerControl(optimizer = "bobyqa"))
b2_fit2 <- lmer(cj ~ 1 + fbcond*difflevel + (1 + difflevel | sub), data = data_block2,
                 control = lmerControl(optimizer = "bobyqa"))
anova(b2_fit1,b2_fit2) # Best: fbcond slope
#
b3_fit1 <- lmer(cj ~ 1 + fbcond*difflevel + (1 + fbcond | sub), data = data_block3,
                 control = lmerControl(optimizer = "bobyqa"))
b3_fit2 <- lmer(cj ~ 1 + fbcond*difflevel + (1 + difflevel | sub), data = data_block3,
                 control = lmerControl(optimizer = "bobyqa"))
anova(b3_fit1,b3_fit2)

# outputs from winning models: 
anova(b1_fit1)
anova(b2_fit1)
anova(b3_fit1)



##-------------------------------------------------
## FOLLOW UP: Relationship with postcheck answers -
##-------------------------------------------------

# FIX NEEDED (reproducibility report, Sept 2026): depends on the coded postcheck
# columns flagged above (SubjectiveInfluenceFb/FeedbackCredibility/ManipulationAwareness),
# which are not available in the replication materials - see the comment near
# postchecks$SubjectiveInfluenceFb further up. Guarded so the rest of the script runs.
maindata_withpost <- merge(maindata,postchecks, by="sub")


fitFC <- lmer(cj ~ 1 + fbcond*difflevel*FeedbackCredibility +
                (1 +difflevel+fbcond | sub), data = maindata_withpost,
              control = lmerControl(optimizer = "bobyqa"))
anova(fitFC)
plot(allEffects(fitFC)) # FIX: was `fit11`, undefined here (see main fit11->fit3 fix above); this plot is of the model just fit


fitMA <- lmer(cj ~ 1 + fbcond*difflevel*ManipulationAwareness +
                (1 +difflevel+fbcond | sub), data = maindata_withpost,
              control = lmerControl(optimizer = "bobyqa"))
anova(fitMA)
plot(allEffects(fitMA))


fitSI <- lmer(cj ~ 1 + fbcond*difflevel*SubjectiveInfluenceFb +
                (1 +difflevel+fbcond | sub), data = maindata_withpost,
              control = lmerControl(optimizer = "bobyqa"))
anova(fitSI)

plot(allEffects(fitSI))





