#======================================================================================================================#
# Script created by Mark Christie and Mike Ford (NOAA): contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 in 2013
# This script:  Simulates data and examines the effect of p-values vs sample sizes for RRS studies
# Usage notes:  Run in chunks

# Modified on 11/7/14 by Kyle Shedd to accomodate different stray rates and sampling proportions of offspring and
# computer power by t.test, permutaion test, and negative binomial GLM

# Updated on 11/10/14 by Kyle Shedd to change the mu.n (i.e. the mean R/S) from 8 (got from Christie et al. 2014, Chinook)
# to 5.39, which appears to be the average R/S for wild PWS pink salmon from BY 1960-2011 (missing BY 2009). These are rough
# data taken from Appendix A5 and A6 of FM 11-07, the escapement goal review with some newer info from recent AMRss

# Updated on 11/12/14 by Kyle Shedd to change the mu.n to 2. We are sampling escapement, not the total run. Thus, while
# the average R/S for a wild PWS pink may be 5.39, the escapement is managed to be constant (in theory) so the R/S should
# be ~ 2 for replacement (assuming monogamous mating)

# Updated 11/14/14 by Kyle shedd to change mu.n to 1 (replacement for a 50% stray rate stream)

# Also removed the low low stray rate (0.05)
#======================================================================================================================#
# Source files, import packages, set working directory, initialize variables

ls()
rm(list=ls(all=TRUE))
setwd("V:/WORK/Pink/AHRG/Parentage simulations/Mark Christie")

# Which run is this?
run=5


# Define stray rates to model
stray_rates=c(0.15,0.5)

# Define sample proportion of offspring
sample_prop_off_rates=c((1:6)/6)



##########################

# Functions to help
nbnom_variance=function(mu,size){
  invisible(mu+((mu^2)/size))
}

size2=function(mu1,size1,mu2){
  var1=nbnom_variance(mu1,size1)
  invisible((mu2^2)/(var1-mu2))
}

#library(coin)
library(MASS)
#library(lattice)

##### start of simulations to determine POWER ####################################################################################################
ptm <- proc.time()

rrs.values <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
n.spawners <- seq(from = 100, to = 3000 , by = 100)

stray = rep(stray_rates,each=length(sample_prop_off_rates))[run]
sample_prop_off = rep(sample_prop_off_rates,length(stray_rates))[run]

# number of trials to get power
trials = 2000                                           # number of independent simulations to do

mu.n = 1  # THESE SIMULATIONS ARE TOTALLY DEPENDENT ON THE MU AND SIZE OF THE NATURAL DISTRIBUTION OF R/S
size.n = 0.95
var.n = nbnom_variance(mu.n,size.n)

OUT <- NULL
for (r in 1:length(rrs.values)){
  errs <- rrs.values[r]
  mu.h = mu.n*errs
  size.h =  size2(mu.n,size.n,mu.h)
  
  for (z in 1:length(n.spawners)) {
    n.h = n.spawners[z]*stray  # number of hatchery F1 spawners
    n.w = n.spawners[z]*(1-stray)  # number of wild F1 spanwers
    
    OUT2 <- NULL
    for(i in 1:trials){
            
      NWoffspring=rnbinom(n=n.w,size=size.n,mu=mu.n)   #take a sample of offspring based on this distribution
      NHoffspring=rnbinom(n=n.h,size=size.h,mu=mu.h)

# figure out how to sample a proportion of the offspring, what will family sizes look like
      n.NWoffspring=sum(NWoffspring) # number off WILD offspring produced
      n.NWsampled=round(n.NWoffspring*sample_prop_off)
      NWsampled=sort(sample(1:n.NWoffspring,n.NWsampled,replace=FALSE)) # individual WILD offspring that were sampled (produced * sample rate)

      cumNWoffspring=c(0,cumsum(NWoffspring))
      NWassigned=sapply(seq_along(NWoffspring),function(i){length(NWsampled[(NWsampled <= cumNWoffspring[i+1]) & (NWsampled > cumNWoffspring[i])])})

#cumNWoffspring=cbind(c(0,cumsum(NWoffspring)[1:n.w-1]),cumsum(NWoffspring))
#NWassigned=apply(cumNWoffspring,1,function(row){length(NWsampled[(NWsampled <= row[2]) & (NWsampled > row[1])])})

      n.NHoffspring=sum(NHoffspring) # number off HATCHERY offspring produced
      n.NHsampled=round(n.NHoffspring*sample_prop_off)
      NHsampled=sort(sample(1:n.NHoffspring,n.NHsampled,replace=FALSE)) # individual HATCHERY offspring that were sampled (produced * sample rate)

      cumNHoffspring=c(0,cumsum(NHoffspring))
      NHassigned=sapply(seq_along(NHoffspring),function(i){length(NHsampled[(NHsampled <= cumNHoffspring[i+1]) & (NHsampled > cumNHoffspring[i])])})

#cumNHoffspring=cbind(c(0,cumsum(NHoffspring)[1:n.h-1]),cumsum(NHoffspring))
#NHassigned=apply(cumNHoffspring,1,function(row){length(NHsampled[(NHsampled <= row[2]) & (NHsampled > row[1])])})

#################################################################
# Do a one-sided test
# output should have t.test, permutation test, and Nbinom

      mydata=data.frame(nOff=c(NWassigned,NHassigned),Origin=c(rep("W",n.w),rep("H",n.h)))

# permutation test (1-tail)
      reps=10000
      true_diff=mean(NWassigned)-mean(NHassigned)
      Wfams=replicate(reps,sum(sample(c(NWassigned,NHassigned),size=n.w,replace=FALSE)))
      Hfams=sum(NWassigned,NHassigned)-Wfams
      rsW=Wfams/n.w
      rsH=Hfams/n.h
      diffs=rsW-rsH
      perm_1tail_pvalue=round(sum(diffs>=true_diff)/reps,4)

      #perm_1tail_pvalue2=round(as.numeric(pvalue(oneway_test(nOff~Origin, data=mydata, distribution=approximate(B=10000), alternative="less"))),4)
# negative binomial GLM (1-tail)
      fit <- glm.nb(nOff~Origin, data=mydata, init.theta=1, link=log)
      nbGLM_1tail_pvalue=round(pnorm(summary(fit)$coefficients[2,3],lower.tail=FALSE),4)

# t.test (1-tail)
      test <- t.test(NWassigned, NHassigned, alternative="greater")
      ttest_1tail_pvalue=round(test$p.value,4)

# building output
      out <- cbind(errs, n.h, n.w, n.NHsampled, n.NWsampled, perm_1tail_pvalue, nbGLM_1tail_pvalue, ttest_1tail_pvalue)
      OUT2 <- rbind(OUT2, out)
    }

# compiling to parameter level
    perm_n.sig <- length(which(OUT2[, "perm_1tail_pvalue"] < 0.05))
    perm_pwr   <- perm_n.sig/trials

    nbGLM_n.sig <- length(which(OUT2[, "nbGLM_1tail_pvalue"] < 0.05))
    nbGLM_pwr   <- nbGLM_n.sig/trials

    ttest_n.sig <- length(which(OUT2[, "ttest_1tail_pvalue"] < 0.05))
    ttest_pwr   <- ttest_n.sig/trials


    out2 <- OUT2[1, ]
    out2[6] <- perm_pwr
    out2[7] <- nbGLM_pwr
    out2[8] <- ttest_pwr


    OUT <- rbind(OUT, out2) 
  }
  
  
}  
proc.time() - ptm

rownames(OUT)=c(1:dim(OUT)[1])
colnames(OUT)=c("RRS","N.h.parents","N.w.parents","N.h.offspring","N.w.offspring","Perm Power","nbGLM Power","ttest Power")

write.table(OUT, paste("mu",mu.n,"/simulation_results_stray",stray,"_prop_",round(sample_prop_off,3),"_trials_",trials,"_muRSn_",mu.n,".txt",sep=''), col.names = TRUE, sep="\t")
