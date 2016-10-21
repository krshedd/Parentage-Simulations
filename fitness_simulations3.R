#======================================================================================================================#
# Script created by Mark Christie and Mike Ford (NOAA): contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 in 2013
# This script:  Simulates data and examines the effect of p-values vs sample sizes for RRS studies
# Usage notes:  Run in chunks
#======================================================================================================================#
# Source files, import packages, set working directory, initialize variables
setwd("V:/WORK/Pink/AHRG/Parentage simulations/Mark Christie")

#library(coin)
#library(MASS)
#library(lattice)
#======================================================================================================================#

##### start of simulations ###########

rrs.values <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
n.spawners <- seq(from = 10, to = 400 , by = 10)
trials = 5000                                           # number of independent simulations to do


OUT <- NULL
for (r in 1:length(rrs.values)){
  errs <- rrs.values[r]

  for (z in 1:length(n.spawners)) {
    n.h = n.spawners[z]  # number of hatchery F1 spawners
    n.w = n.spawners[z]  # number of wild F1 spanwers
    
    OUT2 <- NULL
    for(i in 1:trials){
      NWassigned=sample(rnbinom(3000,.95,0.1),n.w,replace=FALSE)   #take a sample of offspring based on this distribution
      NHassigned=sample(rnbinom(3000,.95,0.1),n.h)*errs
      test <- t.test(NWassigned, NHassigned)
      out <- cbind(errs, n.h, n.w, test$p.value)
      OUT2 <- rbind(OUT2, out)
    }
    
    n.sig <- length(which(OUT2[, 4] < 0.05))
    pwr   <- n.sig/length(OUT2[, 4])
    out2 <- OUT2[1, ]
    out2[4] <- pwr
    OUT <- rbind(OUT, out2) 
  }
  
    
}  


write.table(OUT, "simulation_results.txt", col.names = FALSE, sep="\t", append = TRUE)

############################################################################################
# Function to calculate negative binomial variance parameter from "size" term
nbnom_variance=function(mu,size){
  invisible(mu+((mu^2)/size))
}

# What is variance from Christie nbinom distribution?
nbnom_variance(8,0.95) # 75.37
# What would variance be for same "size" term, but different mu?
nbnom_variance(4,0.95) # 20.84

# Function to determine the appropriate "size" parameter to be used with mu2 (new) to have the same variance as a separate distribution (mu1, size 1)
size2=function(mu1,size1,mu2){
  var1=nbnom_variance(mu1,size1)
  invisible((mu2^2)/(var1-mu2))
}

# Determine appropriate "size" parameter to give same variance as Christie distribution
size2(8,0.95,4)

# Verify that this "size" term gives the same variance as the Christie distribution, just with smaller mu.
nbnom_variance(4,0.224) # 

# Histogram of Christie distribution
hist(rnbinom(n=3000,size=0.95,mu=8),xlim=c(0,60),ylim=c(0,1600),breaks=200)
# Histogram of new distribution with smaller mu, but same Christie variance
hist(rnbinom(n=3000,size=size2(8,0.95,4),mu=4),xlim=c(0,60),ylim=c(0,1600),breaks=200)




############################################################################################

ls()
rm(list=ls(all=TRUE))
setwd("V:/WORK/Pink/AHRG/Parentage simulations/Mark Christie")


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

##### start of simulations ###########
ptm <- proc.time()

rrs.values <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
n.spawners <- seq(from = 100, to = 3000 , by = 100)
stray = 0.5
sample_prop_off = 0.5

# number of trials to get power
trials = 2                                           # number of independent simulations to do

mu.n = 8  # THESE SIMULATIONS ARE TOTALLY DEPENDENT ON THE MU AND SIZE OF THE NATURAL DISTRIBUTION OF R/S
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
      
      
      n.NHoffspring=sum(NHoffspring) # number off HATCHERY offspring produced
      n.NHsampled=round(n.NHoffspring*sample_prop_off)
      NHsampled=sort(sample(1:n.NHoffspring,n.NHsampled,replace=FALSE)) # individual HATCHERY offspring that were sampled (produced * sample rate)
      
      cumNHoffspring=c(0,cumsum(NHoffspring))
      NHassigned=sapply(seq_along(NHoffspring),function(i){length(NHsampled[(NHsampled <= cumNHoffspring[i+1]) & (NHsampled > cumNHoffspring[i])])})
      
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

write.table(OUT, paste("simulation_results_stray",stray,"_prop_",sample_prop_off,"_trials_",trials,".txt",sep=''), col.names = FALSE, sep="\t", append = TRUE)

str(OUT)
plot(OUT[,6]~apply(OUT,1,function(x){sum(x[2],x[3])}),pch=16,col=rep(c(1,2,3,4,5,6),each=30),cex=2)