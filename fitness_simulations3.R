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

# Updated 11/14/14 by Kyle Shedd to change mu.n to 1 (replacement for a 50% stray rate stream)

# Updated 11/14/14 by Kyle Shedd to change mu.n to 5.39 (average R/S for PWS pinks), then cull so that 1500 matings result
# in progeny + new strays = 3000

# Also removed the low low stray rate (0.05)

# Updated 11/18/14 by Kyle Shedd to change rrs values and proportion F1 values to better fill out power plot

# Updates 4/6/15 by Kyle Shedd to optimize for speed by vectorizing! (i.e. kill rbind and sapply and for loops)

# Updated 4/15/15 by Kyle Shedd to hold the number of parents sampled constant, but vary the distribution of RS

# Updated 4/21/15 by Kyle Shedd to change var.n = var.h to size.n = size.h after talking to Jim Jasper. Holding variance 
# constant seems unrealistic in the reduction of dispersion (size).

# Updated 4/21/15 by Kyle Shedd to skip trial iterations in for loop when mydata$nOff == 0 (i.e. when mu.n and size.n are
# very small and the number of parents caught are very small such as Paddy 2013).
#======================================================================================================================#
# Source files, import packages, set working directory, initialize variables

ls()
rm(list = ls(all = TRUE))
setwd("V:/WORK/Pink/AHRG/Parentage simulations/Mark Christie")

##########################

library(coin)
library(MASS)
library(beepr)
#library(lattice)

##### start of simulations to determine POWER ####################################################################################################
ptm <- proc.time()


## Stream specific otolith results for females
# Create a matrix of all stream/year/n.w/n.h values and then call them by row. Capiche?
parents <- read.table(file = "StreamSpecific_All_Female_2013_2014.txt", header = TRUE, row.names = NULL, as.is = TRUE, sep = "\t")
i = 10

stream <- parents[i, 1]
year <- parents[i, 2]
n.w <- parents[i, 3]
n.h <- parents[i, 4]

# number of independent simulations to do
trials <- 2000

# R will spit out a .txt file for each F1 prop, but you can change how many run in each instance to speed things along
sample_prop_off <- c(1/20, 1/10, 1/6, 1/3, 1/2, 2/3, 5/6, 1)[1:8]

rrs.values <- 0.5 # seq(from = 0.5, to = 1, by = 0.05)
mu.n.values <- c(0.25, seq(from = 0.5, to = 5, by = 0.5)) # R/female S (i.e. R/S / 2) # THESE SIMULATIONS ARE TOTALLY DEPENDENT ON THE MU AND SIZE OF THE NATURAL DISTRIBUTION OF R/S
size.values <- c(1, 2, 5, 10) #0.95

#ptm <- proc.time()
OUT <- matrix(data = NA, nrow = length(rrs.values)*length(mu.n.values)*length(size.values), ncol = 10, dimnames = list(seq(length(rrs.values)*length(mu.n.values)*length(size.values)), c("RRS", "mu.n", "size.n", "N.h.parents", "N.w.parents", "N.h.offspring", "N.w.offspring", "Perm Power", "nbGLM Power", "ttest Power")))


#### Stopped HERE ####

for(sample_prop_off in c(1/20, 1/10, 1/6, 1/3, 1/2, 2/3, 5/6, 1)){
  
  
  for (r in seq_along(rrs.values)) {
    errs <- rrs.values[r]
  
    for (m in seq_along(mu.n.values)) {
      mu.n <- mu.n.values[m]
      mu.h <- mu.n*rrs.values

      for(s in seq_along(size.values)){
        size.n <- size.values[s]
              
        OUT2 <- matrix(data = NA, nrow = trials, ncol = 10, dimnames = list(1:trials, c("errs", "mu.n", "size.n", "n.h", "n.w", "n.NHsampled", "n.NWsampled", "perm_1tail_pvalue", "nbGLM_1tail_pvalue", "ttest_1tail_pvalue")))

        for(i in 1:trials) {
            
          NWoffspring <- rnbinom(n = n.w, size = size.n, mu = mu.n)   #take a sample of offspring based on this distribution
          NHoffspring <- rnbinom(n = n.h, size = size.n, mu = mu.h)
          Noffspring <- sum(NWoffspring, NHoffspring) # how many offspring produced given mu.n and RRS
      
# figure out how to sample a proportion of the offspring, what will family sizes look like
      # Wild
          n.NWoffspring <- sum(NWoffspring) # number off WILD offspring produced
          n.NWsampled <- round(n.NWoffspring*sample_prop_off)
          NWsampled <- sort(sample(x = 1:n.NWoffspring, size = n.NWsampled, replace = FALSE)) # individual WILD offspring that were sampled (produced * sample rate)

          NWassigned <- rep(NA, n.w)
          Wdams.sampled <- rep(x = c(1:n.w), NWoffspring)[NWsampled]
          NWassigned[which(!seq(n.w) %in% Wdams.sampled)] <- 0
          NWassigned[as.numeric(names(table(Wdams.sampled)))] <- as.vector(table(Wdams.sampled))

      # Hatchery
          n.NHoffspring <- sum(NHoffspring) # number off HATCHERY offspring produced
          n.NHsampled <- round(n.NHoffspring*sample_prop_off)
          NHsampled <- sort(sample(x = 1:n.NHoffspring, size = n.NHsampled, replace = FALSE)) # individual HATCHERY offspring that were sampled (produced * sample rate)

          NHassigned <- rep(NA, n.h)
          Hdams.sampled <- rep(x = c(1:n.h), NHoffspring)[NHsampled]
          NHassigned[which(!seq(n.h) %in% Hdams.sampled)] <- 0
          NHassigned[as.numeric(names(table(Hdams.sampled)))] <- as.vector(table(Hdams.sampled))


#################################################################
# Do a one-sided test
# output should have t.test, permutation test, and Nbinom

          mydata <- data.frame(nOff = c(NWassigned, NHassigned), Origin = c(rep("W", n.w), rep("H", n.h)))

          if(sum(mydata$nOff) == 0) {
            next
          }
# permutation test (1-tail)
          perm_1tail_pvalue <- round(x = as.numeric(pvalue(oneway_test(nOff~Origin, data = mydata, distribution = approximate(B = 10000), alternative = "less"))), digits = 4)

        #reps <- 10000
        #true_diff <- mean(NWassigned)-mean(NHassigned)
        #Wfams <- replicate(reps, sum(sample(c(NWassigned, NHassigned), size = n.w, replace = FALSE)))
        #Hfams <- sum(NWassigned, NHassigned)-Wfams
        #rsW <- Wfams/n.w
        #rsH <- Hfams/n.h
        #diffs <- rsW-rsH
        #perm_1tail_pvalue <- round(sum(diffs> = true_diff)/reps, 4)

        #perm_1tail_pvalue2 <- round(as.numeric(pvalue(oneway_test(nOff~Origin, data = mydata, distribution = approximate(B = 10000), alternative = "less"))), 4)
# negative binomial GLM (1-tail)
          fit <- glm.nb(nOff~Origin, data = mydata, init.theta = 1, link = log)
          nbGLM_1tail_pvalue <- round(pnorm(summary(fit)$coefficients[2, 3], lower.tail = FALSE), 4)

# t.test (1-tail)
          test <- t.test(NWassigned, NHassigned, alternative = "greater")
          ttest_1tail_pvalue <- round(test$p.value, 4)

# building output
          OUT2[i, ] <- c(errs, mu.n, size.n, n.h, n.w, n.NHsampled, n.NWsampled, perm_1tail_pvalue, nbGLM_1tail_pvalue, ttest_1tail_pvalue)

        }

        OUT[sum((r-1)*length(mu.n.values)*length(size.values), (m-1)*length(size.values), s), ] <- c(OUT2[1, c(1:5)], mean(OUT2[, "n.NHsampled"], na.rm = TRUE), mean(OUT2[, "n.NWsampled"], na.rm = TRUE), sum(OUT2[, "perm_1tail_pvalue"] < 0.05, na.rm = TRUE) / sum(!is.na(OUT2[, "perm_1tail_pvalue"])), sum(OUT2[, "nbGLM_1tail_pvalue"] < 0.05, na.rm = TRUE) / sum(!is.na(OUT2[, "nbGLM_1tail_pvalue"])), sum(OUT2[, "ttest_1tail_pvalue"] < 0.05, na.rm = TRUE) / sum(!is.na(OUT2[, "ttest_1tail_pvalue"])))
  
      }

    }
  
  
  }

  write.table(OUT, paste("StreamSpecific_ConstSize", "/simulation_results_stream_", stream, "_", year, "_prop_", round(sample_prop_off, 3), "_trials_", trials, ".txt", sep = ''), col.names = TRUE, sep = "\t")

}
proc.time() - ptm; beep(8)