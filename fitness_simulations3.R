#======================================================================================================================#
# Script created by Mark Christie and Mike Ford (NOAA): contact at Redpath.Christie@gmail.com
# Script created in version R 3.0.1 in 2013
# This script:  Simulates data and examines the effect of p-values vs sample sizes for RRS studies
# Usage notes:  Run in chunks
#======================================================================================================================#
# Source files, import packages, set working directory, initialize variables
setwd("C:/output")

#library(coin)
#library(MASS)
#library(lattice)
#======================================================================================================================#

##### start of simulations ###########

rrs.values <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
n.spawners <- seq(from = 10, to = 400 , by = 10)
trials = 5000                                           # number of independent simulations to do
mu = 1               


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
