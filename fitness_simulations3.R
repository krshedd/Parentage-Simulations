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
#======================================================================================================================#
# Source files, import packages, set working directory, initialize variables

ls()
rm(list=ls(all=TRUE))
setwd("V:/WORK/Pink/AHRG/Parentage simulations/Mark Christie")

# Which run is this?
run=1


# Define stray rates to model
stray_rates=c(0.05,0.15,0.5)

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

mu.n = 2  # THESE SIMULATIONS ARE TOTALLY DEPENDENT ON THE MU AND SIZE OF THE NATURAL DISTRIBUTION OF R/S
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

write.table(OUT, paste("simulation_results_stray",stray,"_prop_",round(sample_prop_off,3),"_trials_",trials,"_muRSn_",mu.n,".txt",sep=''), col.names = TRUE, sep="\t")

#str(OUT)
#####################################################################################################################

# Create a function to plot results like Christie et al. 2014
# For testing
# file="mu2/simulation_results_stray0.5_prop_0.167_trials_2000_muRSn_2.txt";stray="0.50";prop="1/6";type="perm"
powerplot=function(file,stray,prop,type){
  setwd("V:/WORK/Pink/AHRG/Parentage simulations/Mark Christie")
  OUT=read.table(file)
  ifelse(type=="perm",assign("y",OUT[,6]),ifelse(type=="nbGLM",assign("y",OUT[,7]),ifelse(type=="ttest",assign("y",OUT[,8]),stop("HOSER!!! type must = 'perm', 'nbGLM', or 'ttest"))))
  ifelse(type=="perm",assign("maintxt","Permutation Test"),ifelse(type=="nbGLM",assign("maintxt","Negative Binomial GLM"),assign("maintxt","T-test")))
  par(mar=c(5.1,5.1,3.1,2.1))
  plot(0,xlim=c(0,3000),ylim=c(0,1.4),xlab="Total Number of Parents Sampled (Single Sex)",ylab="Power",main=maintxt,cex.lab=1.5,cex.main=2,xaxt="n",yaxt="n",pch=16,cex=2)
  abline(h=0.8,lwd=5)
  #points(y~apply(OUT,1,function(x){sum(x[2],x[3])}),pch=16,col=rep(c(1,2,3,4,5,6),each=30),cex=2)
  rrs=unique(OUT[,1])
  n.rrs=length(rrs)
  parnts=unique(apply(OUT,1,function(x){sum(x[2],x[3])}))
  n.parnts=length(parnts)
  x=c(0,parnts)
  for(i in 1:n.rrs){
    assign(paste("y",i,sep="_"),c(0,y[(1+(i-1)*n.parnts):(i*n.parnts)]))
    points(get(paste("y",i,sep="_"))~x,pch=16,cex=2,col=i)
    lines(get(paste("y",i,sep="_"))~x,lwd=3,col=i)
  }
  axis(1,cex.axis=1.5)
  axis(2,at=seq(from=0,to=1,by=0.2),cex.axis=1.5)
  legend("topleft",legend=paste("RRS = ",c("0.50","0.60","0.70","0.80","0.90","0.95"),sep=""),bty="n",col=c(1:6),pch=16,cex=1.2,y.intersp = 0.8)
  text(x=3000,y=1.4,paste("Stray rate = ",stray,sep=""),cex=1.5,pos=2)
  text(x=3000,y=1.3,paste("Proportion offspring =   ",prop,sep=""),cex=1.5,pos=2)
}

## MU.N=8
# 50% STRAY
powerplot(file="mu8/simulation_results_stray0.5_prop_0.167_trials_2000.txt",stray="0.50",prop="1/6",type="perm")
powerplot(file="mu8/simulation_results_stray0.5_prop_0.333_trials_2000.txt",stray="0.50",prop="2/6",type="perm")
powerplot(file="mu8/simulation_results_stray0.5_prop_0.5_trials_2000.txt",stray="0.50",prop="3/6",type="perm")
powerplot(file="mu8/simulation_results_stray0.5_prop_0.667_trials_2000.txt",stray="0.50",prop="4/6",type="perm")
powerplot(file="mu8/simulation_results_stray0.5_prop_0.833_trials_2000.txt",stray="0.50",prop="5/6",type="perm")
powerplot(file="mu8/simulation_results_stray0.5_prop_1_trials_2000.txt",stray="0.50",prop="6/6",type="perm")

# 15% STRAY
powerplot(file="mu8/simulation_results_stray0.15_prop_0.167_trials_2000.txt",stray="0.15",prop="1/6",type="perm")
powerplot(file="mu8/simulation_results_stray0.15_prop_0.333_trials_2000.txt",stray="0.15",prop="2/6",type="perm")
powerplot(file="mu8/simulation_results_stray0.15_prop_0.5_trials_2000.txt",stray="0.15",prop="3/6",type="perm")
powerplot(file="mu8/simulation_results_stray0.15_prop_0.667_trials_2000.txt",stray="0.15",prop="4/6",type="perm")
powerplot(file="mu8/simulation_results_stray0.15_prop_0.833_trials_2000.txt",stray="0.15",prop="5/6",type="perm")
powerplot(file="mu8/simulation_results_stray0.15_prop_1_trials_2000.txt",stray="0.15",prop="6/6",type="perm")

# 5% STRAY
powerplot(file="mu8/simulation_results_stray0.05_prop_0.167_trials_2000.txt",stray="0.05",prop="1/6",type="perm")
powerplot(file="mu8/simulation_results_stray0.05_prop_0.333_trials_2000.txt",stray="0.05",prop="2/6",type="perm")
powerplot(file="mu8/simulation_results_stray0.05_prop_0.5_trials_2000.txt",stray="0.05",prop="3/6",type="perm")
powerplot(file="mu8/simulation_results_stray0.05_prop_0.667_trials_2000.txt",stray="0.05",prop="4/6",type="perm")
powerplot(file="mu8/simulation_results_stray0.05_prop_0.833_trials_2000.txt",stray="0.05",prop="5/6",type="perm")
powerplot(file="mu8/simulation_results_stray0.05_prop_1_trials_2000.txt",stray="0.05",prop="6/6",type="perm")


## MU.N=2
# 50% STRAY
powerplot(file="mu2/simulation_results_stray0.5_prop_0.167_trials_2000_muRSn_2.txt",stray="0.50",prop="1/6",type="perm")
powerplot(file="mu2/simulation_results_stray0.5_prop_0.333_trials_2000_muRSn_2.txt",stray="0.50",prop="2/6",type="perm")
powerplot(file="mu2/simulation_results_stray0.5_prop_0.5_trials_2000_muRSn_2.txt",stray="0.50",prop="3/6",type="perm")
powerplot(file="mu2/simulation_results_stray0.5_prop_0.667_trials_2000_muRSn_2.txt",stray="0.50",prop="4/6",type="perm")
powerplot(file="mu2/simulation_results_stray0.5_prop_0.833_trials_2000_muRSn_2.txt",stray="0.50",prop="5/6",type="perm")
powerplot(file="mu2/simulation_results_stray0.5_prop_1_trials_2000_muRSn_2.txt",stray="0.50",prop="6/6",type="perm")

# 15% STRAY
powerplot(file="mu2/simulation_results_stray0.15_prop_0.167_trials_2000_muRSn_2.txt",stray="0.15",prop="1/6",type="perm")
powerplot(file="mu2/simulation_results_stray0.15_prop_0.333_trials_2000_muRSn_2.txt",stray="0.15",prop="2/6",type="perm")
powerplot(file="mu2/simulation_results_stray0.15_prop_0.5_trials_2000_muRSn_2.txt",stray="0.15",prop="3/6",type="perm")
powerplot(file="mu2/simulation_results_stray0.15_prop_0.667_trials_2000_muRSn_2.txt",stray="0.15",prop="4/6",type="perm")
powerplot(file="mu2/simulation_results_stray0.15_prop_0.833_trials_2000_muRSn_2.txt",stray="0.15",prop="5/6",type="perm")
powerplot(file="mu2/simulation_results_stray0.15_prop_1_trials_2000_muRSn_2.txt",stray="0.15",prop="6/6",type="perm")

## MU.N=1
# 50% STRAY
powerplot(file="mu1/simulation_results_stray0.5_prop_0.167_trials_2000_muRSn_1.txt",stray="0.50",prop="1/6",type="perm")
powerplot(file="mu1/simulation_results_stray0.5_prop_0.333_trials_2000_muRSn_1.txt",stray="0.50",prop="2/6",type="perm")
powerplot(file="mu1/simulation_results_stray0.5_prop_0.5_trials_2000_muRSn_1.txt",stray="0.50",prop="3/6",type="perm")
powerplot(file="mu1/simulation_results_stray0.5_prop_0.667_trials_2000_muRSn_1.txt",stray="0.50",prop="4/6",type="perm")
powerplot(file="mu1/simulation_results_stray0.5_prop_0.833_trials_2000_muRSn_1.txt",stray="0.50",prop="5/6",type="perm")
powerplot(file="mu1/simulation_results_stray0.5_prop_1_trials_2000_muRSn_1.txt",stray="0.50",prop="6/6",type="perm")

# 15% STRAY
powerplot(file="mu1/simulation_results_stray0.15_prop_0.167_trials_2000_muRSn_1.txt",stray="0.15",prop="1/6",type="perm")
powerplot(file="mu1/simulation_results_stray0.15_prop_0.333_trials_2000_muRSn_1.txt",stray="0.15",prop="2/6",type="perm")
powerplot(file="mu1/simulation_results_stray0.15_prop_0.5_trials_2000_muRSn_1.txt",stray="0.15",prop="3/6",type="perm")
powerplot(file="mu1/simulation_results_stray0.15_prop_0.667_trials_2000_muRSn_1.txt",stray="0.15",prop="4/6",type="perm")
powerplot(file="mu1/simulation_results_stray0.15_prop_0.833_trials_2000_muRSn_1.txt",stray="0.15",prop="5/6",type="perm")
powerplot(file="mu1/simulation_results_stray0.15_prop_1_trials_2000_muRSn_1.txt",stray="0.15",prop="6/6",type="perm")
#####################################################################################################################

# LOESS DOES A BAD JOB ESTIMATING X, see below, unless span is really low (0.1 vs. default of 0.75)
power_level=0.8
type="perm"

setwd("V:/WORK/Pink/AHRG/Parentage simulations/Mark Christie")
filenames <- list.files("mu2", pattern="*.txt", full.names=TRUE)
OUT=NULL
OUT2=NULL
for(f in 1:length(filenames)){
  OUT=read.table(filenames[f])
  ifelse(type=="perm",assign("y",OUT[,6]),ifelse(type=="nbGLM",assign("y",OUT[,7]),ifelse(type=="ttest",assign("y",OUT[,8]),stop("HOSER!!! type must = 'perm', 'nbGLM', or 'ttest"))))
  rrs=unique(OUT[,1])
  n.rrs=length(rrs)
  parnts=unique(apply(OUT,1,function(x){sum(x[2],x[3])}))
  n.parnts=length(parnts)
  x=c(0,parnts)
  newdat=seq(0,max(parnts),by=10)

  for(i in 1:n.rrs){
    assign(paste("y",i,sep="_"),c(0,y[(1+(i-1)*n.parnts):(i*n.parnts)]))
    fit=loess(get(paste("y",i,sep="_"))~x,span=ifelse(rrs[i]<0.7,0.1,0.3),weights=c(5,rep(1,length(x)-1)))
    OUT2=rbind(OUT2,cbind("mu"=as.numeric(unlist(strsplit(strsplit(filenames[f],"muRSn_")[[1]][2],".txt"))[1]),
                          "stray"=as.numeric(unlist(strsplit(strsplit(filenames[f],"stray")[[1]][2],"_prop"))[1]),
                          "pF1"=as.numeric(unlist(strsplit(strsplit(filenames[f],"prop_")[[1]][2],"_trials"))[1]),
                          "rrs"=rrs[i],
                          "nF0"=ifelse(min(abs(predict(fit,newdata=newdat)-power_level),na.rm=TRUE)>0.02,NA,newdat[which.min(abs(predict(fit,newdata=newdat)-power_level))])
    ))
  }
} 

# ,weights=c(5,rep(1,length(x)-1))

rm(list=setdiff(ls(), c("OUT2","filenames")))

data=OUT2[which(OUT2[,"stray"]==0.50),]

F0s=sort(unique(data[,5]))
F1s=sort(unique(data[,3]))

mat=matrix(nrow=length(F1s),ncol=length(F0s),dimnames=list(F1s,F0s))

for(i in 1:length(F1s)){
  for(j in 1:length(F0s)){
    mat[i,j]=ifelse(length(which(OUT2[,5]==F0s[j] & OUT2[,3]==F1s[i]))==0,NA,pmin(OUT2[which(OUT2[,5]==F0s[j] & OUT2[,3]==F1s[i]),"rrs"]))
  }
}


# not working, need to deal with NAs
filled.contour(x=as.numeric(colnames(mat)),y=as.numeric(rownames(mat)),z=t(mat))



#####################################################################################################################
power_level=0.8
y=c(0.000,0.900,0.975,0.992,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000)
x=seq(0,3000,by=100)
fit=loess(y~x,span=0.2,weights=c(5,rep(1,length(x)-1)))
plot(y~x)
lines(predict(fit)~x, col='red', lwd=2)
newdat[which.min(abs(predict(fit,newdata=newdat)-power_level))]
abline(v=newdat[which.min(abs(predict(fit,newdata=newdat)-power_level))],lwd=3)
#####################################################################################################################







# How to determine the # of total parents (hatchery + wild of a single sex) to sample for power=0.8
power_level=0.8
x=seq(0,3000,by=100)
y_4=c(0.0000,0.0285,0.0595,0.0975,0.1415,0.2065,0.2465,0.2800,0.3485,0.4200,0.4500,0.5190,0.5585,0.6190,0.6450,0.6835,0.7270,0.7425,0.7705,0.8025,0.8355,0.8455,0.8570,0.8990,0.9070,0.9225,0.9280,0.9370,0.9500,0.9475,0.9710)
fit=loess(y_4~x, weights=c(5,rep(1,length(x)-1)))
plot(y_4~x)
lines(predict(fit)~x, col='red', lwd=2)
newdat=seq(0,3000,by=10)
x_power=newdat[which.min(abs(predict(fit,newdata=newdat)-power_level))]
x_power
abline(v=x_power,lwd=3)
abline(h=power_level,lwd=3)


#3D plot examples
x1=rnorm(100);x2=rpois(100,lambda=4)
logist=function(x1,x2){
y=1/(1+exp((???1)*(???3+0.6*x1+.5*x2)))}
par(bg="white")
x1r=range(x1);x1seq=seq(x1r[1],x1r[2],
length=30)
x2r=range(x2);x2seq=seq(x2r[1],x2r[2],
length=30)
z=outer(x1seq,x2seq,logist)
persp(x=x1seq,y=x2seq,z=z,theta=???30,
zlim=c(???0.2,1.2),col=2)


library(scatterplot3d)
x=rnorm(80);y=rpois(80,l=7)
z=3+1.1*x+0.4*y+15*rnorm(80)
s3d=scatterplot3d(x,y,z)



x1=x2=seq(???10,10,length=51)
dens=matrix(dmvnorm(expand.grid(x1,x2),sigma
=rbind(c(3,2),c(2,3))),ncol=length(x1))


# Create a function to plot results in 3D
# For testing
# mu="2";stray="0.5";type="perm"
threedpowerplot=function(mu,stray,rrs,type){
  setwd("V:/WORK/Pink/AHRG/Parentage simulations/Mark Christie")
  z=NULL
  props=c("0.167","0.333","0.5","0.667","0.833","1")
  for(prop in props){
    OUT=read.table(file=paste("mu",mu,"/simulation_results_stray",stray,"_prop_",prop,"_trials_2000_muRSn_",mu,".txt",sep=''))
    ifelse(type=="perm",assign("y",OUT[which(OUT[,1]==rrs),6]),ifelse(type=="nbGLM",assign("y",OUT[which(OUT[,1]==rrs),7]),ifelse(type=="ttest",assign("y",OUT[which(OUT[,1]==rrs),8]),stop("HOSER!!! type must = 'perm', 'nbGLM', or 'ttest"))))
    ifelse(type=="perm",assign("maintxt","Permutation Test"),ifelse(type=="nbGLM",assign("maintxt","Negative Binomial GLM"),assign("maintxt","T-test")))
    y=c(0,y)
    z=c(z,y)
  }
  x=rep(seq(0,3000,by=100),6)
  y=rep(seq(1/6,1,by=1/6),each=31)
  
  s3d=scatterplot3d(x,y,z,pch=16,type="h",scale.y=1,angle=120,cex.symbols=2,zlim=c(0.8,1),xlim=c(0,3000),box=FALSE,
                    xlab="Total # Parents (Single Sex)",ylab="Proportion of Offspring Sampled)",zlab="Power")
}

threedpowerplot(mu=1,stray=0.50,rrs=0.50,type="perm")


  dens=NULL
  for(i in 1:6){
    dens=rbind(dens,z[(1+(i-1)*31):(i*31)])
  }
  
  x1=seq(0,3000,by=100)
  y1=seq(1/6,1,by=1/6)
  for(i in length(x1):1){
    s3d$points3d(rep(x1[i],length(y1)),y1,dens[,i],type="l")
  }
  
  for(i in length(y1):1){
    s3d$points3d(x1,rep(y1[i],length(x1)),dens[i,],type="l")
  }
  
  s3d$plane3d(Intercept=c(0.8,0,0))


#####################################################################################################################
# Which has the most power?
perm.nbGLM=apply(OUT,1,function(row){row[6]-row[7]})
hist(perm.nbGLM)
summary(perm.nbGLM)
# So nbGLM has more power than permutation test ALWAYS
perm.ttest=apply(OUT,1,function(row){row[6]-row[8]})
hist(perm.ttest)
summary(perm.ttest)
# And t-test has more power than permutation test pretty much always
nbGLM.ttest=apply(OUT,1,function(row){row[7]-row[8]})
hist(nbGLM.ttest)
summary(nbGLM.ttest)
# And nbGLM has more power than the t-test, which makes sense
apply(OUT,1,function(row){max(row[6:8])})

##### start of simulations to determine RRS point estimate 95% CI ####################################################################################################





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

##### start of simulations ####################################################################################################
ptm <- proc.time()

errs = 0.8

# Number of parents sampled
n.spawners <- seq(from = 100, to = 3000 , by = 100)

# Define stray rates to model
stray_rates=c(0.05,0.15,0.5)

# Define sample proportion of offspring
sample_prop_off_rates=c((1:6)/6)

mu.n = 8  # THESE SIMULATIONS ARE TOTALLY DEPENDENT ON THE MU AND SIZE OF THE NATURAL DISTRIBUTION OF R/S, 8 was original from Christie et al. 2014
size.n = 0.95
var.n = nbnom_variance(mu.n,size.n)

OUT2=NULL

mu.h = mu.n*errs
size.h =  size2(mu.n,size.n,mu.h)

for(stray in stray_rates){
  for (z in 1:length(n.spawners)) {
    n.h = n.spawners[z]*stray  # number of hatchery F1 spawners
    n.w = n.spawners[z]*(1-stray)  # number of wild F1 spanwers
  
    NWoffspring=rnbinom(n=n.w,size=size.n,mu=mu.n)   #take a sample of offspring based on this distribution
    NHoffspring=NWoffspring*errs
    
    #################################################################
    # Do a one-sided test
    # output should have t.test, permutation test, and Nbinom
    
    mydata=data.frame(nOff=c(NWoffspring,NHoffspring),Origin=c(rep("W",n.w),rep("H",n.h)))
    
    # negative binomial GLM (1-tail)
    fit <- glm.nb(nOff~Origin, data=mydata, init.theta=1, link=log)
    nbGLM_RRS=round(1/exp(summary(fit)$coefficients[2,1]),4)
    nbGLM_5CI=round(1/exp(summary(fit)$coefficients[2,1]+1.645*summary(fit)$coefficients[2,2]),4)
    nbGLM_95CI=round(1/exp(summary(fit)$coefficients[2,1]-1.645*summary(fit)$coefficients[2,2]),4)
  
    # building output
    out <- cbind(errs, n.h, n.w, n.NHsampled, n.NWsampled, nbGLM_RRS, nbGLM_5CI, nbGLM_95CI)
    OUT2 <- rbind(OUT2, out)
  }
}
    
    
proc.time() - ptm

rownames(OUT)=c(1:dim(OUT)[1])
colnames(OUT)=c("RRS","N.h.parents","N.w.parents","N.h.offspring","N.w.offspring","Perm Power","nbGLM Power","ttest Power")

write.table(OUT, paste("simulation_results_stray",stray,"_prop_",round(sample_prop_off,3),"_trials_",trials,"_muRSn_",mu.n,".txt",sep=''), col.names = TRUE, sep="\t")
