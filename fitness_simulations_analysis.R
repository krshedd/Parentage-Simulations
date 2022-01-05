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
#setwd("F:/2015 AHRP April Seattle")
source("V:/WORK/Pink/AHRG/Parentage simulations/R Functions/HW_Functions.R")    # Kyle's functions for parentage analysis
#source("F:/2015 AHRP April Seattle/R Functions/HW_Functions.R")
opar=par()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Create a function to plot results like Christie et al. 2014 ####
# For testing
# file="mu2/simulation_results_stray0.5_prop_0.167_trials_2000_muRSn_2.txt";stray="0.50";prop="1/6";type="perm"
powerplot=function(file,stray,prop,type){
  setwd("V:/WORK/Pink/AHRG/Parentage simulations/Mark Christie")
  xmax=3000
  power=0.8
  mu.fam=as.numeric(unlist(strsplit(strsplit(file,"mu")[[1]][2],"/si"))[1])
  OUT=read.table(file)
  ifelse(type=="perm",assign("y",OUT[,6]),ifelse(type=="nbGLM",assign("y",OUT[,7]),ifelse(type=="ttest",assign("y",OUT[,8]),stop("HOSER!!! type must = 'perm', 'nbGLM', or 'ttest"))))
  ifelse(type=="perm",assign("maintxt","Permutation Test"),ifelse(type=="nbGLM",assign("maintxt","Negative Binomial GLM"),assign("maintxt","T-test")))
  par(mar=c(5.1,5.1,3.1,2.1))
  par(bg=colors()[c(356)])
  plot(0,xlim=c(0,xmax),ylim=c(0,1.4),xlab="Dams sampled",ylab="Power",main=maintxt,cex.lab=1.5,cex.main=2,xaxt="n",yaxt="n",pch=16,cex=2)
  #xlab="Total Number of Parents Sampled (Single Sex)"
  abline(h=power,lwd=5)
  #points(y~apply(OUT,1,function(x){sum(x[2],x[3])}),pch=16,col=rep(c(1,2,3,4,5,6),each=30),cex=2)
  rrs=unique(OUT[,1])
  n.rrs=length(rrs)
  parnts=unique(apply(OUT,1,function(x){sum(x[2],x[3])}))
  n.parnts=length(parnts)
  x=c(0,parnts)
  #newx=seq(0,max(parnts),by=10)
  for(i in 1:n.rrs){
    assign(paste("y",i,sep="_"),c(0,y[(1+(i-1)*n.parnts):(i*n.parnts)]))
    points(get(paste("y",i,sep="_"))~x,pch=16,cex=2,col=i)
    #lines(get(paste("y",i,sep="_"))~x,lwd=3,col=i)
    fit=loess(get(paste("y",i,sep="_"))~x,span=ifelse(rrs[i]<0.7,0.1,0.3),weights=c(5,rep(1,length(x)-1))) #<0.7,0.1,0.3
    lines(predict(fit,x)~x,col=i,lwd=3)
    #lines(predict(fit,newx)~newx,col=i,lwd=3)
  }
  axis(1,cex.axis=1.5)
  axis(2,at=seq(from=0,to=1,by=0.2),cex.axis=1.5,las=2)
  legend("topleft",legend=paste("RRS = ",rrs,sep=""),bty="n",col=c(1:n.rrs),pch=16,cex=1,y.intersp = 0.8) # cex=1.2
  text(x=xmax,y=1.4,paste("Stray rate = ",stray,sep=""),cex=1.5,pos=2)
  text(x=xmax,y=1.3,paste("Proportion offspring =   ",prop,sep=""),cex=1.5,pos=2)
  text(x=xmax,y=1.2,substitute(paste(mu," = ",m,sep=''),list(m=mu.fam)),cex=1.5,pos=2)
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

# Newer
powerplot(file="mu2/simulation_results_stray0.1_prop_0.05_trials_2000_muRSn_2_RRS_0.01.txt",stray="0.10",prop="1/20",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.05_trials_2000_muRSn_2_RRS_0.05.txt",stray="0.10",prop="1/20",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.05_trials_2000_muRSn_2_RRS_0.1.txt",stray="0.10",prop="1/20",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.05_trials_2000_muRSn_2_RRS_0.2.txt",stray="0.10",prop="1/20",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.05_trials_2000_muRSn_2_RRS_0.3.txt",stray="0.10",prop="1/20",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.05_trials_2000_muRSn_2_RRS_0.4.txt",stray="0.10",prop="1/20",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.05_trials_2000_muRSn_2_RRS_0.5.txt",stray="0.10",prop="1/20",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.05_trials_2000_muRSn_2_RRS_0.55.txt",stray="0.10",prop="1/20",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.05_trials_2000_muRSn_2_RRS_0.6.txt",stray="0.10",prop="1/20",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.05_trials_2000_muRSn_2_RRS_0.62.txt",stray="0.10",prop="1/20",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.05_trials_2000_muRSn_2_RRS_0.64.txt",stray="0.10",prop="1/20",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.05_trials_2000_muRSn_2_RRS_0.66.txt",stray="0.10",prop="1/20",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.05_trials_2000_muRSn_2_RRS_0.68.txt",stray="0.10",prop="1/20",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.05_trials_2000_muRSn_2_RRS_0.7.txt",stray="0.10",prop="1/20",type="perm")

powerplot(file="mu2/simulation_results_stray0.1_prop_0.1_trials_2000_muRSn_2.txt",stray="0.10",prop="1/10",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.167_trials_2000_muRSn_2.txt",stray="0.10",prop="1/6",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.333_trials_2000_muRSn_2.txt",stray="0.10",prop="1/3",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.5_trials_2000_muRSn_2.txt",stray="0.10",prop="1/2",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.667_trials_2000_muRSn_2.txt",stray="0.10",prop="2/3",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.833_trials_2000_muRSn_2.txt",stray="0.10",prop="5/6",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_1_trials_2000_muRSn_2.txt",stray="0.10",prop="1",type="perm")

powerplot(file="mu2/simulation_results_stray0.1_prop_0.167_trials_2000_muRSn_2_RRS_0.72.txt",stray="0.10",prop="1/6",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.333_trials_2000_muRSn_2_RRS_0.72.txt",stray="0.10",prop="1/3",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.333_trials_2000_muRSn_2_RRS_0.75.txt",stray="0.10",prop="1/3",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.333_trials_2000_muRSn_2_RRS_0.76.txt",stray="0.10",prop="1/3",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.333_trials_2000_muRSn_2_RRS_0.77.txt",stray="0.10",prop="1/3",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.5_trials_2000_muRSn_2_RRS_0.75.txt",stray="0.10",prop="1/2",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.5_trials_2000_muRSn_2_RRS_0.77.txt",stray="0.10",prop="1/2",type="perm")
powerplot(file="mu2/simulation_results_stray0.1_prop_0.5_trials_2000_muRSn_2_RRS_0.78.txt",stray="0.10",prop="1/2",type="perm")

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Find intercept at 0.8 power ####
# Add line segment to show how many F0 need to be sampled...

file="mu2/simulation_results_stray0.5_prop_0.667_trials_2000_muRSn_2.txt";stray="0.50";prop="6/6";type="perm"

setwd("V:/WORK/Pink/AHRG/Parentage simulations/Mark Christie")
xmax=3000
power=0.8
mu.fam=as.numeric(unlist(strsplit(strsplit(file,"mu")[[1]][2],"/si"))[1])
OUT=read.table(file)
ifelse(type=="perm",assign("y",OUT[,6]),ifelse(type=="nbGLM",assign("y",OUT[,7]),ifelse(type=="ttest",assign("y",OUT[,8]),stop("HOSER!!! type must = 'perm', 'nbGLM', or 'ttest"))))
ifelse(type=="perm",assign("maintxt","Permutation Test"),ifelse(type=="nbGLM",assign("maintxt","Negative Binomial GLM"),assign("maintxt","T-test")))
rrs=unique(OUT[,1])
n.rrs=length(rrs)
parnts=unique(apply(OUT,1,function(x){sum(x[2],x[3])}))
n.parnts=length(parnts)
x=c(0,parnts)
newx=seq(0,max(parnts),by=10)

i=4
assign(paste("y",i,sep="_"),c(0,y[(1+(i-1)*n.parnts):(i*n.parnts)]))
fit=loess(get(paste("y",i,sep="_"))~x,span=ifelse(rrs[i]<0.7,0.1,0.3),weights=c(5,rep(1,length(x)-1))) #<0.7,0.1,0.3

power_level=0.8

npars=newx[which.min(abs(predict(fit,newdata=newx)-power_level))]

segments(npars,-0.5,npars,0.785,lwd=5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create matrix for color contour ####
# LOESS DOES A BAD JOB ESTIMATING X, see below, unless span is really low (0.1 vs. default of 0.75)
opar=par()



par(opar)
power_level=0.80
type="perm"
stray=0.1 # CHANGE THIS ####
mu.fam=2

setwd("V:/WORK/Pink/AHRG/Parentage simulations/Mark Christie")
filenames <- list.files(paste("mu",mu.fam,sep=""), pattern=paste("stray",stray,"_",sep=""), full.names=TRUE)
OUT=NULL
OUT2=NULL
for(f in 1:length(filenames)){
  OUT=read.table(filenames[f])
  ifelse(type=="perm",assign("y",OUT[,6]),ifelse(type=="nbGLM",assign("y",OUT[,7]),ifelse(type=="ttest",assign("y",OUT[,8]),stop("HOSER!!! type must = 'perm', 'nbGLM', or 'ttest"))))
  #if (y==0) next
  rrs=unique(OUT[,1])
  n.rrs=length(rrs)
  parnts=unique(apply(OUT,1,function(x){sum(x[2],x[3])}))
  n.parnts=length(parnts)
  x=c(0,parnts)
  newdat=seq(0,max(parnts),by=10)

  for(i in 1:n.rrs){
    assign(paste("y",i,sep="_"),c(0,y[(1+(i-1)*n.parnts):(i*n.parnts)]))
    fit=loess(get(paste("y",i,sep="_"))~x,span=ifelse(rrs[i]<0.7,0.1,0.3),weights=c(5,rep(1,length(x)-1))) #<0.7,0.1,0.3
    OUT2=rbind(OUT2,cbind("mu"=as.numeric(unlist(strsplit(strsplit(filenames[f],"muRSn_")[[1]][2],""))[1]),
                          "stray"=as.numeric(unlist(strsplit(strsplit(filenames[f],"stray")[[1]][2],"_prop"))[1]),
                          "pF1"=as.numeric(unlist(strsplit(strsplit(filenames[f],"prop_")[[1]][2],"_trials"))[1]),
                          "rrs"=rrs[i],
                          "nF0"=ifelse(min(abs(predict(fit,newdata=newdat)-power_level),na.rm=TRUE)>0.02,NA,newdat[which.min(abs(predict(fit,newdata=newdat)-power_level))])
    ))
  }
} 

# ,weights=c(5,rep(1,length(x)-1))

rm(list=setdiff(ls(), c("OUT2","filenames","stray","power_level","type","opar","mu.fam")))

data=OUT2[which(OUT2[,"stray"]==stray),]

F0s=sort(unique(data[,5]))
F1s=sort(unique(data[,3]))

mat=matrix(nrow=length(F1s),ncol=length(F0s),dimnames=list(F1s,F0s))

for(i in 1:length(F1s)){
  for(j in 1:length(F0s)){
    mat[i,j]=ifelse(length(which(data[,5]==F0s[j] & data[,3]==F1s[i]))==0,NA,pmin(data[which(data[,5]==F0s[j] & data[,3]==F1s[i]),"rrs"]))
  }
}

t(mat)

newmat=NULL
invisible(ifelse(stray==0.5, assign("i_span",c(0.3,0.4,0.3,0.3,0.34,0.3,0.25,0.4)),ifelse(stray==0.15, assign("i_span",c(0.6,0.5,0.6,0.4,0.3,0.35,0.3,0.3)), ifelse(stray==0.1, assign("i_span",c(0.6,0.7,0.7,0.6,0.7,0.5,0.5,0.5)), "HOSER: NEED TO SPECIFY i_span for each STRAY!!!"))))
for(i in 1:ncol(t(mat))){
  plot(t(mat)[,i]~rownames(t(mat)),xlim=c(0,3000),ylim=c(0,1),pch=16,xlab="F0s",ylab="RRS with power = 0.8",cex=2,main=paste("p F1 sampled",colnames(t(mat)))[i])
  #fit=lm(t(mat)[,i]~poly(as.numeric(rownames(t(mat))),4,raw=TRUE)) # poly isn't a good fit
  fit=loess(t(mat)[,i]~as.numeric(rownames(t(mat))),span=i_span[i])
  newmat=cbind(newmat,predict(fit,data.frame(x=as.numeric(rownames(t(mat))))))
  lines(as.numeric(rownames(t(mat))),predict(fit,data.frame(x=as.numeric(rownames(t(mat))))),col=4,lwd=3)
}
dimnames(newmat)=dimnames(t(mat))

#### Refining i_span ####
#i=5
#plot(t(mat)[,i]~rownames(t(mat)),xlim=c(0,3000),ylim=c(0,1),pch=16,xlab="F0s",ylab="RRS with power = 0.8",cex=2,main=paste("p F1 sampled",colnames(t(mat)))[i])
#fit=loess(t(mat)[,i]~as.numeric(rownames(t(mat))),span=0.8)
#lines(as.numeric(rownames(t(mat))),predict(fit,data.frame(x=as.numeric(rownames(t(mat))))),col=4,lwd=3)

# not working, need to deal with NAs
# created "newmat" that interpolates between data points using loess, this does much better!

#filled.contour(x=as.numeric(rownames(newmat)),y=as.numeric(colnames(newmat)),z=newmat)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Create color contour ####
## Making color contour plot of power based on # of parents (sex) sampled and F1 sampling proportion
## This is dependent on the stray rate and mean family size (i.e. the mu of the negative binomial for family size)

# png(file=paste("Figures/Single parent RRS mu_",mu.fam,"_stray_",stray,"_power_",power_level,".png",sep=""),width=724,height=654)

par(opar)
par(mar=c(5.1,6.1,3.1,6.1))
par(bg=colors()[c(356)])
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
zmin=0.5 #min(newmat,na.rm=TRUE) #0.50
zmax=1

filled.contour3(as.numeric(rownames(newmat)),y=as.numeric(colnames(newmat)),newmat,color=jet.colors,
                xlim=c(0,1500),ylim=c(0,1),zlim=c(zmin,zmax),xlab="",ylab="",cex.lab=2,
                axes=FALSE)
mtext(expression(paste("F"[1]," sampling proportion",sep="")),side=2,line=3.5,cex=2,las=0)
mtext(side=1,"Dams sampled",cex=2,line=3)
plot.axes = { axis(1, at=seq(0,3000,by=500), labels=seq(0,3000,by=500),cex.axis=1.5);
              axis(2, at=c(0,as.numeric(colnames(newmat))),labels=c("0","1/20","1/10","1/6","1/3","1/2","2/3","5/6","1"),las=1,cex.axis=1.5) }

#mtext(side=3,paste("Stray rate = ",stray,"; mu = ",mu.fam,"; Power = ",power_level,sep=""),cex=2,line=1)
mtext(side=3,substitute(paste("Stray rate = ",s,"; ",mu," = ",m,"; Power = ",p,sep=""),list(s=stray,m=mu.fam,p=power_level)),cex=2,line=1)
##methods
#points(1080,2/3,pch=16,cex=2);segments(0,2/3,1080,2/3,lwd=5);segments(1080,0,1080,2/3,lwd=5);text(1080,0.025,"1080",cex=2,pos=4);text(1080,2/3,"RRS = 0.8",cex=2,pos=3)
##how much power?
nsamp=637;year=2013;segments(nsamp/2,0,nsamp/2,1,lwd=5);text(nsamp/2,0.2,year,cex=3,pos=ifelse(nsamp<1000,4,2))
nsamp=260;year=2014;segments(nsamp/2,0,nsamp/2,1,lwd=5);text(nsamp/2,0.6,year,cex=3,pos=ifelse(nsamp<1500,4,2))

#for(i in 1:length(colnames(newmat))){
  #x=as.numeric(names(t(mat)[,i])[!is.na(t(mat)[,i])])
  #y=rep(as.numeric(colnames(t(mat))[i]),length(x))
  #points(x,y,pch=16,cex=2)
#}

#for(i in 1:length(colnames(newmat))){
  #for(j in 1:length(rownames(newmat))){
    #text(as.numeric(rownames(newmat)[j]),(as.numeric(colnames(newmat)[i])+0.01),ifelse(is.na(t(mat)[j,i])==TRUE,"",paste(t(mat)[j,i])),cex=1,srt=45,pos=3)
  #}
#}

zmin=0.5
zmax=1

par(new = "TRUE",plt = c(0.87,0.92,0.25,0.85),las = 1,cex.axis = 1.5)
filled.legend(as.numeric(rownames(newmat)),y=as.numeric(colnames(newmat)),newmat,color=jet.colors,xlab="",ylab="",xlim=c(min(xintercepts),max(xintercepts)),ylim=c(min(slopes),max(slopes)),zlim=c(zmin,zmax))
filled.legend(as.numeric(rownames(newmat)),y=as.numeric(colnames(newmat)),newmat,color=jet.colors,xlab="",ylab="",xlim=c(min(xintercepts),max(xintercepts)),ylim=c(min(slopes),max(slopes)),zlim=c(zmin,zmax))
mtext("RRS",cex=2,padj=-0.3)

# dev.off()

## change zmin to 0.5 to make comparable to previous figures (i.e. Tech Doc 5)
## also change xlim=c(0,1500) (1500 dams = 3000 total fish sampled) to be comparable to previous figures (i.e. Tech Doc 5)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# How much done? ####
setwd("V:/WORK/Pink/AHRG/Parentage simulations/Mark Christie/StreamSpecific_ConstSize_RRS_0.8")  # size.n = size.h # _RRS_0.8/RRS_0.5_Accident
getwd()
pats <- unlist(lapply(sort(c("Admiralty", "Erb", "Fish", "Gilmour", "HoganBay", "Paddy", "Prospect", "Sawmill", "Spring", "Stockdale")), function(stream) paste(stream, c(2013, 2014), sep="_")))
done <- setNames(object = lapply(pats, function(lst) list.files(path = getwd(), pattern = lst)), nm = pats)
lapply(done, length)

#### Stream-specific ####
# Visualize stream-specific power under different distributions of RS 
# (assuming that W and H  have the same variance, just different means)!!!!

# 3D scatterplot ####
library(scatterplot3d)
setwd("V:/WORK/Pink/AHRG/Parentage simulations/Mark Christie/StreamSpecific")  # var.n = var.h
setwd("V:/WORK/Pink/AHRG/Parentage simulations/Mark Christie/StreamSpecific_ConstSize")  # size.n = size.h

setwd("V:/Analysis/5_Coastwide/Multispecies/Alaska Hatchery Research Program/Parentage Simulations/Simulations based on Christie et al. 2015 Box 2/StreamSpecific_ConstSize")

stream=c("Spring", "Stockdale", "Gilmour", "HoganBay", "Erb", "Paddy", "Fish", "Admiralty", "Prospect", "Sawmill")[7]
year=2017

par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mfrow = c(1, 1))

for(prop in c(0.05, 0.1, 0.167, 0.333, 0.5, 0.667, 0.833, 1)){
  out <- read.table(file=grep(pattern = paste(year, "_prop_", prop, "_", sep=""), x = list.files(path=getwd(), pattern=stream, full.names=TRUE), value = TRUE))
  
  s3d <- scatterplot3d(out[,"mu.n"], out[,"size.n"], out[,"Perm.Power"], pch=16, type="h", scale.y=1, angle=120, cex.symbols=2,
                       zlim=c(0,1), xlim=c(0,5), box=FALSE, xlab="mean Natural-origin RS", ylab="size parameter Natural-origin",
                       zlab="Power", main=paste(stream, " RRS = 0.5, Prop F1 = ", prop, sep=""), highlight.3d=TRUE)
}


# Create matrix
newmat <- matrix(data = out[,c("Perm.Power")], nrow = length(unique(out[, "mu.n"])), ncol = length(unique(out[, "size.n"])), dimnames = list(unique(out[, "mu.n"]), unique(out[, "size.n"])), byrow = TRUE)


# Color contour plot ####
#par(opar)
par(mar=c(5.1,6.1,3.1,6.1))
par(mfrow = c(1, 1))

par(bg=colors()[c(356)])
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
zmin=0
zmax=1

filled.contour3(x = as.numeric(rownames(newmat)),y = as.numeric(colnames(newmat)), z = newmat, color = jet.colors,
                xlim= c (0, max(as.numeric(rownames(newmat)))), ylim = c(0, max(as.numeric(colnames(newmat)))), zlim = c(zmin,zmax),
                xlab = "Mu.N RS", ylab = "Dispersion", main = paste(stream, "F1 sample prop", prop, sep = " "), cex.lab = 2, axes = TRUE)
contour(x = as.numeric(rownames(newmat)), y = as.numeric(colnames(newmat)), z = newmat, levels = 0.8, add = TRUE)

par(new = "TRUE",plt = c(0.87,0.92,0.25,0.85),las = 1,cex.axis = 1.5)
filled.legend(as.numeric(rownames(newmat)),y=as.numeric(colnames(newmat)),newmat,color=jet.colors,xlab="",ylab="",xlim=c(min(xintercepts),max(xintercepts)),ylim=c(min(slopes),max(slopes)),zlim=c(zmin,zmax))
filled.legend(as.numeric(rownames(newmat)),y=as.numeric(colnames(newmat)),newmat,color=jet.colors,xlab="",ylab="",xlim=c(min(xintercepts),max(xintercepts)),ylim=c(min(slopes),max(slopes)),zlim=c(zmin,zmax))
mtext("Power",cex=2,padj=-0.3)



# Change format to x = mu.n and y = prop F1 ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine all output into big output
### Main Dir
# setwd("V:/WORK/Pink/AHRG/Parentage simulations/Mark Christie")
# setwd("F:/2015 AHRP April Seattle")
setwd("V:/Analysis/5_Coastwide/Multispecies/Alaska Hatchery Research Program/Parentage Simulations/Simulations based on Christie et al. 2015 Box 2")

source("H:/R Source Scripts/Functions.GCL_KS.R")

#### January 4, 2022 ####

library(scatterplot3d)
rm(list = ls())
source("functions/filled.contour3.R")
setwd("output/rrs_0.8/male")

#### Target Dir
# setwd("StreamSpecific")  # var.n = var.h
# setwd("StreamSpecific_ConstSize")  # size.n = size.h
# setwd("StreamSpecific_ConstSize_RRS_0.8")  # size.n = size.h

#stream <- c("Spring", "Stockdale", "Gilmour", "HoganBay", "Erb", "Paddy", "Fish", "Admiralty", "Prospect", "Sawmill")[1]
stream <- c("Erb", "HoganBay", "Paddy", "Spring", "Stockdale", "Gilmour", "Admiralty", "Fish", "Prospect", "Sawmill")[10]
year <- 2017
props <- c(0.05, 0.1, 0.167, 0.333, 0.5, 0.667, 0.833, 1)
kSize <- 1
# for(kSize in c(1, 2, 5, 10)){
  
OUT <-
  setNames(object = sapply(props, function(prop) {
    read.table(file = grep(
      pattern = paste(year, "_prop_", prop, "_", sep = ""),
      x = list.files(
        path = getwd(),
        pattern = stream,
        full.names = TRUE
      ),
      value = TRUE
    ))
  }, simplify = FALSE),
  nm = props)

## Isolate a specific "size.n"
# different size.n values
unique(OUT[[1]][, "size.n"])  # if NA shows up, there is a problem

# Check for NAs
lapply(OUT, function (prop) {prop[which(is.na(prop[, "size.n"])), ]})
#OUT[["0.05"]][4, 1:3] <- c(0.5, 0.25, 10)  # Fix NAs

newlst <- lapply(OUT, function(prop) {prop[which(prop[,"size.n"] == kSize), ]})

newmat <- sapply(newlst, function(prop) {prop[, "Perm.Power"]})
rownames(newmat) <- unique(newlst[[1]][,"mu.n"])

par(mar = c(5.1, 6.1, 5.1, 6.1))
par(mfrow = c(1, 1))
par(bg = colors()[c(356)])
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
zmin <- 0
zmax <- 1

filled.contour3(
  x = as.numeric(rownames(newmat)),
  y = as.numeric(colnames(newmat)),
  z = newmat,
  color = jet.colors,
  xlim = c(0, 2),  # c(0, 5), changing to max 2 for 2017 power plots
  ylim = c(0, max(as.numeric(colnames(
    newmat
  )))),
  zlim = c(zmin, zmax),
  xlab = "",
  ylab = "",
  cex.lab = 2,
  axes = FALSE
)

mtext(
  expression(paste("F"[1], " sampling proportion", sep = "")),
  side = 2,
  line = 3.5,
  cex = 2,
  las = 0
)

mtext(
  side = 1,
  "Mean RS Natural-origin",
  cex = 2,
  line = 3
)

mtext(
  side = 3,
  paste(stream, " ", year, sep = ""),
  cex = 2,
  line = 3
)

mtext(
  side = 3,
  paste("Dispersion = ", kSize, "; N = ", OUT[[1]][1, "N.w.parents"], "; H = ", OUT[[1]][1, "N.h.parents"], sep = ""),
  cex = 2,
  line = 1
)

plot.axes = {
  axis(1, at = seq(0, 2, by = 0.5), cex.axis = 1.5)  # seq(0, 5, by = 1)
  
  axis(
    2,
    at = c(0, as.numeric(colnames(newmat))),
    labels = c("0", "1/20", "1/10", "1/6", "1/3", "1/2", "2/3", "5/6", "1"),
    las = 1,
    cex.axis = 1.5
  )
}

#points(x = rep(x = as.numeric(rownames(newmat)), each = ncol(newmat)), y = rep(x = as.numeric(colnames(newmat)), times = nrow(newmat)), cex = 2, pch = 16)  # add points to show how interpolations are done
contour(
  x = as.numeric(rownames(newmat)),
  y = as.numeric(colnames(newmat)),
  z = newmat,
  levels = 0.8,
  add = TRUE,
  lwd = 5,
  labcex = 2.5,
  method = "edge"
)

par(
  new = "TRUE",
  plt = c(0.87, 0.92, 0.20, 0.80),
  las = 1,
  cex.axis = 1.5
)
filled.legend(
  as.numeric(rownames(newmat)),
  y = as.numeric(colnames(newmat)),
  newmat,
  color = jet.colors,
  xlab = "",
  ylab = "",
  xlim = c(min(xintercepts), max(xintercepts)),
  ylim = c(min(slopes), max(slopes)),
  zlim = c(zmin, zmax)
)
filled.legend(
  as.numeric(rownames(newmat)),
  y = as.numeric(colnames(newmat)),
  newmat,
  color = jet.colors,
  xlab = "",
  ylab = "",
  xlim = c(min(xintercepts), max(xintercepts)),
  ylim = c(min(slopes), max(slopes)),
  zlim = c(zmin, zmax)
)
mtext("Power",
      cex = 2,
      padj = -0.3,
      adj = 0.2)
# }  # kSize
setwd("../../..")















# Testing
prop <- props[3]
read.table(file = grep(pattern = paste(year, "_prop_", prop, "_", sep=""), x = list.files(path=getwd(), pattern=stream, full.names=TRUE), value = TRUE))


# Fit lowess to each sample prop
par(mar=c(5.1, 5.1, 4.1, 2.1))
par(mfrow = c(1, 1))
apply(newmat, 2, function(col) {plot(col ~ as.numeric(rownames(newmat)), cex = 2, pch = 16, type = "p", ylim = c(0, 1), xlab = "Mu.N", ylab = "Power", cex.lab = 2)})


# Smooth out data ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NOTE: I would need to do this for 4 size.n's for all 10 streams for 2 years (i.e. 80 times)
smoothmat=NULL
i_span <- c(0.3, 0.6, 0.7, 0.7, 0.8, 0.7, 0.9, 0.7) #invisible(ifelse(stray==0.5, assign("i_span",c(0.3,0.4,0.3,0.3,0.34,0.3,0.25,0.4)),ifelse(stray==0.15, assign("i_span",c(0.6,0.5,0.6,0.4,0.3,0.35,0.3,0.3)), ifelse(stray==0.1, assign("i_span",c(0.6,0.7,0.7,0.6,0.7,0.5,0.5,0.5)), "HOSER: NEED TO SPECIFY i_span for each STRAY!!!"))))
for(i in 1:ncol(newmat)){
  plot(newmat[, i]~as.numeric(rownames(newmat)), ylim = c(0, 1), pch = 16, xlab = "Mu.N", ylab = "p F1 sampled", cex = 2, main = paste("Size.N", colnames(newmat))[i])
  #fit = lm(newmat[, i]~poly(as.numeric(rownames(newmat)), 4, raw = TRUE)) # poly isn't a good fit
  fit = loess(newmat[, i]~as.numeric(rownames(newmat)), span = i_span[i])
  smoothmat = cbind(smoothmat, predict(fit, data.frame(x = as.numeric(rownames(newmat)))))
  lines(as.numeric(rownames(newmat)), predict(fit, data.frame(x = as.numeric(rownames(newmat)))), col = 4, lwd = 3)
}
dimnames(smoothmat) = dimnames(newmat)

# Refining i_span
i=1
plot(newmat[, i] ~ as.numeric(rownames(newmat)), cex = 2, pch = 16, type = "p", xlab = "Mu.N", ylab = "Power", cex.lab = 2)
fit=loess(newmat[,i]~as.numeric(rownames(newmat)), span = 0.3) #,span=i_span[i]
lines(as.numeric(rownames(newmat)),predict(fit,data.frame(x=as.numeric(rownames(newmat)))),col=4,lwd=3)

# Remove negative numbers
smoothmat[which(smoothmat < 0)] <- 0

# Plot smooth ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mar=c(5.1,6.1,3.1,6.1))
par(mfrow = c(1, 1))
par(bg=colors()[c(356)])
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
zmin=0
zmax=1

filled.contour3(x = as.numeric(rownames(smoothmat)),y = as.numeric(colnames(smoothmat)), z = smoothmat, color = jet.colors,
                xlim= c (0, max(as.numeric(rownames(smoothmat)))), ylim = c(0, max(as.numeric(colnames(smoothmat)))), zlim = c(zmin,zmax),
                xlab = "", ylab = "", cex.lab = 2, axes = FALSE)

mtext(expression(paste("F"[1]," sampling proportion",sep="")),side=2,line=3.5,cex=2,las=0)
mtext(side=1,"Mean RS Natural-origin",cex=2,line=3)
plot.axes = { axis(1, at=seq(0,8,by=2), cex.axis=1.5);
              axis(2, at=c(0,as.numeric(colnames(smoothmat))),labels=c("0","1/20","1/10","1/6","1/3","1/2","2/3","5/6","1"),las=1,cex.axis=1.5) }

#points(x = rep(x = as.numeric(rownames(smoothmat)), each = ncol(smoothmat)), y = rep(x = as.numeric(colnames(smoothmat)), times = nrow(smoothmat)), cex = 2, pch = 16, add = TRUE)
contour(x = as.numeric(rownames(smoothmat)), y = as.numeric(colnames(smoothmat)), z = smoothmat, levels = 0.8, add = TRUE)

par(new = "TRUE",plt = c(0.87,0.92,0.25,0.85),las = 1,cex.axis = 1.5)
filled.legend(as.numeric(rownames(smoothmat)),y=as.numeric(colnames(smoothmat)),smoothmat,color=jet.colors,xlab="",ylab="",xlim=c(min(xintercepts),max(xintercepts)),ylim=c(min(slopes),max(slopes)),zlim=c(zmin,zmax))
filled.legend(as.numeric(rownames(smoothmat)),y=as.numeric(colnames(smoothmat)),smoothmat,color=jet.colors,xlab="",ylab="",xlim=c(min(xintercepts),max(xintercepts)),ylim=c(min(slopes),max(slopes)),zlim=c(zmin,zmax))
mtext("Power",cex=2,padj=-0.3)



# Show histograms of different distributions ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NOTE: ON 4/21/15 Jim recommends holding size constant, not  variance (i.e. size.n = size.h, not var.n = var.h)

nbnom_variance <- function(mu, size){
  invisible(mu+((mu^2)/size))
}

size2 <- function(mu1, size1, mu2){
  var1 <- nbnom_variance(mu1, size1)
  invisible((mu2^2)/(var1-mu2))
}

rrs.values <- 0.5
mu.n.values <- c(0.25, seq(from = 0.5, to = 5, by = 0.5)) # R/female S (i.e. R/S / 2) # THESE SIMULATIONS ARE TOTALLY DEPENDENT ON THE MU AND SIZE OF THE NATURAL DISTRIBUTION OF R/S
size.n.values <- c(1, 2, 5, 10) #0.95

par(mar = c(2.1, 2.1, 0, 0))
par(oma = c(5.6, 5.6, 0, 0))
par(mfrow = c(length(size.n.values), length(mu.n.values)))
par(bg=colors()[c(356)])

for(size.n in size.n.values){
    
  for (mu.n in mu.n.values) {
    mu.h <- mu.n*rrs.values
    var.n <- nbnom_variance(mu.n, size.n)
    size.h <-  size2(mu.n, size.n, mu.h) 

    plot(table(rnbinom(n = 1000, size = size.n, mu = mu.n)), col="red", type = "h", xlab = "", ylab = "", main = "", bty = "n", lwd = 3, xlim = c(0, 25), ylim = c(0, 800), axes = FALSE)
    points(table(rnbinom(n = 1000, size = size.n, mu = mu.h) + 0.35), col="blue", type = "h", lwd = 3)
    legend(x = -2, y = ifelse(size.n == 1 & mu.n == 0.25, 800, 0), legend = c("Hatchery", "Natural"), fill = c("blue", "red"), cex = 2, bty = "n")
    mtext(text = ifelse(size.n == 10, mu.n, ""), side = 1, line = 3.25, cex = 2)
    mtext(text = ifelse(mu.n == 0.25, size.n, ""), side = 2, line = 2, cex = 2)
    axis(side = 1, cex.axis = 0.8)
    axis(side = 2, labels = seq(from = 0, to = 0.8, by = 0.1), at = seq(from = 0, to = 800, by = 100))
    
  }
}
mtext(text = "Mu of Natural-origin", side = 1, outer = TRUE, cex=3, line = 4)
mtext(text = "Dispersion of Natural-origin", side = 2, outer = TRUE, cex=3, line = 2.5)


# Show single histograms of different distributions for April Meeting ####~~~~~
require(coin)

# Input variables: TODO
kRRS <- 0.8
size.n <- 5
mu.n <- 4
mu.h <- mu.n * kRRS
k.Ymax <- 200
k.nH <- 200 # 200
k.nN <- 500 # 500
sample_prop_off <- 0.5

# Ploting formatting
par(bg=colors()[c(356)])
par(mfrow = c(1, 1))
par(mar=c(5.1,7.1,2.1,2.1))

# Determine distribution based on Inputs
NWoffspring=rnbinom(n=k.nN,size=size.n,mu=mu.n)   #take a sample of offspring based on this distribution
NHoffspring=rnbinom(n=k.nH,size=size.n,mu=mu.h)
Noffspring=sum(NWoffspring,NHoffspring) # how many offspring produced given mu.n and RRS
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
mu.n <- mean(NWassigned); mu.h <- mean(NHassigned); kRRS <- mu.h / mu.n; kRRS

# Make Plot
plot(table(NWassigned), col = 2, type = "h", xlab = "Family Size", ylab = "", main = "", bty = "n", lwd = 7, xlim = c(0, 25), ylim = c(0, k.Ymax), axes = FALSE, cex.lab = 2.5)
points(table(NHassigned + 0.35), col = 4, type = "h", lwd = 7)
mtext(side=2,"Number of Families",cex=2.5,line=5)
axis(side = 1, cex.axis = 2, lwd = 3)
axis(side = 2, at = seq(from = 0, to = k.Ymax, by = k.Ymax / 4), cex.axis = 2, lwd = 3)
segments(x0 = mu.h, x1 = mu.h, y0 = 0, y1 = k.Ymax, lwd = 10, lty = 3, col = "darkblue")
text(mu.n + 1, 0.6 * k.Ymax, paste("RS = ", round(mu.h, 2), "; n = ", k.nH, sep = ""), cex = 2.5, col = 4, pos = 4)
segments(x0 = mu.n, x1 = mu.n, y0 = 0, y1 = k.Ymax, lwd = 10, lty = 3, col = "darkred")
text(mu.n + 1, 0.5 * k.Ymax, paste("RS = ", round(mu.n, 2), "; n = ", k.nN, sep = ""), cex = 2.5, col = 2, pos = 4)
text(10, 0.4 * k.Ymax, paste("RRS = ", round(kRRS, 2), sep = ""), col = 1, cex = 3, pos = 4)
text(10, 0.3 * k.Ymax, paste("prop F1 = ", round(sample_prop_off, 2), sep = ""), col = 1, cex = 2.5, pos = 4)
legend("topright",legend=c("Hatchery","Natural"),col=c(4,2),lwd=10,bty="n",cex=3)

# Permutation test
mean(NHassigned)/mean(NWassigned)
mydata_samp=data.frame(nOff=c(NWassigned,NHassigned),Origin=c(rep("W",k.nN),rep("H",k.nH)))
k.Pval <- round(as.numeric(pvalue(oneway_test(nOff~Origin, data=mydata_samp, distribution=approximate(B=10000), alternative="less"))),4)
text(10, 0.2 * k.Ymax, paste("p = ", round(k.Pval, 3), sep = ""), col = 1, cex = 3, pos = 4)

  


##### Plots for Pink and Chum #################################################

setwd("V:/DOC/Power point presentations/Pink and Chum Workshop/2015 Richmond BC/Plots")

# How to compare negative binomail distributions
par(bg=colors()[c(356)])
par(mar=c(5.1,7.1,2.1,2.1))

hist(rnbinom(n=200,size=0.95,mu=2),breaks=20,col=2,xlim=c(0,15),right=FALSE)

nat=table(rnbinom(n=200,size=1,mu=2))
nat=nat/sum(nat)

hat=table(rnbinom(n=50,size=1,mu=1))
dimnames(hat)=list(as.character(as.numeric(unlist(dimnames(hat)))+0.3))
hat=hat/sum(hat)

nat;hat

plot(nat,type="h",col=2,lwd=10,bty="n",xlim=c(0,14),ylim=c(0,max(hat)),axes=FALSE,xlab="Family Size",ylab="",cex.lab=2.5)
points(hat,type="h",col=4,lwd=10)
axis(side=1,cex.axis=2,lwd=3,at=seq(0,14,by=2))
axis(side=2,,cex.axis=2,lwd=3,at=seq(0,round(x=max(hat),digits=1),by=0.1),labels=paste(seq(0,round(x=max(hat),digits=1),by=0.1)*100,"%",sep=""),las=2)
mtext(side=2,"Percent of Families",cex=2.5,line=5)
legend("topright",legend=c("Hatchery","Natural"),col=c(4,2),lwd=10,bty="n",cex=3)
segments(x0=1,x1=1,y0=0,y1=0.45,lwd=10,lty=3,col="darkblue");text(2,0.25,"RS = 1",cex=2.5,col=4,pos=4)
segments(x0=2,x1=2,y0=0,y1=0.45,lwd=10,lty=3,col="darkred");text(2,0.2,"RS = 2",cex=2.5,col=2,pos=4)
text(7,0.225,"RRS = 0.5",col=1,cex=3,pos=4)



## How does sampling only a proportion of offspring affect results?
require(coin)
n.h=6
n.w=52
sample_prop_off=0.5
NWoffspring=rnbinom(n=n.w,size=1,mu=2)   #take a sample of offspring based on this distribution
NHoffspring=rnbinom(n=n.h,size=1,mu=1)

Noffspring=sum(NWoffspring,NHoffspring) # how many offspring produced given mu.n and RRS

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


mean(NHoffspring)/mean(NWoffspring)
mean(NHassigned)/mean(NWassigned)

mydata_samp=data.frame(nOff=c(NWassigned,NHassigned),Origin=c(rep("W",n.w),rep("H",n.h)))
mydata_all=data.frame(nOff=c(NWoffspring,NHoffspring),Origin=c(rep("W",n.w),rep("H",n.h)))

round(as.numeric(pvalue(oneway_test(nOff~Origin, data=mydata_samp, distribution=approximate(B=10000), alternative="less"))),4)
round(as.numeric(pvalue(oneway_test(nOff~Origin, data=mydata_all, distribution=approximate(B=10000), alternative="less"))),4)

# Proportion sampled
nat=table(NWassigned)/sum(table(NWassigned))
hat=table(NHassigned)/sum(table(NHassigned))
dimnames(hat)=list(as.character(as.numeric(unlist(dimnames(hat)))+0.3))

plot(nat,type="h",col=2,lwd=10,bty="n",xlim=c(0,14),ylim=c(0, 0.8),axes=FALSE,xlab="Family Size",ylab="",cex.lab=2.5)
points(hat,type="h",col=4,lwd=10)
axis(side=1,cex.axis=2,lwd=3,at=seq(0,14,by=2))
axis(side=2,cex.axis=2,lwd=3,at=seq(0,round(x=0.8,digits=1),by=0.1),labels=paste(seq(0,round(x=0.8,digits=1),by=0.1)*100,"%",sep=""),las=2)
mtext(side=2,"Percent of Families",cex=2.5,line=5)
legend("topright",legend=c("Hatchery","Natural"),col=c(4,2),lwd=10,bty="n",cex=3)
RSh=mean(NHassigned)
RSw=mean(NWassigned)
segments(x0=RSh,x1=RSh,y0=0,y1=max(hat),lwd=10,lty=3,col="darkblue");text(2,0.4,paste("RS = ",round(RSh,2),sep=""),cex=2.5,col=4,pos=4)
segments(x0=RSw,x1=RSw,y0=0,y1=max(hat),lwd=10,lty=3,col="darkred");text(2,0.3,paste("RS = ",round(RSw,2),sep=""),cex=2.5,col=2,pos=4)
text(7.25,0.35,paste("RRS = ",round(RSh/RSw,2),sep=""),col=1,cex=3,pos=4)

# All sampled
nat=table(NWoffspring)/sum(table(NWoffspring))
hat=table(NHoffspring)/sum(table(NHoffspring))
dimnames(hat)=list(as.character(as.numeric(unlist(dimnames(hat)))+0.3))

plot(nat,type="h",col=2,lwd=10,bty="n",xlim=c(0,14),ylim=c(0, 0.8),axes=FALSE,xlab="Family Size",ylab="",cex.lab=2.5)
points(hat,type="h",col=4,lwd=10)
axis(side=1,cex.axis=2,lwd=3,at=seq(0,14,by=2))
axis(side=2,cex.axis=2,lwd=3,at=seq(0,round(x=0.8,digits=1),by=0.1),labels=paste(seq(0,round(x=0.8,digits=1),by=0.1)*100,"%",sep=""),las=2)
mtext(side=2,"Percent of Families",cex=2.5,line=5)
legend("topright",legend=c("Hatchery","Natural"),col=c(4,2),lwd=10,bty="n",cex=3)
RSh=mean(NHoffspring)
RSw=mean(NWoffspring)
segments(x0=RSh,x1=RSh,y0=0,y1=max(hat),lwd=10,lty=3,col="darkblue");text(2,0.4,paste("RS = ",round(RSh,2),sep=""),cex=2.5,col=4,pos=4)
segments(x0=RSw,x1=RSw,y0=0,y1=max(hat),lwd=10,lty=3,col="darkred");text(2,0.3,paste("RS = ",round(RSw,2),sep=""),cex=2.5,col=2,pos=4)
text(7.25,0.35,paste("RRS = ",round(RSh/RSw,2),sep=""),col=1,cex=3,pos=4)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## John H. Clark's idea 4/28/15 ####
# Instead of RRS compare the BY stray rate (proportion) to the proportion of H-origin fish in the return (of the "natural-origin" fish, what proportion had hatchery parents)

require(coin)

# Input variables: TODO
kRRS <- 0.5
size.n <- 5
mu.n <- 3
mu.h <- mu.n * kRRS
k.Ymax <- 200
k.nH <- 30 # 200
k.nN <- 200 # 500
sample_prop_off <- 0.25

# Ploting formatting
par(bg=colors()[c(356)])
par(mfrow = c(1, 1))
par(mar=c(5.1,7.1,2.1,2.1))

# Determine distribution based on Inputs
NWoffspring=rnbinom(n=k.nN,size=size.n,mu=mu.n)   #take a sample of offspring based on this distribution
NHoffspring=rnbinom(n=k.nH,size=size.n,mu=mu.h)
Noffspring=sum(NWoffspring,NHoffspring) # how many offspring produced given mu.n and RRS
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
mu.n <- mean(NWassigned); mu.h <- mean(NHassigned); kRRS <- mu.h / mu.n; kRRS

# Make Plot
plot(table(NWassigned), col = 2, type = "h", xlab = "Family Size", ylab = "", main = "", bty = "n", lwd = 7, xlim = c(0, 25), ylim = c(0, k.Ymax), axes = FALSE, cex.lab = 2.5)
points(table(NHassigned + 0.35), col = 4, type = "h", lwd = 7)
mtext(side=2,"Number of Families",cex=2.5,line=5)
axis(side = 1, cex.axis = 2, lwd = 3)
axis(side = 2, at = seq(from = 0, to = k.Ymax, by = k.Ymax / 4), cex.axis = 2, lwd = 3)
segments(x0 = mu.h, x1 = mu.h, y0 = 0, y1 = k.Ymax, lwd = 10, lty = 3, col = "darkblue")
text(mu.n + 1, 0.6 * k.Ymax, paste("RS = ", round(mu.h, 2), "; n = ", k.nH, sep = ""), cex = 2.5, col = 4, pos = 4)
segments(x0 = mu.n, x1 = mu.n, y0 = 0, y1 = k.Ymax, lwd = 10, lty = 3, col = "darkred")
text(mu.n + 1, 0.5 * k.Ymax, paste("RS = ", round(mu.n, 2), "; n = ", k.nN, sep = ""), cex = 2.5, col = 2, pos = 4)
text(10, 0.4 * k.Ymax, paste("RRS = ", round(kRRS, 2), sep = ""), col = 1, cex = 3, pos = 4)
text(10, 0.3 * k.Ymax, paste("prop F1 = ", round(sample_prop_off, 2), sep = ""), col = 1, cex = 2.5, pos = 4)
legend("topright",legend=c("Hatchery","Natural"),col=c(4,2),lwd=10,bty="n",cex=3)

# Permutation test
mean(NHassigned)/mean(NWassigned)
mydata_samp=data.frame(nOff=c(NWassigned,NHassigned),Origin=c(rep("W",k.nN),rep("H",k.nH)))
k.Pval <- round(as.numeric(pvalue(oneway_test(nOff~Origin, data=mydata_samp, distribution=approximate(B=10000), alternative="less"))),4)
text(10, 0.2 * k.Ymax, paste("p = ", round(k.Pval, 3), sep = ""), col = 1, cex = 3, pos = 4)

# Proportion test
k.prop.Pval <- prop.test(x = c(k.nH, n.NHsampled), n = c(sum(k.nH, k.nN), sum(n.NHsampled, n.NWsampled)), alternative = "greater")$p.value

# Proportions
matrix(data =  round(c(k.nH / sum(k.nH, k.nN), n.NHsampled / sum(n.NHsampled, n.NWsampled)), 3), ncol = 2, 
       dimnames = list("Proportion", c("Parents", "Offspring")))

# Fitness drop and p-value
matrix(data = round(c(kRRS, (n.NHsampled / sum(n.NHsampled, n.NWsampled)) / (k.nH / sum(k.nH, k.nN)), k.Pval, k.prop.Pval), 3),
       nrow = 2, ncol = 2, byrow = TRUE, dimnames = list(c("Fitness", "P-value"), c("Permutation", "Proportion")))


## Try the haploid genome ####


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
power_level=0.8
y=c(0.000,0.900,0.975,0.992,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000)
x=seq(0,3000,by=100)
fit=loess(y~x,span=0.2,weights=c(5,rep(1,length(x)-1)))
plot(y~x)
lines(predict(fit)~x, col='red', lwd=2)
newdat[which.min(abs(predict(fit,newdata=newdat)-power_level))]
abline(v=newdat[which.min(abs(predict(fit,newdata=newdat)-power_level))],lwd=3)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~







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


# 3D plot examples ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
threedpowerplot(mu=2,stray=0.50,rrs=0.80,type="perm")


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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

##### start of simulations to determine RRS point estimate 95% CI ####~~~~~~~~~



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

##### start of simulations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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















### Fixing span errors
ifelse(colnames(t(mat))[i]<0.15,0.6,0.3)

i=8
 
plot(t(mat)[,i]~rownames(t(mat)),xlim=c(0,3000),ylim=c(0,1),pch=16,xlab="F0s",ylab="RRS with power = 0.8",cex=2,main=paste("p F1 sampled",colnames(t(mat)))[i])
#fit=lm(t(mat)[,i]~poly(as.numeric(rownames(t(mat))),4,raw=TRUE)) # poly isn't a good fit
fit=loess(t(mat)[,i]~as.numeric(rownames(t(mat))),span=0.3)
newmat=cbind(newmat,predict(fit,data.frame(x=as.numeric(rownames(t(mat))))))
lines(as.numeric(rownames(t(mat))),predict(fit,data.frame(x=as.numeric(rownames(t(mat))))),col=4,lwd=3)


newmat=NULL
i_span=c(0.3,0.4,0.3,0.3,0.34,0.3,0.25,0.4)
for(i in 1:length(colnames(t(mat)))){
  plot(t(mat)[,i]~rownames(t(mat)),xlim=c(0,3000),ylim=c(0,1),pch=16,xlab="F0s",ylab="RRS with power = 0.8",cex=2,main=paste("p F1 sampled",colnames(t(mat)))[i])
  #fit=lm(t(mat)[,i]~poly(as.numeric(rownames(t(mat))),4,raw=TRUE)) # poly isn't a good fit
  fit=loess(t(mat)[,i]~as.numeric(rownames(t(mat))),span=i_span[i])
  newmat=cbind(newmat,predict(fit,data.frame(x=as.numeric(rownames(t(mat))))))
  lines(as.numeric(rownames(t(mat))),predict(fit,data.frame(x=as.numeric(rownames(t(mat))))),col=4,lwd=3)
}
dimnames(newmat)=dimnames(t(mat))
