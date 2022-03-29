library(mvtnorm)
library(extraDistr)
library(gtools)
library(msm)
library(parallel);
cl <- makeCluster(getOption("cl.cores", 4))

source("uniformization.R")
source("semi_markov_MCMC.R")
clusterEvalQ(cl, library(expm))


data=read.csv("destavola.csv")
data=data[,c(1,3,2,4)]

twoway3.q <- rbind(c(-0.25, 0.25, 0), c(0.1, -0.3, 0.2), c(0, 0, 0))
RateM.init=msm( state ~ time, subject=id, data = data,qmatrix = twoway3.q,deathexact=FALSE)$Qmatrices$baseline

out=gibbs.msm.SM(data=data,NMCMC=1100,RateM.init=10*RateM.init,f=0.001,g=0.001,abs=3,death=0,SM=1)

outSM=out[-(1000:11000),c(2,3,4,6,10,11)]

mean=apply(outSM,FUN=mean,MAR=2)
st.dev=apply(outSM,FUN=sd,MAR=2)
q1=apply(outSM,FUN=quantile,MAR=2,prob=c(0.025))
q2=apply(outSM,FUN=quantile,MAR=2,prob=c(0.975))

ris=cbind(mean,st.dev,q1,q2)
print(ris)