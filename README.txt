The function gibbs.msm.SM in the file semi_markov_MCMC.R performs the Metropolis within Gibbs algorithm for discretely observed semi-Markov processes described in the paper "Bayesian inference for discretely observed continuous time multi-state models" by R. Barone and A. Tancredi

Usage:

 gibbs.msm.SM(data,NMCMC,RateM.init,f,g,abs,death,SM)

Arguments:

data: a data frame with the state sequence and the observation times for all the observed units. The data frame must report also an identification number for each unit/patient, and a column indicating the first observation for each unit. See the file destavola.csv

NMCMC: the number of MCMC iterations

RateM.init: The initial matrix for the rate parameters

f: the shape parameter for the gamma prior of the parameters gamma_i
 
g: the rate parameter for the gamma prior of the parameters gamma_i

abs: the absorbing state. Select 0 if you do not have an absorbing state

death: Indicator variable for an exactly observed death (1) or not exactly observed (0)

SM: SM=1 fits the semi-Markov model, SM=0 fits the Markov model


To reproduce part of the result described in Section 4.1 use 

source("breast_cancer.R")




