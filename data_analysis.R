##################################MinT#################################################
##data 

dat<-subset(mint, Substrate!="Carbon Free" & Substrate!="Avicel")[,c(1:5,9:15, 18)]
summary(dat)

##libraries
library(deSolve)
library(dplyr)
library(FME)
library(reshape)
#######################################################################################
##############################First order decay########################################
#######################################################################################


#different substrates used in the experiment have slightly different concentrations at t=0
#but I would like to fit one function across substrates to estimate only one parameter
#I will try to create a function, which could do that
#the function use one model (first order decay) and apply it to 3 different substrates with
#3 different initial carbon concentrations. It combines the output into one data frame. 
#This data frame can be used with the modCost function, which can be subsequently used to
#estimate parameters across all data


#this is the function
first_order<-function(X, pars, t){
  
  #first order decay function
  deriv<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      dC<--k*C
      
      return(list(c(dC), r=k*C))
      
    })
  }
  
  #function to apply acros 3 different initial substrate concentrations
  base_function<-function(Ci){
    
    
    out<-ode(y=c(C=Ci), parms=pars, times=t, func=deriv)
    return(as.data.frame(out))}
  
  #results
  res<-lapply(X=X, FUN = "base_function") %>% bind_cols()
  
  return(res)

}


#this is the output
#X is the vector of initial carbon concentrations, pars and t are obvious
test<-first_order(X=c(100,200,300), pars = c(k=0.1), t=seq(0,120))




#now, I would like to see if I can estimate the parameter k across 3 different substrates
#with 3 different initial substrate concentrations

#I will create testing data data set
#since I plan to calibrate the model against respiration rate (r) I will generate
#data frame with two columns - time and r
ex<-first_order(X=c(100,200,300), pars = c(k=0.01), t=seq(0,120, by=10))[,c(1,3, 6, 9)]


#add some error
ex$r<-ex$r+rnorm(nrow(ex),mean=0, sd=0.01)
ex$r1<-ex$r1+rnorm(nrow(ex),mean=0, sd=0.02)
ex$r2<-ex$r2+rnorm(nrow(ex),mean=0, sd=0.03)


#create cost function
cost_function<-function(pars){
  
  out<-first_order(X=c(100,200,300), pars = pars, t=seq(0,120))
  cost<-modCost(model = out, obs = ex)
  
  return(cost)
  
}

#it calculates lot of things, which can be minimized 
cost_function(pars = c(k=0.01))


#it is minimized by finding the best parameter. This one is shown in summary output
mcmc<-modMCMC(f=cost_function, p=c(k=0.1), niter=500)
summary(mcmc)#the output shows k to be 0.0122, which is not bad


#let's hope it will work with real data
############################################################################################
Obs_decay<-dat[,c(2, 5, 12)]


Obs_decay<-reshape(Obs_decay, idvar = "Time", timevar = "Substrate", direction = "wide")

#Glucose = r, Cellobiose = r1, Mix = r2 
colnames(Obs_decay)<-c("time", "r", "r1", "r2")

#create cost function
cost_decay_all<-function(pars){
  
  out<-first_order(X=c(33.3*0.75+6.98,35.06*0.75+6.98, 22.6*0.75+6.98), pars = pars, t=seq(0,130))
  cost<-modCost(model = out, obs = Obs_decay)
  
  return(cost)
  
}

#estimate the parameter k across all substrates
mcmc_decay_all<-modMCMC(f=cost_decay_all, p=c(k=0.01), niter=10000)
summary(mcmc_decay_all)

#calculate the logLik and show the number of parameters
logLik_calc<-function(cost, pars){
  
      x<-as.data.frame(cost(pars = pars)$residuals)
      
      ll<-(-1/2/sd(x$obs)^2*sum((x$obs-x$mod)^2))
      npar=length(pars)
      
      return(c(logLik=ll, npar=npar))


}


#calculation
logLik_calc(cost=cost_decay_all, pars = c(k=summary(mcmc_decay_all)[6,1]))
