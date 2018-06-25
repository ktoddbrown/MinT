##################################MinT#################################################
##data 

dat<-subset(mint, Substrate!="Carbon Free" & Substrate!="Avicel")[,c(1:5,9:15, 18)]
summary(dat)

##libraries
library(deSolve)
library(dplyr)
library(FME)
library(reshape)
library(reshape2)
library(ggplot2)
library(foreach)
library(doParallel)
library(data.table)

##ggplot theme
theme_min<-theme(axis.text.x=element_text(vjust=0.2, size=14, colour="black"),
                 axis.text.y=element_text(hjust=0.2, size=14, colour="black"),
                 axis.title=element_text(size=14, colour="black"),
                 axis.line=element_line(size=0.5, colour="black"),
                 strip.text=element_text(size=16, face="bold"),
                 axis.ticks=element_line(size=1, colour="black"),
                 axis.ticks.length=unit(-0.05, "cm"),
                 panel.background=element_rect(colour="black", fill="white"),
                 panel.grid=element_line(linetype=0),
                 legend.text=element_text(size=14, colour="black"),
                 legend.title=element_text(size=14, colour="black"),
                 legend.position=c("right"),
                 legend.key.size=unit(1, "cm"),
                 strip.background=element_rect(fill="grey98", colour="black"),
                 legend.key=element_rect(fill="white", size=1.2),
                 legend.spacing=unit(0.5, "cm"),
                 plot.title=element_text(size=28, face="bold", hjust=-0.05))


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
mcmc<-modMCMC(f=cost_function, p=c(k=0.1), niter=5000)
summary(mcmc)#the output shows k to be 0.0104, which is not bad


#let's hope it will work with real data
############################################################################################
Obs_decay<-dat[,c(2, 5, 12)]

m1<-merge(Obs_decay[Obs_decay$Substrate=="Glucose", c(1,3)], 
      Obs_decay[Obs_decay$Substrate=="Cellobiose", c(1,3)],all = T, by="Time")

m2<-merge(m1, 
          Obs_decay[Obs_decay$Substrate=="Mix", c(1,3)],all = T, by="Time")

#Glucose = r, Cellobiose = r1, Mix = r2 
colnames(m2)<-c("time", "r", "r1", "r2")

#create cost function
cost_decay_all<-function(pars){
  
  out<-first_order(X=c(33.3*0.75+6.98, 35.06*0.75+6.98, 22.6*0.75+6.98), pars = pars, t=seq(0,130))
  cost<-modCost(model = out, obs = m2)
  
  return(cost)
  
}

#estimate the parameter k across all substrates
mcmc_decay_all<-modMCMC(f=cost_decay_all, p=c(k=0.01), niter=10000)
summary(mcmc_decay_all)

#calculate the logLik and show the number of parameters
logLik_calc<-function(cost, pars){
  
      x<-as.data.frame(cost(pars = pars)$residuals)
      mu<-mean(x$obs)
      variance<-sd(x$obs)^2
      
      
      ll<--1*sum((x$obs-x$mod)^2)/2/variance
      npar=length(pars)
      
      SSmodel<-sum((x$obs-x$mod)^2)
      SSdata<-sum((x$obs-mean(x$obs))^2)
      
      rsq=1-(SSmodel/SSdata)
      
      return(c(logLik=ll, npar=npar, rsq=rsq))


}


#calculation
logLik_calc(cost=cost_decay_all, pars = c(k=summary(mcmc_decay_all)[6,1]))

#Storing results
decay_all_ll<-logLik_calc(cost=cost_decay_all, pars = c(k=summary(mcmc_decay_all)[6,1]))

#Figure
plot_decay_all<-cost_decay_all(pars = c(k=summary(mcmc_decay_all)[6,1]))$residuals
ggplot(plot_decay_all, aes(obs, mod))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  ylab(expression(paste("Predicted respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Observed respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  ggtitle("First order decay function")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))
  

#########################################################################################
##First order decay for each substrate separately

decay_substrates<-function(data){
  
  #define 3 different data sets
  dat<-setDT(dat)[, id := .GRP, by = .(Substrate)]
  
  #first order decay function
  deriv<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      dC<--k*C
      
      return(list(c(dC), r=k*C))
      
    })
  }
  
  #initial carbon concentration
  cinit=c(33.3*0.75+6.98, 35.06*0.75+6.98, 22.6*0.75+6.98)
  
  
  #parameter estimation function
  estim<-function(data, cinit){
    
    Obs_dat<-data[, c(2,12)]
    colnames(Obs_dat)<-c("time", "r")
    
  
  
  #cost function
  cost_function<-function(pars){
    
    out<-as.data.frame(ode(y=c(C=cinit), parms = pars, t=seq(0,130), func = deriv))
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
  }
  
  res<-modMCMC(f=cost_function, p=c(k=0.01), niter=10000)
  
  return(res)
  
  }
  
  #parameter estimation
  res<-foreach(i=1:length(unique(dat$id)), .combine=list, .multicombine = TRUE,
               .packages=c("FME")) %dopar% {
                 
                 estim(data=dat[dat$id==i,], cinit = cinit[i])
                 
                 }
  
  #paraneter
  parameters<-data.frame(Substrate=c("Glucose", "Cellobiose", "Mix"),
                         k50=c(summary(res[[1]])[6,1],
                               summary(res[[2]])[6,1],
                               summary(res[[3]])[6,1]),
                         k25=c(summary(res[[1]])[5,1],
                               summary(res[[2]])[5,1],
                               summary(res[[3]])[5,1]),
                         k75=c(summary(res[[1]])[7,1],
                               summary(res[[2]])[7,1],
                               summary(res[[3]])[7,1]))
  

  #combining results
  cost_function2<-function(pars, data, cinit){
    
    Obs_dat2<-data[, c(2,12)]
    colnames(Obs_dat2)<-c("time", "r")
    
    
    out<-as.data.frame(ode(y=c(C=cinit), parms = pars, t=seq(0,130), func = deriv))
    cost<-modCost(model = out, obs = Obs_dat2)
    
    return(cost)
    
  }
  
  
  obs<-numeric()
  mod<-numeric()
  
  for(i in 1:3){
    
    obs<-append(obs, cost_function2(pars = c(k=parameters[i,2]), data=dat[dat$id==i,], cinit=cinit[i])$residuals$obs)
    mod<-append(mod, cost_function2(pars = c(k=parameters[i,2]), data=dat[dat$id==i,], cinit=cinit[i])$residuals$mod)
  }
  
  OvP<-data.frame(obs, mod, Substrate=c(rep("Glucose", 91),
                                        rep("Cellobiose", 96),
                                        rep("Mix", 95)))
  
  #logLik calculation
  mu<-mean(obs)
  variance<-sd(obs)^2
  
  
  ll<--1*sum((obs-mod)^2)/2/variance
  
  
  SSmodel<-sum((obs-mod)^2)
  SSdata<-sum((obs-mean(obs))^2)
  
  rsq=1-(SSmodel/SSdata)
  
  likelihood<-c(logLik=ll, npar=3, rsq=rsq)
  
  al<-list()
  al$parameters<-parameters
  al$OvP<-OvP
  al$likelihood<-likelihood
  
  
  
  return(al)
  
  }

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

decay_substrates_results<-decay_substrates(dat)

stopImplicitCluster()


#showing results
decay_substrates_results$likelihood
decay_substrates_results$parameters

#storing results
decay_substrates_ll<-decay_substrates_results$likelihood


#Figure
plot_decay_sbustrates<-decay_substrates_results$OvP
ggplot(plot_decay_sbustrates, aes(obs, mod))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  ylab(expression(paste("Predicted respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Observed respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  labs(title="First order decay function \n for different substrates")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

#Models comparison???
#I am not sure if it is correct
pchisq(-2*(as.numeric(decay_all_ll[1])-as.numeric(decay_substrates_ll[1])), df=2)



#########################################################################################
##First order decay for each structure separately

decay_structures<-function(data){
  
  #define 3 different data sets
  dat<-setDT(data)[, id := .GRP, by = .(Structure)]
  
  #function I need
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
  
  
  #parameter estimation function
  estim<-function(data){
    
  Obs_decay<-data[,c(2, 5, 12)]
  
  m1<-merge(Obs_decay[Obs_decay$Substrate=="Glucose", c(1,3)], 
            Obs_decay[Obs_decay$Substrate=="Cellobiose", c(1,3)],all = T, by="Time")
  
  m2<-merge(m1, 
            Obs_decay[Obs_decay$Substrate=="Mix", c(1,3)],all = T, by="Time")
  
  #Glucose = r, Cellobiose = r1, Mix = r2 
  colnames(m2)<-c("time", "r", "r1", "r2")
  
  #create cost function
  cost_decay_structure<-function(pars){
    
    out<-first_order(X=c(33.3*0.75+6.98, 35.06*0.75+6.98, 22.6*0.75+6.98), pars = pars, t=seq(0,130))
    cost<-modCost(model = out, obs = m2)
    
    return(cost)
    
  }
  
  mcmc<-modMCMC(f=cost_decay_structure, p=c(k=0.01), niter=10000)
  
  return(mcmc)
  
  }
  
  
  #parameter estimation
  res<-foreach(i=1:length(unique(dat$id)), .combine=list, .multicombine = TRUE,
               .packages=c("FME", "dplyr")) %dopar% {estim(data=dat[dat$id==i,])}
  
  #paraneter
  parameters<-data.frame(Substrate=c("Broth", "Glass wool", "Mixed glass"),
                         k50=c(summary(res[[1]])[6,1],
                               summary(res[[2]])[6,1],
                               summary(res[[3]])[6,1]),
                         k25=c(summary(res[[1]])[5,1],
                               summary(res[[2]])[5,1],
                               summary(res[[3]])[5,1]),
                         k75=c(summary(res[[1]])[7,1],
                               summary(res[[2]])[7,1],
                               summary(res[[3]])[7,1]))
  
  #combining results
  cost_decay_structure2<-function(pars, data){
    
    Obs_decay<-data[,c(2, 5, 12)]
    
    m1<-merge(Obs_decay[Obs_decay$Substrate=="Glucose", c(1,3)], 
              Obs_decay[Obs_decay$Substrate=="Cellobiose", c(1,3)],all = T, by="Time")
    
    m2<-merge(m1, 
              Obs_decay[Obs_decay$Substrate=="Mix", c(1,3)],all = T, by="Time")
    
    #Glucose = r, Cellobiose = r1, Mix = r2 
    colnames(m2)<-c("time", "r", "r1", "r2")
    
    out<-first_order(X=c(33.3*0.75+6.98, 35.06*0.75+6.98, 22.6*0.75+6.98), pars = pars, t=seq(0,130))
    cost<-modCost(model = out, obs = m2)
    
    return(cost)
    
  }
  
  obs<-numeric()
  mod<-numeric()
  
  for(i in 1:3){
    
    obs<-append(obs, cost_decay_structure2(pars = c(k=parameters[i,2]), data=dat[dat$id==i,])$residuals$obs)
    mod<-append(mod, cost_decay_structure2(pars = c(k=parameters[i,2]), data=dat[dat$id==i,])$residuals$mod)
  }
  
  OvP<-data.frame(obs, mod)
  
  #logLik calculation
  mu<-mean(obs)
  variance<-sd(obs)^2
  
  
  ll<--1*sum((obs-mod)^2)/2/variance
  
  
  SSmodel<-sum((obs-mod)^2)
  SSdata<-sum((obs-mean(obs))^2)
  
  rsq=1-(SSmodel/SSdata)
  
  likelihood<-c(logLik=ll, npar=3, rsq=rsq)
  
  al<-list()
  al$parameters<-parameters
  al$OvP<-OvP
  al$likelihood<-likelihood
  
  
  
  return(al)
  
  
}

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

decay_structures_results<-decay_structures(dat)

stopImplicitCluster()

#showing results
decay_structures_results$likelihood
decay_structures_results$parameters

#storing results
decay_structures_ll<-decay_structures_results$likelihood


#Figure
plot_decay_structures<-decay_structures_results$OvP
ggplot(plot_decay_structures, aes(obs, mod))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  ylab(expression(paste("Predicted respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Observed respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  labs(title="First order decay function \n for different structures")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

#Models comparison???
#I am not sure if it is correct
pchisq(-2*(as.numeric(decay_all_ll[1])-as.numeric(decay_substrates_ll[1])), df=2)

#########################################################################################
#First order decay for all combinations of Substrate and Structure 

decay_unique<-function(data){
  
  #define 3*3 different data sets
  dat<-setDT(data)[, id := .GRP, by = .(Substrate, Structure)]
  
  #first order decay function
  deriv<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      dC<--k*C
      
      return(list(c(dC), r=k*C))
      
    })
  }
  
  #parameter estimation function
  estim<-function(data){
    
    Obs_dat<-data[, c(2,12)]
    colnames(Obs_dat)<-c("time", "r")
    cinit<-as.numeric(data[1,"DOCinit"])*0.75+6.98
    
    
    #cost function
    cost_function<-function(pars){
      
      out<-as.data.frame(ode(y=c(C=cinit), parms = pars, t=seq(0,130), func = deriv))
      cost<-modCost(model = out, obs = Obs_dat)
      
      return(cost)
      
    }
  
    res<-modMCMC(f=cost_function, p=c(k=0.01), niter=10000)
    res$Substrate<-data[1,"Substrate"]
    res$Structure<-data[1,"Structure"]
    
    
    return(res)
    
  }
  
  #parameter estimation
  res<-foreach(i=1:length(unique(dat$id)), .combine=list, .multicombine = TRUE,
               .packages=c("FME")) %dopar% {
                 
                 estim(data=dat[dat$id==i,])
                 
               }
  
  
  #extract parameters
  k50<-numeric()
  k25<-numeric()
  k75<-numeric()
  Substrate<-character()
  Structure<-character()
  
  for(i in 1:9){
    
    k50<-append(k50, summary(res[[i]])[6,1])
    k25<-append(k25, summary(res[[i]])[5,1])
    k75<-append(k75, summary(res[[i]])[7,1])
    
    Substrate<-append(Substrate, res[[i]]$Substrate)
    Structure<-append(Structure, res[[i]]$Structure)
  }
  
  parameters<-data.frame(Substrate, Structure, k50, k25, k75)
  
  
  #OvP
  
  obs<-numeric()
  mod<-numeric()
  
  cost_function2<-function(pars, data){
    
    Obs_dat<-data[, c(2,12)]
    colnames(Obs_dat)<-c("time", "r")
    cinit<-as.numeric(data[1,"DOCinit"])*0.75+6.98
    
    out<-as.data.frame(ode(y=c(C=cinit), parms = pars, t=seq(0,130), func = deriv))
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
  }
  
  for(i in 1:9){
    
    obs<-append(obs, cost_function2(pars = c(k=as.numeric(parameters[i,3])), data=dat[dat$id==i,])$residuals$obs)
    mod<-append(mod, cost_function2(pars = c(k=as.numeric(parameters[i,3])), data=dat[dat$id==i,])$residuals$mod)
  }
  
  
  OvP<-data.frame(obs, mod)
  
  #logLik calculation
  mu<-mean(obs)
  variance<-sd(obs)^2
  
  
  ll<--1*sum((obs-mod)^2)/2/variance
  
  
  SSmodel<-sum((obs-mod)^2)
  SSdata<-sum((obs-mean(obs))^2)
  
  rsq=1-(SSmodel/SSdata)
  
  likelihood<-c(logLik=ll, npar=9, rsq=rsq)
  
  al<-list()
  al$parameters<-parameters
  al$OvP<-OvP
  al$likelihood<-likelihood
  
  
  
  return(al)
  
}

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

decay_unique_results<-decay_unique(dat)

stopImplicitCluster()

#showing results
decay_unique_results$likelihood
decay_unique_results$parameters

#storing results
decay_unique_ll<-decay_unique_results$likelihood


#Figure
plot_decay_unique<-decay_unique_results$OvP
ggplot(plot_decay_unique, aes(obs, mod))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  ylab(expression(paste("Predicted respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Observed respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  labs(title="First order decay function \n for different combinations of substrates and structures")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

#Models comparison???
#I am not sure if it is correct
pchisq(-2*(as.numeric(decay_all_ll[1])-as.numeric(decay_substrates_ll[1])), df=2)

decay_all_ll
decay_substrates_ll
decay_structures_ll
decay_unique_ll

###############################################################################################
#####################################Monod growth##############################################
###############################################################################################

#across all substrates and structures
monod_all<-function(data){
  
  #this is the monod growth function fitted across different substrates 
  #with different initial carbon concentrations
  monod<-function(X, pars, t){
    
    #monod growth function
    deriv<-function(time, state, pars){
      
      with(as.list(c(state, pars)),{
        
        dCmic<--k*Cmic+CUE*Vmax*Cmic*C/(Km+C)
        dC<-k*Cmic-Vmax*Cmic*C/(Km+C)
        
        return(list(c(dCmic, dC), r=(1-CUE)*Vmax*Cmic*C/(Km+C)))
        
      })
    }
    
    #function to apply acros 3 different initial substrate concentrations
    base_function<-function(Ci){
      
      
      out<-ode(y=c(Cmic=6.98, C=Ci), parms=pars, times=t, func=deriv)
      return(as.data.frame(out))}
    
    #results
    res<-lapply(X=X, FUN = "base_function") %>% bind_cols()
    
    return(res)
    
  }
  
  #create cost function
  estim<-function(data){
    
    Obs_dat<-data[,c(2, 5, 8, 12)]
    
    m1<-merge(Obs_dat[Obs_dat$Substrate=="Glucose", c(1,3,4)], 
              Obs_dat[Obs_dat$Substrate=="Cellobiose", c(1,3,4)],all = T, by="Time")
    
    m2<-merge(m1, 
              Obs_dat[Obs_dat$Substrate=="Mix", c(1,3,4)],all = T, by="Time")
    
    colnames(m2)<-c("time", "Cmic", "r", "Cmic1", "r1", "Cmic2", "r2")
    
  
  cost_function<-function(pars){
    
    
    out<-monod(X=c(33.3*0.75, 35.06*0.75, 22.6*0.75), pars = pars, t=seq(0,130))
    cost<-modCost(model = out, obs = m2, weight = "mean")
    
    return(cost)
    
  }
  
  res<-modMCMC(f=cost_function, p=c(Vmax=0.1, Km=3, CUE=0.8, k=0.01),
               lower=c(Vmax=1e-4, Km=1e-4, CUE=1e-2, k=1e-5),
               upper=c(Vmax=1e4, Km=1e4, CUE=0.999, k=1e5),niter=10000)
  
  return(res)
  
  }
  
  res<-estim(data)
  
  parameters<-summary(res)[6,]
  
  
  
  
  #OvP for respiration
  obs_r<-numeric()
  mod_r<-numeric()
  
  cost_r<-function(pars, data){
    
    Obs_dat<-data[,c(2, 5, 12)]
    
    m1<-merge(Obs_dat[Obs_dat$Substrate=="Glucose", c(1,3)], 
              Obs_dat[Obs_dat$Substrate=="Cellobiose", c(1,3)],all = T, by="Time")
    
    m2<-merge(m1, 
              Obs_dat[Obs_dat$Substrate=="Mix", c(1,3)],all = T, by="Time")
    
    colnames(m2)<-c("time", "r", "r1", "r2")
    
    out<-monod(X=c(33.3*0.75, 35.06*0.75, 22.6*0.75), pars = pars, t=seq(0,130))
    cost<-modCost(model = out, obs = m2)
    
    return(cost)
    
  }
  
  obs_r<-append(obs_r, cost_r(pars=summary(res)[6,], data=dat)$residuals$obs)
  mod_r<-append(mod_r, cost_r(pars=summary(res)[6,], data=dat)$residuals$mod)
  
  OvP_r<-data.frame(obs_r, mod_r)
  
  #logLik calculation
  mu_r<-mean(obs_r)
  variance_r<-sd(obs_r)^2
  
  
  ll_r<--1*sum((obs_r-mod_r)^2)/2/variance_r
  
  
  SSmodel_r<-sum((obs_r-mod_r)^2)
  SSdata_r<-sum((obs_r-mean(obs_r))^2)
  
  rsq_r=1-(SSmodel_r/SSdata_r)
  
  likelihood_r<-c(logLik=ll_r, npar=4, rsq=rsq_r)
  
  #OvP for biomass
  obs_Cmic<-numeric()
  mod_Cmic<-numeric()
  
  cost_Cmic<-function(pars, data){
    
    Obs_dat<-data[,c(2, 5, 8)]
    
    m1<-merge(Obs_dat[Obs_dat$Substrate=="Glucose", c(1,3)], 
              Obs_dat[Obs_dat$Substrate=="Cellobiose", c(1,3)],all = T, by="Time")
    
    m2<-merge(m1, 
              Obs_dat[Obs_dat$Substrate=="Mix", c(1,3)],all = T, by="Time")
    
    colnames(m2)<-c("time", "Cmic", "Cmic1", "Cmic2")
    
    out<-monod(X=c(33.3*0.75, 35.06*0.7, 22.6*0.75), pars = pars, t=seq(0,130))
    cost<-modCost(model = out, obs = m2)
    
    return(cost)
    
  }
  
  obs_Cmic<-append(obs_Cmic, cost_Cmic(pars=summary(res)[6,], data=dat)$residuals$obs)
  mod_Cmic<-append(mod_Cmic, cost_Cmic(pars=summary(res)[6,], data=dat)$residuals$mod)
  
  OvP_Cmic<-data.frame(obs_Cmic, mod_Cmic)
  
  #logLik calculation
  mu_Cmic<-mean(obs_Cmic)
  variance_Cmic<-sd(obs_Cmic)^2
  
  
  ll_Cmic<--1*sum((obs_Cmic-mod_Cmic)^2)/2/variance_Cmic
  
  
  SSmodel_Cmic<-sum((obs_Cmic-mod_Cmic)^2)
  SSdata_Cmic<-sum((obs_Cmic-mean(obs_Cmic))^2)
  
  rsq_Cmic=1-(SSmodel_Cmic/SSdata_Cmic)
  
  likelihood_Cmic<-c(logLik=ll_Cmic, npar=4, rsq=rsq_Cmic)
  
  al<-list()
  
  al$parameters<-parameters
  al$OvP_r<-OvP_r
  al$ll_r<-likelihood_r
  al$OvP_Cmic<-OvP_Cmic
  al$ll_Cmic<-likelihood_Cmic
  
  return(al)
}


monod_all_results<-monod_all(dat)

#showing results
monod_all_results$ll_r
monod_all_results$ll_Cmic
monod_all_results$parameters

#storing results
monod_all_ll<-monod_all_results$ll_r


#Figures
plot_monod_all<-monod_all_results$OvP_r
ggplot(plot_monod_all, aes(obs_r, mod_r))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  ylab(expression(paste("Predicted respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Observed respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  labs(title="Monod growth function")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

monod_all_Cmic<-monod_all_results$OvP_Cmic
ggplot(monod_all_Cmic, aes(obs_Cmic, mod_Cmic))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,30))+
  scale_y_continuous(limits = c(0,30))+
  ylab(expression(paste("Predicted ", C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Monod growth function")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

#Models comparison???
#I am not sure if it is correct
pchisq(-2*(as.numeric(decay_all_ll[1])-as.numeric(decay_substrates_ll[1])), df=2)

decay_all_ll
decay_substrates_ll
decay_structures_ll
decay_unique_ll
monod_all_ll


##########################################################################################
#for each structure separately

monod_structures<-function(data){
  
  #define 3 different data sets
  dat<-setDT(data)[, id := .GRP, by = .(Structure)]
  
  #this is the monod growth function fitted across different substrates 
  #with different initial carbon concentrations
  monod<-function(X, pars, t){
    
    #monod growth function
    deriv<-function(time, state, pars){
      
      with(as.list(c(state, pars)),{
        
        dCmic<--k*Cmic+CUE*Vmax*Cmic*C/(Km+C)
        dC<-k*Cmic-Vmax*Cmic*C/(Km+C)
        
        return(list(c(dCmic, dC), r=(1-CUE)*Vmax*Cmic*C/(Km+C)))
        
      })
    }
    
    #function to apply acros 3 different initial substrate concentrations
    base_function<-function(Ci){
      
      
      out<-ode(y=c(Cmic=6.98, C=Ci), parms=pars, times=t, func=deriv)
      return(as.data.frame(out))}
    
    #results
    res<-lapply(X=X, FUN = "base_function") %>% bind_cols()
    
    return(res)
    
  }
  
  #create cost function
  estim<-function(data){
    
    Obs_dat<-data[,c(2, 5, 8, 12)]
    
    m1<-merge(Obs_dat[Obs_dat$Substrate=="Glucose", c(1,3,4)], 
              Obs_dat[Obs_dat$Substrate=="Cellobiose", c(1,3,4)],all = T, by="Time")
    
    m2<-merge(m1, 
              Obs_dat[Obs_dat$Substrate=="Mix", c(1,3,4)],all = T, by="Time")
    
    colnames(m2)<-c("time", "Cmic", "r", "Cmic1", "r1", "Cmic2", "r2")
    
    
    cost_function<-function(pars){
      
      
      out<-monod(X=c(33.3*0.75, 35.06*0.75, 22.6*0.75), pars = pars, t=seq(0,130))
      cost<-modCost(model = out, obs = m2, weight = "mean")
      
      return(cost)
      
    }
    
    res<-modMCMC(f=cost_function, p=c(Vmax=0.1, Km=3, CUE=0.8, k=0.01),
                 lower=c(Vmax=1e-4, Km=1e-4, CUE=1e-2, k=1e-5),
                 upper=c(Vmax=1e4, Km=1e4, CUE=0.999, k=1e5),niter=10000)
    
    res$Structure<-data[1, "Structure"]
    
    return(res)
    
  }
  
  #parameter estimation
  res<-foreach(i=1:length(unique(dat$id)), .combine=list, .multicombine = TRUE,
               .packages=c("FME", "dplyr")) %dopar% {
                 
                 estim(data=dat[dat$id==i,])
                 
               }
  
  
  
  
  parameters<-rbind(as.data.frame(summary(res[[1]])[6,]), 
                    as.data.frame(summary(res[[2]])[6,]),
                    as.data.frame(summary(res[[3]])[6,]))
  parameters$Structure<-c(res[[1]]$Structure,
                          res[[2]]$Structure,
                          res[[3]]$Structure)
  
  parameters.l<-rbind(as.data.frame(summary(res[[1]])[5,]), 
                      as.data.frame(summary(res[[2]])[5,]),
                      as.data.frame(summary(res[[3]])[5,]))
  colnames(parameters.l)<-c("Vmax.l", "Km.l", "CUE.l", "k.l")
  
  parameters.u<-rbind(as.data.frame(summary(res[[1]])[7,]), 
                      as.data.frame(summary(res[[2]])[7,]),
                      as.data.frame(summary(res[[3]])[7,]))
  colnames(parameters.u)<-c("Vmax.u", "Km.u", "CUE.u", "k.u")
  
  
  parameters.all<-cbind(parameters, parameters.l, parameters.u)
  
  
  #OvP for respiration
  obs_r<-numeric()
  mod_r<-numeric()
  
  cost_r<-function(pars, data){
    
    Obs_dat<-data[,c(2, 5, 12)]
    
    m1<-merge(Obs_dat[Obs_dat$Substrate=="Glucose", c(1,3)], 
              Obs_dat[Obs_dat$Substrate=="Cellobiose", c(1,3)],all = T, by="Time")
    
    m2<-merge(m1, 
              Obs_dat[Obs_dat$Substrate=="Mix", c(1,3)],all = T, by="Time")
    
    colnames(m2)<-c("time", "r", "r1", "r2")
    
    out<-monod(X=c(33.3*0.75, 35.06*0.75, 22.6*0.75), pars = pars, t=seq(0,130))
    cost<-modCost(model = out, obs = m2)
    
    return(cost)
    
  }
  
  for(i in 1:3){
  
  obs_r<-append(obs_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
  mod_r<-append(mod_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
  
  }
  
  OvP_r<-data.frame(obs_r, mod_r)
  
  #logLik calculation
  mu_r<-mean(obs_r)
  variance_r<-sd(obs_r)^2
  
  
  ll_r<--1*sum((obs_r-mod_r)^2)/2/variance_r
  
  
  SSmodel_r<-sum((obs_r-mod_r)^2)
  SSdata_r<-sum((obs_r-mean(obs_r))^2)
  
  rsq_r=1-(SSmodel_r/SSdata_r)
  
  likelihood_r<-c(logLik=ll_r, npar=12, rsq=rsq_r)
  
  #OvP for biomass
  obs_Cmic<-numeric()
  mod_Cmic<-numeric()
  
  cost_Cmic<-function(pars, data){
    
    Obs_dat<-data[,c(2, 5, 8)]
    
    m1<-merge(Obs_dat[Obs_dat$Substrate=="Glucose", c(1,3)], 
              Obs_dat[Obs_dat$Substrate=="Cellobiose", c(1,3)],all = T, by="Time")
    
    m2<-merge(m1, 
              Obs_dat[Obs_dat$Substrate=="Mix", c(1,3)],all = T, by="Time")
    
    colnames(m2)<-c("time", "Cmic", "Cmic1", "Cmic2")
    
    out<-monod(X=c(33.3*0.75, 35.06*0.75, 22.6*0.75), pars = pars, t=seq(0,130))
    cost<-modCost(model = out, obs = m2)
    
    return(cost)
    
  }
  
  for(i in 1:3){
  obs_Cmic<-append(obs_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
  mod_Cmic<-append(mod_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
  }
  
  OvP_Cmic<-data.frame(obs_Cmic, mod_Cmic)
  
  #logLik calculation
  mu_Cmic<-mean(obs_Cmic)
  variance_Cmic<-sd(obs_Cmic)^2
  
  
  ll_Cmic<--1*sum((obs_Cmic-mod_Cmic)^2)/2/variance_Cmic
  
  
  SSmodel_Cmic<-sum((obs_Cmic-mod_Cmic)^2)
  SSdata_Cmic<-sum((obs_Cmic-mean(obs_Cmic))^2)
  
  rsq_Cmic=1-(SSmodel_Cmic/SSdata_Cmic)
  
  likelihood_Cmic<-c(logLik=ll_Cmic, npar=12, rsq=rsq_Cmic)
  
  al<-list()
  
  al$parameters<-parameters.all
  al$OvP_r<-OvP_r
  al$ll_r<-likelihood_r
  al$OvP_Cmic<-OvP_Cmic
  al$ll_Cmic<-likelihood_Cmic
  
  return(al)
}


no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

monod_structures_results<-monod_structures(dat)

stopImplicitCluster()

#showing results
monod_structures_results$ll_r
monod_structures_results$ll_Cmic
monod_structures_results$parameters

#storing results
monod_structures_ll<-monod_structures_results$ll_r


#Figures
plot_monod_structures<-monod_structures_results$OvP_r
ggplot(plot_monod_structures, aes(obs_r, mod_r))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  ylab(expression(paste("Predicted respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Observed respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  labs(title="Monod growth function \n for different structures")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

monod_structures_Cmic<-monod_structures_results$OvP_Cmic
ggplot(monod_structures_Cmic, aes(obs_Cmic, mod_Cmic))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,30))+
  scale_y_continuous(limits = c(0,30))+
  ylab(expression(paste("Predicted ", C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Monod growth function /n for different structures")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

#Models comparison???
#I am not sure if it is correct
pchisq(-2*(as.numeric(decay_all_ll[1])-as.numeric(decay_substrates_ll[1])), df=2)

decay_all_ll
decay_substrates_ll
decay_structures_ll
decay_unique_ll
monod_all_ll
monod_structures_ll

###########################################################################################
#for each substrate separately

monod_substrates<-function(data){
  
  #define 3 different data sets
  dat<-setDT(data)[, id := .GRP, by = .(Substrate)]
  
    #monod growth function
    deriv<-function(time, state, pars){
      
      with(as.list(c(state, pars)),{
        
        dCmic<--k*Cmic+CUE*Vmax*Cmic*C/(Km+C)
        dC<-k*Cmic-Vmax*Cmic*C/(Km+C)
        
        return(list(c(dCmic, dC), r=(1-CUE)*Vmax*Cmic*C/(Km+C)))
        
      })
    }
    
    
  #create cost function
  estim<-function(data){
    
    Obs_dat<-data[,c(2, 8, 12)]
    colnames(Obs_dat)<-c("time", "Cmic", "r")
    
    cinit<-as.numeric(data[1, "DOCinit"])*0.75
    
    cost_function<-function(pars){
      
      
      out<-as.data.frame(ode(y=c(Cmic=6.98, C=cinit), parms=pars, times=seq(0,130), func=deriv))
      cost<-modCost(model = out, obs = Obs_dat, weight = "mean")
      
      return(cost)
      
    }
    
    res<-modMCMC(f=cost_function, p=c(Vmax=0.1, Km=3, CUE=0.8, k=0.01),
                 lower=c(Vmax=1e-4, Km=1e-4, CUE=1e-2, k=1e-5),
                 upper=c(Vmax=1e4, Km=1e4, CUE=0.999, k=1e5),niter=10000)
    
    res$Substrate<-data[1, "Substrate"]
    
    return(res)
    
  }
  
  #parameter estimation
  res<-foreach(i=1:length(unique(dat$id)), .combine=list, .multicombine = TRUE,
               .packages=c("FME", "dplyr")) %dopar% {
                 
                 estim(data=dat[dat$id==i,])
                 
               }
  
  
  
  
  parameters<-rbind(as.data.frame(summary(res[[1]])[6,]), 
                    as.data.frame(summary(res[[2]])[6,]),
                    as.data.frame(summary(res[[3]])[6,]))
  parameters$Substrate<-c(res[[1]]$Substrate,
                          res[[2]]$Substrate,
                          res[[3]]$Substrate)
  
  parameters.l<-rbind(as.data.frame(summary(res[[1]])[5,]), 
                      as.data.frame(summary(res[[2]])[5,]),
                      as.data.frame(summary(res[[3]])[5,]))
  colnames(parameters.l)<-c("Vmax.l", "Km.l", "CUE.l", "k.l")
  
  parameters.u<-rbind(as.data.frame(summary(res[[1]])[7,]), 
                      as.data.frame(summary(res[[2]])[7,]),
                      as.data.frame(summary(res[[3]])[7,]))
  colnames(parameters.u)<-c("Vmax.u", "Km.u", "CUE.u", "k.u")
  
  
  parameters.all<-cbind(parameters, parameters.l, parameters.u)
  
  
  #OvP for respiration
  obs_r<-numeric()
  mod_r<-numeric()
  
  cost_r<-function(pars, data){
    
    Obs_dat<-data[,c(2, 12)]
    
    colnames(Obs_dat)<-c("time", "r")
    
    cinit<-as.numeric(data[1,"DOCinit"])*0.75
    
    out<-out<-as.data.frame(ode(y=c(Cmic=6.98, C=cinit), parms=pars, times=seq(0,130), func=deriv))
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
  }
  
  for(i in 1:3){
    
    obs_r<-append(obs_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
    mod_r<-append(mod_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
    
  }
  
  OvP_r<-data.frame(obs_r, mod_r)
  
  #logLik calculation
  mu_r<-mean(obs_r)
  variance_r<-sd(obs_r)^2
  
  
  ll_r<--1*sum((obs_r-mod_r)^2)/2/variance_r
  
  
  SSmodel_r<-sum((obs_r-mod_r)^2)
  SSdata_r<-sum((obs_r-mean(obs_r))^2)
  
  rsq_r=1-(SSmodel_r/SSdata_r)
  
  likelihood_r<-c(logLik=ll_r, npar=12, rsq=rsq_r)
  
  #OvP for biomass
  obs_Cmic<-numeric()
  mod_Cmic<-numeric()
  
  cost_Cmic<-function(pars, data){
    
    Obs_dat<-data[,c(2, 8)]
    
    colnames(Obs_dat)<-c("time", "Cmic")
    cinit<-as.numeric(data[1, "DOCinit"])*0.75
    
    out<-as.data.frame(ode(y=c(Cmic=6.98, C=cinit), parms=pars, times=seq(0,130), func=deriv))
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
  }
  
  for(i in 1:3){
    obs_Cmic<-append(obs_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
    mod_Cmic<-append(mod_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
  }
  
  OvP_Cmic<-data.frame(obs_Cmic, mod_Cmic)
  
  #logLik calculation
  mu_Cmic<-mean(obs_Cmic)
  variance_Cmic<-sd(obs_Cmic)^2
  
  
  ll_Cmic<--1*sum((obs_Cmic-mod_Cmic)^2)/2/variance_Cmic
  
  
  SSmodel_Cmic<-sum((obs_Cmic-mod_Cmic)^2)
  SSdata_Cmic<-sum((obs_Cmic-mean(obs_Cmic))^2)
  
  rsq_Cmic=1-(SSmodel_Cmic/SSdata_Cmic)
  
  likelihood_Cmic<-c(logLik=ll_Cmic, npar=12, rsq=rsq_Cmic)
  
  al<-list()
  
  al$parameters<-parameters.all
  al$OvP_r<-OvP_r
  al$ll_r<-likelihood_r
  al$OvP_Cmic<-OvP_Cmic
  al$ll_Cmic<-likelihood_Cmic
  
  return(al)
}

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

monod_substrates_results<-monod_substrates(dat)

stopImplicitCluster()

#showing results
monod_substrates_results$ll_r
monod_substrates_results$ll_Cmic
monod_substrates_results$parameters

#storing results
monod_substrates_ll<-monod_substrates_results$ll_r


#Figures
plot_monod_substrates<-monod_substrates_results$OvP_r
ggplot(plot_monod_substrates, aes(obs_r, mod_r))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  ylab(expression(paste("Predicted respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Observed respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  labs(title="Monod growth function \n for different substrates")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

monod_substrates_Cmic<-monod_substrates_results$OvP_Cmic
ggplot(monod_substrates_Cmic, aes(obs_Cmic, mod_Cmic))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,30))+
  scale_y_continuous(limits = c(0,30))+
  ylab(expression(paste("Predicted ", C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Monod growth function \n for different substrates")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

#Models comparison???
#I am not sure if it is correct
pchisq(-2*(as.numeric(decay_all_ll[1])-as.numeric(decay_substrates_ll[1])), df=2)

decay_all_ll
decay_substrates_ll
decay_structures_ll
decay_unique_ll
monod_all_ll
monod_structures_ll
monod_substrates_ll


###########################################################################################
#for all combinations of substrate and structure

monod_unique<-function(data){
  
  #define 3 different data sets
  dat<-setDT(data)[, id := .GRP, by = .(Substrate, Structure)]
  
  #monod growth function
  deriv<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      dCmic<--k*Cmic+CUE*Vmax*Cmic*C/(Km+C)
      dC<-k*Cmic-Vmax*Cmic*C/(Km+C)
      
      return(list(c(dCmic, dC), r=(1-CUE)*Vmax*Cmic*C/(Km+C)))
      
    })
  }
  
  
  #create cost function
  estim<-function(data){
    
    Obs_dat<-data[,c(2, 8, 12)]
    colnames(Obs_dat)<-c("time", "Cmic", "r")
    
    cinit<-as.numeric(data[1, "DOCinit"])*0.75
    
    cost_function<-function(pars){
      
      
      out<-as.data.frame(ode(y=c(Cmic=6.98, C=cinit), parms=pars, times=seq(0,130), func=deriv))
      cost<-modCost(model = out, obs = Obs_dat, weight = "mean")
      
      return(cost)
      
    }
    
    res<-modMCMC(f=cost_function, p=c(Vmax=0.1, Km=3, CUE=0.8, k=0.01),
                 lower=c(Vmax=1e-4, Km=1e-4, CUE=1e-2, k=1e-5),
                 upper=c(Vmax=1e4, Km=1e4, CUE=0.999, k=1e5),niter=10000)
    
    res$Substrate<-data[1, "Substrate"]
    res$Structure<-data[1, "Structure"]
    
    return(res)
    
  }
  
  #parameter estimation
  res<-foreach(i=1:length(unique(dat$id)), .combine=list, .multicombine = TRUE,
               .packages=c("FME", "dplyr")) %dopar% {
                 
                 estim(data=dat[dat$id==i,])
                 
               }
  
  
  parameters<-rbind(as.data.frame(summary(res[[1]])[6,]), 
                    as.data.frame(summary(res[[2]])[6,]),
                    as.data.frame(summary(res[[3]])[6,]),
                    as.data.frame(summary(res[[4]])[6,]), 
                    as.data.frame(summary(res[[5]])[6,]),
                    as.data.frame(summary(res[[6]])[6,]),
                    as.data.frame(summary(res[[7]])[6,]), 
                    as.data.frame(summary(res[[8]])[6,]),
                    as.data.frame(summary(res[[9]])[6,]))
  
  parameters$Substrate<-c(res[[1]]$Substrate,
                          res[[2]]$Substrate,
                          res[[3]]$Substrate,
                          res[[4]]$Substrate,
                          res[[5]]$Substrate,
                          res[[6]]$Substrate,
                          res[[7]]$Substrate,
                          res[[8]]$Substrate,
                          res[[9]]$Substrate)
  
  parameters.l<-rbind(as.data.frame(summary(res[[1]])[5,]), 
                      as.data.frame(summary(res[[2]])[5,]),
                      as.data.frame(summary(res[[3]])[5,]),
                      as.data.frame(summary(res[[4]])[5,]), 
                      as.data.frame(summary(res[[5]])[5,]),
                      as.data.frame(summary(res[[6]])[5,]),
                      as.data.frame(summary(res[[7]])[5,]), 
                      as.data.frame(summary(res[[8]])[5,]),
                      as.data.frame(summary(res[[9]])[5,]))
  
  colnames(parameters.l)<-c("Vmax.l", "Km.l", "CUE.l", "k.l")
  
  parameters.u<-rbind(as.data.frame(summary(res[[1]])[7,]), 
                      as.data.frame(summary(res[[2]])[7,]),
                      as.data.frame(summary(res[[3]])[7,]),
                      as.data.frame(summary(res[[4]])[7,]), 
                      as.data.frame(summary(res[[5]])[7,]),
                      as.data.frame(summary(res[[6]])[7,]),
                      as.data.frame(summary(res[[7]])[7,]), 
                      as.data.frame(summary(res[[8]])[7,]),
                      as.data.frame(summary(res[[9]])[7,]))
  colnames(parameters.u)<-c("Vmax.u", "Km.u", "CUE.u", "k.u")
  
  
  parameters.all<-cbind(parameters, parameters.l, parameters.u)
  
  
  #OvP for respiration
  obs_r<-numeric()
  mod_r<-numeric()
  
  cost_r<-function(pars, data){
    
    Obs_dat<-data[,c(2, 12)]
    
    colnames(Obs_dat)<-c("time", "r")
    
    cinit<-as.numeric(data[1,"DOCinit"])*0.75
    
    out<-out<-as.data.frame(ode(y=c(Cmic=6.98, C=cinit), parms=pars, times=seq(0,130), func=deriv))
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
  }
  
  for(i in 1:9){
    
    obs_r<-append(obs_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
    mod_r<-append(mod_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
    
  }
  
  OvP_r<-data.frame(obs_r, mod_r)
  
  #logLik calculation
  mu_r<-mean(obs_r)
  variance_r<-sd(obs_r)^2
  
  
  ll_r<--1*sum((obs_r-mod_r)^2)/2/variance_r
  
  
  SSmodel_r<-sum((obs_r-mod_r)^2)
  SSdata_r<-sum((obs_r-mean(obs_r))^2)
  
  rsq_r=1-(SSmodel_r/SSdata_r)
  
  likelihood_r<-c(logLik=ll_r, npar=36, rsq=rsq_r)
  
  #OvP for biomass
  obs_Cmic<-numeric()
  mod_Cmic<-numeric()
  
  cost_Cmic<-function(pars, data){
    
    Obs_dat<-data[,c(2, 8)]
    
    colnames(Obs_dat)<-c("time", "Cmic")
    cinit<-as.numeric(data[1, "DOCinit"])*0.75
    
    out<-as.data.frame(ode(y=c(Cmic=6.98, C=cinit), parms=pars, times=seq(0,130), func=deriv))
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
  }
  
  for(i in 1:9){
    obs_Cmic<-append(obs_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
    mod_Cmic<-append(mod_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
  }
  
  OvP_Cmic<-data.frame(obs_Cmic, mod_Cmic)
  
  #logLik calculation
  mu_Cmic<-mean(obs_Cmic)
  variance_Cmic<-sd(obs_Cmic)^2
  
  
  ll_Cmic<--1*sum((obs_Cmic-mod_Cmic)^2)/2/variance_Cmic
  
  
  SSmodel_Cmic<-sum((obs_Cmic-mod_Cmic)^2)
  SSdata_Cmic<-sum((obs_Cmic-mean(obs_Cmic))^2)
  
  rsq_Cmic=1-(SSmodel_Cmic/SSdata_Cmic)
  
  likelihood_Cmic<-c(logLik=ll_Cmic, npar=36, rsq=rsq_Cmic)
  
  al<-list()
  
  al$parameters<-parameters.all
  al$OvP_r<-OvP_r
  al$ll_r<-likelihood_r
  al$OvP_Cmic<-OvP_Cmic
  al$ll_Cmic<-likelihood_Cmic
  
  return(al)
}

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

monod_unique_results<-monod_unique(dat)

stopImplicitCluster()

#showing results
monod_unique_results$ll_r
monod_unique_results$ll_Cmic
monod_unique_results$parameters

#storing results
monod_unique_ll<-monod_unique_results$ll_r


#Figures
plot_monod_unique<-monod_unique_results$OvP_r
ggplot(plot_monod_unique, aes(obs_r, mod_r))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  ylab(expression(paste("Predicted respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Observed respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  labs(title="Monod growth function \n for all combinations of substrates and structures")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

monod_unique_Cmic<-monod_unique_results$OvP_Cmic
ggplot(monod_unique_Cmic, aes(obs_Cmic, mod_Cmic))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,30))+
  scale_y_continuous(limits = c(0,30))+
  ylab(expression(paste("Predicted ", C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Monod growth function \n for all combinations of substrates and structures")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

#Models comparison???
#I am not sure if it is correct
pchisq(-2*(as.numeric(decay_all_ll[1])-as.numeric(decay_substrates_ll[1])), df=2)

decay_all_ll
decay_substrates_ll
decay_structures_ll
decay_unique_ll
monod_all_ll
monod_structures_ll
monod_substrates_ll
monod_unique_ll
