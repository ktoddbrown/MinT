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


###############################################################################################
#####################################Microbial - enzyme model##################################
###############################################################################################

#across all substrates and structures
mem_all<-function(data){
  
  #this is the mem model tobe fitted across different substrates 
  #with different initial carbon concentrations
  mem<-function(X, pars, t){
    
    #mem model
    deriv<-function(time, state, pars){
      
      with(as.list(c(state, pars)),{
        
        dCmic<--kmic*Cmic+CUE*Vmax*Cmic*C/(Km+C)-pe*Cmic
        dC<-kmic*Cmic-Vmax*Cmic*C/(Km+C)+ke*E
        dE<-pe*Cmic-ke*E
        
        return(list(c(dCmic, dC, dE), r=(1-CUE)*Vmax*Cmic*C/(Km+C)))
        
      })
    }
    
    #function to apply acros 3 different initial substrate concentrations
    base_function<-function(Ci){
      
      
      out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=t, func=deriv)
      return(as.data.frame(out))}
    
    #results
    res<-lapply(X=X, FUN = "base_function") %>% bind_cols()
    
    return(res)
    
  }
  
  #create cost function
  estim<-function(data){
    
    Obs_dat<-data[,c(2, 5, 8, 12, 9)]
    
    m1<-merge(Obs_dat[Obs_dat$Substrate=="Glucose", c(1,3,4, 5)], 
              Obs_dat[Obs_dat$Substrate=="Cellobiose", c(1,3,4, 5)],all = T, by="Time")
    
    m2<-merge(m1, 
              Obs_dat[Obs_dat$Substrate=="Mix", c(1,3,4, 5)],all = T, by="Time")
    
    colnames(m2)<-c("time", "Cmic", "r","E", "Cmic1", "r1", "E1", "Cmic2", "r2", "E2")
    
    
    cost_function<-function(pars){
      
      
      out<-mem(X=c(33.3*0.75, 35.06*0.75, 22.6*0.75), pars = pars, t=seq(0,130))
      cost<-modCost(model = out, obs = m2, weight = "mean")
      
      return(cost)
      
    }
    
    res<-modMCMC(f=cost_function, p=c(Vmax=0.1, Km=3, CUE=0.8, kmic=0.01, ke=0.01, pe=0.01),
                 lower=c(Vmax=1e-4, Km=1e-4, CUE=1e-2, kmic=1e-5, ke=1e-6, pe=1e-6),
                 upper=c(Vmax=1e4, Km=1e4, CUE=0.999, kmic=1e5, ke=1e6, pe=1e6),niter=10000)
    
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
    
    out<-mem(X=c(33.3*0.75, 35.06*0.75, 22.6*0.75), pars = pars, t=seq(0,130))
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
  
  likelihood_r<-c(logLik=ll_r, npar=6, rsq=rsq_r)
  
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
    
    out<-mem(X=c(33.3*0.75, 35.06*0.7, 22.6*0.75), pars = pars, t=seq(0,130))
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
  
  likelihood_Cmic<-c(logLik=ll_Cmic, npar=6, rsq=rsq_Cmic)
  
  #OvP for enzymes
  obs_E<-numeric()
  mod_E<-numeric()
  
  cost_E<-function(pars, data){
    
    Obs_dat<-data[,c(2, 5, 9)]
    
    m1<-merge(Obs_dat[Obs_dat$Substrate=="Glucose", c(1,3)], 
              Obs_dat[Obs_dat$Substrate=="Cellobiose", c(1,3)],all = T, by="Time")
    
    m2<-merge(m1, 
              Obs_dat[Obs_dat$Substrate=="Mix", c(1,3)],all = T, by="Time")
    
    colnames(m2)<-c("time", "E", "E1", "E2")
    
    out<-mem(X=c(33.3*0.75, 35.06*0.7, 22.6*0.75), pars = pars, t=seq(0,130))
    cost<-modCost(model = out, obs = m2)
    
    return(cost)
    
  }
  
  obs_E<-append(obs_E, cost_E(pars=summary(res)[6,], data=dat)$residuals$obs)
  mod_E<-append(mod_E, cost_E(pars=summary(res)[6,], data=dat)$residuals$mod)
  
  OvP_E<-data.frame(obs_E, mod_E)
  
  #logLik calculation
  mu_E<-mean(obs_E)
  variance_E<-sd(obs_E)^2
  
  
  ll_E<--1*sum((obs_E-mod_E)^2)/2/variance_E
  
  
  SSmodel_E<-sum((obs_E-mod_E)^2)
  SSdata_E<-sum((obs_E-mean(obs_E))^2)
  
  rsq_E=1-(SSmodel_E/SSdata_E)
  
  likelihood_E<-c(logLik=ll_E, npar=6, rsq=rsq_E)
  
  al<-list()
  
  al$parameters<-parameters
  al$OvP_r<-OvP_r
  al$ll_r<-likelihood_r
  al$OvP_Cmic<-OvP_Cmic
  al$ll_Cmic<-likelihood_Cmic
  al$OvP_E<-OvP_E
  al$ll_E<-likelihood_E
  
  return(al)
}


mem_all_results<-mem_all(dat)

#showing results
mem_all_results$ll_r
mem_all_results$ll_Cmic
mem_all_results$ll_E
mem_all_results$parameters

#storing results
mem_all_ll<-mem_all_results$ll_r


#Figures
plot_mem_all<-mem_all_results$OvP_r
ggplot(plot_mem_all, aes(obs_r, mod_r))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  ylab(expression(paste("Predicted respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Observed respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  labs(title="Microbial - Enzyme model")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

mem_all_Cmic<-mem_all_results$OvP_Cmic
ggplot(mem_all_Cmic, aes(obs_Cmic, mod_Cmic))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,30))+
  scale_y_continuous(limits = c(0,30))+
  ylab(expression(paste("Predicted ", C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Microbial - Enzyme model")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

mem_all_E<-mem_all_results$OvP_E
ggplot(mem_all_E, aes(obs_E, mod_E))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,8))+
  scale_y_continuous(limits = c(0,8))+
  ylab(expression(paste("Predicted ", E, " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  E, " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Microbial - Enzyme model")+
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

mem_all_ll

#############################################################################################
#for different structures
mem_structures<-function(data){
  
  dat<-setDT(data)[, id := .GRP, by = .(Structure)]
  
  #this is the mem model tobe fitted across different substrates 
  #with different initial carbon concentrations
  mem<-function(X, pars, t){
    
    #mem model
    deriv<-function(time, state, pars){
      
      with(as.list(c(state, pars)),{
        
        dCmic<--kmic*Cmic+CUE*Vmax*Cmic*C/(Km+C)-pe*Cmic
        dC<-kmic*Cmic-Vmax*Cmic*C/(Km+C)+ke*E
        dE<-pe*Cmic-ke*E
        
        return(list(c(dCmic, dC, dE), r=(1-CUE)*Vmax*Cmic*C/(Km+C)))
        
      })
    }
    
    #function to apply acros 3 different initial substrate concentrations
    base_function<-function(Ci){
      
      
      out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=t, func=deriv)
      return(as.data.frame(out))}
    
    #results
    res<-lapply(X=X, FUN = "base_function") %>% bind_cols()
    
    return(res)
    
  }
  
  #create cost function
  estim<-function(data){
    
    Obs_dat<-data[,c(2, 5, 8, 12, 9)]
    
    m1<-merge(Obs_dat[Obs_dat$Substrate=="Glucose", c(1,3,4, 5)], 
              Obs_dat[Obs_dat$Substrate=="Cellobiose", c(1,3,4, 5)],all = T, by="Time")
    
    m2<-merge(m1, 
              Obs_dat[Obs_dat$Substrate=="Mix", c(1,3,4, 5)],all = T, by="Time")
    
    colnames(m2)<-c("time", "Cmic", "r","E", "Cmic1", "r1", "E1", "Cmic2", "r2", "E2")
    
    
    cost_function<-function(pars){
      
      
      out<-mem(X=c(33.3*0.75, 35.06*0.75, 22.6*0.75), pars = pars, t=seq(0,130))
      cost<-modCost(model = out, obs = m2, weight = "mean")
      
      return(cost)
      
    }
    
    res<-modMCMC(f=cost_function, p=c(Vmax=0.1, Km=3, CUE=0.8, kmic=0.01, ke=0.01, pe=0.01),
                 lower=c(Vmax=1e-4, Km=1e-4, CUE=1e-2, kmic=1e-5, ke=1e-6, pe=1e-6),
                 upper=c(Vmax=1e4, Km=1e4, CUE=0.999, kmic=1e5, ke=1e6, pe=1e6),niter=10000)
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
  colnames(parameters.l)<-c("Vmax.l", "Km.l", "CUE.l", "kmic.l", "ke.l", "pe.l")
  
  parameters.u<-rbind(as.data.frame(summary(res[[1]])[7,]), 
                      as.data.frame(summary(res[[2]])[7,]),
                      as.data.frame(summary(res[[3]])[7,]))
  colnames(parameters.u)<-c("Vmax.u", "Km.u", "CUE.u", "kmic.u", "ke.u", "pe.u")
  
  
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
    
    out<-mem(X=c(33.3*0.75, 35.06*0.75, 22.6*0.75), pars = pars, t=seq(0,130))
    cost<-modCost(model = out, obs = m2)
    
    return(cost)
    
  }
  
  for(i in 1:3){
    
    obs_r<-append(obs_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$obs)
    mod_r<-append(mod_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$mod)
    
  }
  
  OvP_r<-data.frame(obs_r, mod_r)
  
  #logLik calculation
  mu_r<-mean(obs_r)
  variance_r<-sd(obs_r)^2
  
  
  ll_r<--1*sum((obs_r-mod_r)^2)/2/variance_r
  
  
  SSmodel_r<-sum((obs_r-mod_r)^2)
  SSdata_r<-sum((obs_r-mean(obs_r))^2)
  
  rsq_r=1-(SSmodel_r/SSdata_r)
  
  likelihood_r<-c(logLik=ll_r, npar=18, rsq=rsq_r)
  
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
    
    out<-mem(X=c(33.3*0.75, 35.06*0.7, 22.6*0.75), pars = pars, t=seq(0,130))
    cost<-modCost(model = out, obs = m2)
    
    return(cost)
    
  }
  
  for(i in 1:3){
    obs_Cmic<-append(obs_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$obs)
    mod_Cmic<-append(mod_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$mod)
  }
  
  OvP_Cmic<-data.frame(obs_Cmic, mod_Cmic)
  
  #logLik calculation
  mu_Cmic<-mean(obs_Cmic)
  variance_Cmic<-sd(obs_Cmic)^2
  
  
  ll_Cmic<--1*sum((obs_Cmic-mod_Cmic)^2)/2/variance_Cmic
  
  
  SSmodel_Cmic<-sum((obs_Cmic-mod_Cmic)^2)
  SSdata_Cmic<-sum((obs_Cmic-mean(obs_Cmic))^2)
  
  rsq_Cmic=1-(SSmodel_Cmic/SSdata_Cmic)
  
  likelihood_Cmic<-c(logLik=ll_Cmic, npar=18, rsq=rsq_Cmic)
  
  #OvP for enzymes
  obs_E<-numeric()
  mod_E<-numeric()
  
  cost_E<-function(pars, data){
    
    Obs_dat<-data[,c(2, 5, 9)]
    
    m1<-merge(Obs_dat[Obs_dat$Substrate=="Glucose", c(1,3)], 
              Obs_dat[Obs_dat$Substrate=="Cellobiose", c(1,3)],all = T, by="Time")
    
    m2<-merge(m1, 
              Obs_dat[Obs_dat$Substrate=="Mix", c(1,3)],all = T, by="Time")
    
    colnames(m2)<-c("time", "E", "E1", "E2")
    
    out<-mem(X=c(33.3*0.75, 35.06*0.7, 22.6*0.75), pars = pars, t=seq(0,130))
    cost<-modCost(model = out, obs = m2)
    
    return(cost)
    
  }
  
  for(i in 1:3){
    obs_E<-append(obs_E, cost_E(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
    mod_E<-append(mod_E, cost_E(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
  }
  
  OvP_E<-data.frame(obs_E, mod_E)
  
  #logLik calculation
  mu_E<-mean(obs_E)
  variance_E<-sd(obs_E)^2
  
  
  ll_E<--1*sum((obs_E-mod_E)^2)/2/variance_E
  
  
  SSmodel_E<-sum((obs_E-mod_E)^2)
  SSdata_E<-sum((obs_E-mean(obs_E))^2)
  
  rsq_E=1-(SSmodel_E/SSdata_E)
  
  likelihood_E<-c(logLik=ll_E, npar=18, rsq=rsq_E)
  
  al<-list()
  
  al$parameters<-parameters.all
  al$OvP_r<-OvP_r
  al$ll_r<-likelihood_r
  al$OvP_Cmic<-OvP_Cmic
  al$ll_Cmic<-likelihood_Cmic
  al$OvP_E<-OvP_E
  al$ll_E<-likelihood_E
  
  return(al)
}

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

mem_structures_results<-mem_structures(dat)

stopImplicitCluster()

#showing results
mem_structures_results$ll_r
mem_structures_results$ll_Cmic
mem_structures_results$ll_E
mem_structures_results$parameters

#storing results
mem_structures_ll<-mem_structures_results$ll_r


#Figures
plot_mem_structures<-mem_structures_results$OvP_r
ggplot(plot_mem_structures, aes(obs_r, mod_r))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  ylab(expression(paste("Predicted respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Observed respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  labs(title="Microbial - Enzyme model \n for different structures")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

mem_structures_Cmic<-mem_structures_results$OvP_Cmic
ggplot(mem_structures_Cmic, aes(obs_Cmic, mod_Cmic))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,30))+
  scale_y_continuous(limits = c(0,30))+
  ylab(expression(paste("Predicted ", C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Microbial - Enzyme model \n for different structures")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

mem_structures_E<-mem_structures_results$OvP_E
ggplot(mem_structures_E, aes(obs_E, mod_E))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,8))+
  scale_y_continuous(limits = c(0,8))+
  ylab(expression(paste("Predicted ", E, " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  E, " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Microbial - Enzyme model \n for different structures")+
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

mem_all_ll
mem_structures_ll

##############################################################################################
#for different substrates
mem_substrates<-function(data){
  
  dat<-setDT(data)[, id := .GRP, by = .(Substrate)]
  
  #this is the mem model tobe fitted across different substrates 
  #with different initial carbon concentrations
  
  
  #mem model
  deriv<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      dCmic<--kmic*Cmic+CUE*Vmax*Cmic*C/(Km+C)-pe*Cmic
      dC<-kmic*Cmic-Vmax*Cmic*C/(Km+C)+ke*E
      dE<-pe*Cmic-ke*E
      
      return(list(c(dCmic, dC, dE), r=(1-CUE)*Vmax*Cmic*C/(Km+C)))
      
    })
  }
  
  
  
  #create cost function
  estim<-function(data){
    
    Obs_dat<-data[,c(2, 8, 12, 9)]
    
    colnames(Obs_dat)<-c("time", "Cmic", "r","E")
    Ci<-as.numeric(data[1, "DOCinit"])*0.75
    
    cost_function<-function(pars){
      
      
      out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
      cost<-modCost(model = out, obs = Obs_dat, weight = "mean")
      
      return(cost)
      
    }
    
    res<-modMCMC(f=cost_function, p=c(Vmax=0.1, Km=3, CUE=0.8, kmic=0.01, ke=0.01, pe=0.01),
                 lower=c(Vmax=1e-4, Km=1e-4, CUE=1e-2, kmic=1e-5, ke=1e-6, pe=1e-6),
                 upper=c(Vmax=1e4, Km=1e4, CUE=0.999, kmic=1e5, ke=1e6, pe=1e6),niter=10000)
    res$Structure<-data[1, "Substrate"]
    
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
  colnames(parameters.l)<-c("Vmax.l", "Km.l", "CUE.l", "kmic.l", "ke.l", "pe.l")
  
  parameters.u<-rbind(as.data.frame(summary(res[[1]])[7,]), 
                      as.data.frame(summary(res[[2]])[7,]),
                      as.data.frame(summary(res[[3]])[7,]))
  colnames(parameters.u)<-c("Vmax.u", "Km.u", "CUE.u", "kmic.u", "ke.u", "pe.u")
  
  
  parameters.all<-cbind(parameters, parameters.l, parameters.u)
  
  
  
  
  #OvP for respiration
  obs_r<-numeric()
  mod_r<-numeric()
  
  cost_r<-function(pars, data){
    
    Obs_dat<-data[,c(2, 12)]
    
    colnames(Obs_dat)<-c("time", "r")
    Ci<-as.numeric(data[1, "DOCinit"])*0.75
    
    out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
  }
  
  for(i in 1:3){
    
    obs_r<-append(obs_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$obs)
    mod_r<-append(mod_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$mod)
    
  }
  
  OvP_r<-data.frame(obs_r, mod_r)
  
  #logLik calculation
  mu_r<-mean(obs_r)
  variance_r<-sd(obs_r)^2
  
  
  ll_r<--1*sum((obs_r-mod_r)^2)/2/variance_r
  
  
  SSmodel_r<-sum((obs_r-mod_r)^2)
  SSdata_r<-sum((obs_r-mean(obs_r))^2)
  
  rsq_r=1-(SSmodel_r/SSdata_r)
  
  likelihood_r<-c(logLik=ll_r, npar=18, rsq=rsq_r)
  
  #OvP for biomass
  obs_Cmic<-numeric()
  mod_Cmic<-numeric()
  
  cost_Cmic<-function(pars, data){
    
    Obs_dat<-data[,c(2, 8)]
    
    colnames(Obs_dat)<-c("time", "Cmic")
    Ci<-as.numeric(data[1, "DOCinit"])*0.75
    
    out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
    
    
  }
  
  for(i in 1:3){
    obs_Cmic<-append(obs_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$obs)
    mod_Cmic<-append(mod_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$mod)
  }
  
  OvP_Cmic<-data.frame(obs_Cmic, mod_Cmic)
  
  #logLik calculation
  mu_Cmic<-mean(obs_Cmic)
  variance_Cmic<-sd(obs_Cmic)^2
  
  
  ll_Cmic<--1*sum((obs_Cmic-mod_Cmic)^2)/2/variance_Cmic
  
  
  SSmodel_Cmic<-sum((obs_Cmic-mod_Cmic)^2)
  SSdata_Cmic<-sum((obs_Cmic-mean(obs_Cmic))^2)
  
  rsq_Cmic=1-(SSmodel_Cmic/SSdata_Cmic)
  
  likelihood_Cmic<-c(logLik=ll_Cmic, npar=18, rsq=rsq_Cmic)
  
  #OvP for enzymes
  obs_E<-numeric()
  mod_E<-numeric()
  
  cost_E<-function(pars, data){
    
    Obs_dat<-data[,c(2, 9)]
    
    colnames(Obs_dat)<-c("time", "E")
    Ci<-as.numeric(data[1, "DOCinit"])*0.75
    
    out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
  }
  
  for(i in 1:3){
    obs_E<-append(obs_E, cost_E(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
    mod_E<-append(mod_E, cost_E(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
  }
  
  OvP_E<-data.frame(obs_E, mod_E)
  
  #logLik calculation
  mu_E<-mean(obs_E)
  variance_E<-sd(obs_E)^2
  
  
  ll_E<--1*sum((obs_E-mod_E)^2)/2/variance_E
  
  
  SSmodel_E<-sum((obs_E-mod_E)^2)
  SSdata_E<-sum((obs_E-mean(obs_E))^2)
  
  rsq_E=1-(SSmodel_E/SSdata_E)
  
  likelihood_E<-c(logLik=ll_E, npar=18, rsq=rsq_E)
  
  al<-list()
  
  al$parameters<-parameters.all
  al$OvP_r<-OvP_r
  al$ll_r<-likelihood_r
  al$OvP_Cmic<-OvP_Cmic
  al$ll_Cmic<-likelihood_Cmic
  al$OvP_E<-OvP_E
  al$ll_E<-likelihood_E
  
  return(al)
}

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

mem_substrates_results<-mem_substrates(dat)

stopImplicitCluster()

#showing results
mem_substrates_results$ll_r
mem_substrates_results$ll_Cmic
mem_substrates_results$ll_E
mem_substrates_results$parameters

#storing results
mem_substrates_ll<-mem_substrates_results$ll_r


#Figures
plot_mem_substrates<-mem_substrates_results$OvP_r
ggplot(plot_mem_substrates, aes(obs_r, mod_r))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  ylab(expression(paste("Predicted respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Observed respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  labs(title="Microbial - Enzyme model \n for different substrates")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

mem_substrates_Cmic<-mem_substrates_results$OvP_Cmic
ggplot(mem_substrates_Cmic, aes(obs_Cmic, mod_Cmic))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,30))+
  scale_y_continuous(limits = c(0,30))+
  ylab(expression(paste("Predicted ", C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Microbial - Enzyme model \n for different substrates")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

mem_substrates_E<-mem_substrates_results$OvP_E
ggplot(mem_substrates_E, aes(obs_E, mod_E))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,8))+
  scale_y_continuous(limits = c(0,8))+
  ylab(expression(paste("Predicted ", E, " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  E, " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Microbial - Enzyme model \n for different substrates")+
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

mem_all_ll
mem_structures_ll
mem_substrates_ll

#############################################################################################
#for all combinations of substrates and structures
mem_unique<-function(data){
  
  dat<-setDT(data)[, id := .GRP, by = .(Substrate, Structure)]
  
  #this is the mem model tobe fitted across different substrates 
  #with different initial carbon concentrations
  
  
  #mem model
  deriv<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      dCmic<--kmic*Cmic+CUE*Vmax*Cmic*C/(Km+C)-pe*Cmic
      dC<-kmic*Cmic-Vmax*Cmic*C/(Km+C)+ke*E
      dE<-pe*Cmic-ke*E
      
      return(list(c(dCmic, dC, dE), r=(1-CUE)*Vmax*Cmic*C/(Km+C)))
      
    })
  }
  
  
  
  #create cost function
  estim<-function(data){
    
    Obs_dat<-data[,c(2, 8, 12, 9)]
    
    colnames(Obs_dat)<-c("time", "Cmic", "r","E")
    Ci<-as.numeric(data[1, "DOCinit"])*0.75
    
    cost_function<-function(pars){
      
      
      out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
      cost<-modCost(model = out, obs = Obs_dat, weight = "mean")
      
      return(cost)
      
    }
    
    res<-modMCMC(f=cost_function, p=c(Vmax=0.1, Km=3, CUE=0.8, kmic=0.01, ke=0.01, pe=0.01),
                 lower=c(Vmax=1e-4, Km=1e-4, CUE=1e-2, kmic=1e-5, ke=1e-6, pe=1e-6),
                 upper=c(Vmax=1e4, Km=1e4, CUE=0.999, kmic=1e5, ke=1e6, pe=1e6),niter=10000)
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
  parameters$Structure<-c(res[[1]]$Structure,
                          res[[2]]$Structure,
                          res[[3]]$Structure,
                          res[[4]]$Structure,
                          res[[5]]$Structure,
                          res[[6]]$Structure,
                          res[[7]]$Structure,
                          res[[8]]$Structure,
                          res[[9]]$Structure)
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
  colnames(parameters.l)<-c("Vmax.l", "Km.l", "CUE.l", "kmic.l", "ke.l", "pe.l")
  
  parameters.u<-rbind(as.data.frame(summary(res[[1]])[7,]), 
                      as.data.frame(summary(res[[2]])[7,]),
                      as.data.frame(summary(res[[3]])[7,]), 
                      as.data.frame(summary(res[[4]])[7,]), 
                      as.data.frame(summary(res[[5]])[7,]),
                      as.data.frame(summary(res[[6]])[7,]),
                      as.data.frame(summary(res[[7]])[7,]), 
                      as.data.frame(summary(res[[8]])[7,]),
                      as.data.frame(summary(res[[9]])[7,]))
  colnames(parameters.u)<-c("Vmax.u", "Km.u", "CUE.u", "kmic.u", "ke.u", "pe.u")
  
  
  parameters.all<-cbind(parameters, parameters.l, parameters.u)
  
  
  
  
  #OvP for respiration
  obs_r<-numeric()
  mod_r<-numeric()
  
  cost_r<-function(pars, data){
    
    Obs_dat<-data[,c(2, 12)]
    
    colnames(Obs_dat)<-c("time", "r")
    Ci<-as.numeric(data[1, "DOCinit"])*0.75
    
    out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
  }
  
  for(i in 1:9){
    
    obs_r<-append(obs_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$obs)
    mod_r<-append(mod_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$mod)
    
  }
  
  OvP_r<-data.frame(obs_r, mod_r)
  
  #logLik calculation
  mu_r<-mean(obs_r)
  variance_r<-sd(obs_r)^2
  
  
  ll_r<--1*sum((obs_r-mod_r)^2)/2/variance_r
  
  
  SSmodel_r<-sum((obs_r-mod_r)^2)
  SSdata_r<-sum((obs_r-mean(obs_r))^2)
  
  rsq_r=1-(SSmodel_r/SSdata_r)
  
  likelihood_r<-c(logLik=ll_r, npar=54, rsq=rsq_r)
  
  #OvP for biomass
  obs_Cmic<-numeric()
  mod_Cmic<-numeric()
  
  cost_Cmic<-function(pars, data){
    
    Obs_dat<-data[,c(2, 8)]
    
    colnames(Obs_dat)<-c("time", "Cmic")
    Ci<-as.numeric(data[1, "DOCinit"])*0.75
    
    out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
    
    
  }
  
  for(i in 1:9){
    obs_Cmic<-append(obs_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$obs)
    mod_Cmic<-append(mod_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$mod)
  }
  
  OvP_Cmic<-data.frame(obs_Cmic, mod_Cmic)
  
  #logLik calculation
  mu_Cmic<-mean(obs_Cmic)
  variance_Cmic<-sd(obs_Cmic)^2
  
  
  ll_Cmic<--1*sum((obs_Cmic-mod_Cmic)^2)/2/variance_Cmic
  
  
  SSmodel_Cmic<-sum((obs_Cmic-mod_Cmic)^2)
  SSdata_Cmic<-sum((obs_Cmic-mean(obs_Cmic))^2)
  
  rsq_Cmic=1-(SSmodel_Cmic/SSdata_Cmic)
  
  likelihood_Cmic<-c(logLik=ll_Cmic, npar=54, rsq=rsq_Cmic)
  
  #OvP for enzymes
  obs_E<-numeric()
  mod_E<-numeric()
  
  cost_E<-function(pars, data){
    
    Obs_dat<-data[,c(2, 9)]
    
    colnames(Obs_dat)<-c("time", "E")
    Ci<-as.numeric(data[1, "DOCinit"])*0.75
    
    out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
  }
  
  for(i in 1:9){
    obs_E<-append(obs_E, cost_E(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
    mod_E<-append(mod_E, cost_E(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
  }
  
  OvP_E<-data.frame(obs_E, mod_E)
  
  #logLik calculation
  mu_E<-mean(obs_E)
  variance_E<-sd(obs_E)^2
  
  
  ll_E<--1*sum((obs_E-mod_E)^2)/2/variance_E
  
  
  SSmodel_E<-sum((obs_E-mod_E)^2)
  SSdata_E<-sum((obs_E-mean(obs_E))^2)
  
  rsq_E=1-(SSmodel_E/SSdata_E)
  
  likelihood_E<-c(logLik=ll_E, npar=54, rsq=rsq_E)
  
  al<-list()
  
  al$parameters<-parameters.all
  al$OvP_r<-OvP_r
  al$ll_r<-likelihood_r
  al$OvP_Cmic<-OvP_Cmic
  al$ll_Cmic<-likelihood_Cmic
  al$OvP_E<-OvP_E
  al$ll_E<-likelihood_E
  
  return(al)
}

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

mem_unique_results<-mem_unique(dat)

stopImplicitCluster()

#showing results
mem_unique_results$ll_r
mem_unique_results$ll_Cmic
mem_unique_results$ll_E
mem_unique_results$parameters

#storing results
mem_unique_ll<-mem_unique_results$ll_r


#Figures
plot_mem_unique<-mem_unique_results$OvP_r
ggplot(plot_mem_unique, aes(obs_r, mod_r))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  ylab(expression(paste("Predicted respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Observed respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  labs(title="Microbial - Enzyme model \n for different structures and substrates")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

mem_unique_Cmic<-mem_unique_results$OvP_Cmic
ggplot(mem_unique_Cmic, aes(obs_Cmic, mod_Cmic))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,30))+
  scale_y_continuous(limits = c(0,30))+
  ylab(expression(paste("Predicted ", C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Microbial - Enzyme model \n for different structures and substrates")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

mem_unique_E<-mem_unique_results$OvP_E
ggplot(mem_unique_E, aes(obs_E, mod_E))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,8))+
  scale_y_continuous(limits = c(0,8))+
  ylab(expression(paste("Predicted ", E, " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  E, " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Microbial - Enzyme model \n for different structures and substrates")+
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

mem_all_ll
mem_structures_ll
mem_substrates_ll
mem_unique_ll

###############################################################################################
###########################Metabolic Microbial - enzyme model##################################
###############################################################################################

#across all substrates and structures
mmem_all<-function(data){
  
  #this is the mem model tobe fitted across different substrates 
  #with different initial carbon concentrations
  mmem<-function(X, pars, t){
    
    #mmem model
    deriv<-function(time, state, pars){
      
      with(as.list(c(state, pars)),{
        
        #Carbon uptake
        
        Cu=Vmax*C*Cmic/(Km+C)
        
        #defining stochiometric coefficients
        #PO ratio
        #PO=1.5*(1-((Cu/Cmic)/((Cu/Cmic)+Sover)))
        y=2*1.5*Yatp
        
        #x
        #x=(Cu/Cmic)/((Cu/Cmic)+Sprod)
        
        #coefficients
        psi.r=1/((y*x)+(y*0.64*(1-x))+1)
        psi.g=y*x/((y*x)+(y*0.64*(1-x))+1)
        psi.e=0.64*y*(1-x)/((y*x)+(y*0.64*(1-x))+1)
        
        
        dCmic<--kmic*Cmic+Cu*psi.g
        dC<-kmic*Cmic-Cu+ke*E
        dE<-Cu*psi.e-ke*E
        
        return(list(c(dCmic, dC, dE), r=Cu*psi.r))
        
      })
    }
    
    #function to apply acros 3 different initial substrate concentrations
    base_function<-function(Ci){
      
      
      out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=t, func=deriv)
      return(as.data.frame(out))}
    
    #results
    res<-lapply(X=X, FUN = "base_function") %>% bind_cols()
    
    return(res)
    
  }
  
  #create cost function
  estim<-function(data){
    
    Obs_dat<-data[,c(2, 5, 8, 12, 9)]
    
    m1<-merge(Obs_dat[Obs_dat$Substrate=="Glucose", c(1,3,4, 5)], 
              Obs_dat[Obs_dat$Substrate=="Cellobiose", c(1,3,4, 5)],all = T, by="Time")
    
    m2<-merge(m1, 
              Obs_dat[Obs_dat$Substrate=="Mix", c(1,3,4, 5)],all = T, by="Time")
    
    colnames(m2)<-c("time", "Cmic", "r","E", "Cmic1", "r1", "E1", "Cmic2", "r2", "E2")
    
    
    cost_function<-function(pars){
      
      
      out<-mmem(X=c(33.3*0.75, 35.06*0.75, 22.6*0.75), pars = pars, t=seq(0,130))
      cost<-modCost(model = out, obs = m2, weight = "mean")
      
      return(cost)
      
    }
    
    res<-modMCMC(f=cost_function, p=c(Vmax=0.1, Km=3, kmic=0.01, ke=0.01, x=0.7, Yatp=1),
                 lower=c(Vmax=1e-4, Km=1e-4, kmic=1e-5, ke=1e-6, x=0, Yatp=0.01),
                 upper=c(Vmax=1e4, Km=1e4, kmic=1e5, ke=1e6, x=1, Yatp=30),niter=10000)
    
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
    
    out<-mmem(X=c(33.3*0.75, 35.06*0.75, 22.6*0.75), pars = pars, t=seq(0,130))
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
  
  likelihood_r<-c(logLik=ll_r, npar=6, rsq=rsq_r)
  
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
    
    out<-mmem(X=c(33.3*0.75, 35.06*0.7, 22.6*0.75), pars = pars, t=seq(0,130))
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
  
  likelihood_Cmic<-c(logLik=ll_Cmic, npar=6, rsq=rsq_Cmic)
  
  #OvP for enzymes
  obs_E<-numeric()
  mod_E<-numeric()
  
  cost_E<-function(pars, data){
    
    Obs_dat<-data[,c(2, 5, 9)]
    
    m1<-merge(Obs_dat[Obs_dat$Substrate=="Glucose", c(1,3)], 
              Obs_dat[Obs_dat$Substrate=="Cellobiose", c(1,3)],all = T, by="Time")
    
    m2<-merge(m1, 
              Obs_dat[Obs_dat$Substrate=="Mix", c(1,3)],all = T, by="Time")
    
    colnames(m2)<-c("time", "E", "E1", "E2")
    
    out<-mmem(X=c(33.3*0.75, 35.06*0.7, 22.6*0.75), pars = pars, t=seq(0,130))
    cost<-modCost(model = out, obs = m2)
    
    return(cost)
    
  }
  
  obs_E<-append(obs_E, cost_E(pars=summary(res)[6,], data=dat)$residuals$obs)
  mod_E<-append(mod_E, cost_E(pars=summary(res)[6,], data=dat)$residuals$mod)
  
  OvP_E<-data.frame(obs_E, mod_E)
  
  #logLik calculation
  mu_E<-mean(obs_E)
  variance_E<-sd(obs_E)^2
  
  
  ll_E<--1*sum((obs_E-mod_E)^2)/2/variance_E
  
  
  SSmodel_E<-sum((obs_E-mod_E)^2)
  SSdata_E<-sum((obs_E-mean(obs_E))^2)
  
  rsq_E=1-(SSmodel_E/SSdata_E)
  
  likelihood_E<-c(logLik=ll_E, npar=6, rsq=rsq_E)
  
  al<-list()
  
  al$parameters<-parameters
  al$OvP_r<-OvP_r
  al$ll_r<-likelihood_r
  al$OvP_Cmic<-OvP_Cmic
  al$ll_Cmic<-likelihood_Cmic
  al$OvP_E<-OvP_E
  al$ll_E<-likelihood_E
  
  return(al)
}


mmem_all_results<-mmem_all(dat)

#showing results
mmem_all_results$ll_r
mmem_all_results$ll_Cmic
mmem_all_results$ll_E
mmem_all_results$parameters

#storing results
mmem_all_ll<-mmem_all_results$ll_r


#Figures
plot_mmem_all<-mmem_all_results$OvP_r
ggplot(plot_mmem_all, aes(obs_r, mod_r))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  ylab(expression(paste("Predicted respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Observed respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  labs(title="Metabolic Microbial - Enzyme model")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

mmem_all_Cmic<-mmem_all_results$OvP_Cmic
ggplot(mmem_all_Cmic, aes(obs_Cmic, mod_Cmic))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,30))+
  scale_y_continuous(limits = c(0,30))+
  ylab(expression(paste("Predicted ", C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Microbial - Enzyme model \n for different structures and substrates")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

mmem_all_E<-mmem_all_results$OvP_E
ggplot(mmem_all_E, aes(obs_E, mod_E))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,8))+
  scale_y_continuous(limits = c(0,8))+
  ylab(expression(paste("Predicted ", E, " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  E, " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Metabolic Microbial - Enzyme model")+
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

mem_all_ll
mem_structures_ll
mem_substrates_ll
mem_unique_ll

mmem_all_ll

#############################################################################################
#for different structures
mmem_structures<-function(data){
  
  dat<-setDT(data)[, id := .GRP, by = .(Structure)]
  
  #this is the mmem model tobe fitted across different substrates 
  #with different initial carbon concentrations
  mmem<-function(X, pars, t){
    
    #mmem model
    deriv<-function(time, state, pars){
      
      with(as.list(c(state, pars)),{
        
        #Carbon uptake
        
        Cu=Vmax*C*Cmic/(Km+C)
        
        #defining stochiometric coefficients
        #PO ratio
        #PO=2.5*(1-((Cu/Cmic)/((Cu/Cmic)+Sover)))
        y=2*1.5*Yatp
        
        #x
        #x=(Cu/Cmic)/((Cu/Cmic)+Sprod)
        
        #coefficients
        psi.r=1/((y*x)+(y*0.64*(1-x))+1)
        psi.g=y*x/((y*x)+(y*0.64*(1-x))+1)
        psi.e=0.64*y*(1-x)/((y*x)+(y*0.64*(1-x))+1)
        
        
        dCmic<--kmic*Cmic+Cu*psi.g
        dC<-kmic*Cmic-Cu+ke*E
        dE<-Cu*psi.e-ke*E
        
        return(list(c(dCmic, dC, dE), r=Cu*psi.r))
        
      })
    }
    
    #function to apply acros 3 different initial substrate concentrations
    base_function<-function(Ci){
      
      
      out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=t, func=deriv)
      return(as.data.frame(out))}
    
    #results
    res<-lapply(X=X, FUN = "base_function") %>% bind_cols()
    
    return(res)
    
  }
  
  #create cost function
  estim<-function(data){
    
    Obs_dat<-data[,c(2, 5, 8, 12, 9)]
    
    m1<-merge(Obs_dat[Obs_dat$Substrate=="Glucose", c(1,3,4, 5)], 
              Obs_dat[Obs_dat$Substrate=="Cellobiose", c(1,3,4, 5)],all = T, by="Time")
    
    m2<-merge(m1, 
              Obs_dat[Obs_dat$Substrate=="Mix", c(1,3,4, 5)],all = T, by="Time")
    
    colnames(m2)<-c("time", "Cmic", "r","E", "Cmic1", "r1", "E1", "Cmic2", "r2", "E2")
    
    
    cost_function<-function(pars){
      
      
      out<-mmem(X=c(33.3*0.75, 35.06*0.75, 22.6*0.75), pars = pars, t=seq(0,130))
      cost<-modCost(model = out, obs = m2, weight = "mean")
      
      return(cost)
      
    }
    
    res<-modMCMC(f=cost_function, p=c(Vmax=0.1, Km=3, kmic=0.01, ke=0.01, x=0.7, Yatp=1),
                 lower=c(Vmax=1e-4, Km=1e-4, kmic=1e-5, ke=1e-6, x=0, Yatp=0.01),
                 upper=c(Vmax=1e4, Km=1e4, kmic=1e5, ke=1e6, x=1, Yatp=30),niter=10000)
    
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
  colnames(parameters.l)<-c("Vmax.l", "Km.l", "kmic.l", "ke.l", "x.l", "Yatp.l")
  
  parameters.u<-rbind(as.data.frame(summary(res[[1]])[7,]), 
                      as.data.frame(summary(res[[2]])[7,]),
                      as.data.frame(summary(res[[3]])[7,]))
  colnames(parameters.u)<-c("Vmax.u", "Km.u", "kmic.u", "ke.u", "x.u", "Yatp.l")
  
  
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
    
    out<-mmem(X=c(33.3*0.75, 35.06*0.75, 22.6*0.75), pars = pars, t=seq(0,130))
    cost<-modCost(model = out, obs = m2)
    
    return(cost)
    
  }
  
  for(i in 1:3){
    
    obs_r<-append(obs_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$obs)
    mod_r<-append(mod_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$mod)
    
  }
  
  OvP_r<-data.frame(obs_r, mod_r)
  
  #logLik calculation
  mu_r<-mean(obs_r)
  variance_r<-sd(obs_r)^2
  
  
  ll_r<--1*sum((obs_r-mod_r)^2)/2/variance_r
  
  
  SSmodel_r<-sum((obs_r-mod_r)^2)
  SSdata_r<-sum((obs_r-mean(obs_r))^2)
  
  rsq_r=1-(SSmodel_r/SSdata_r)
  
  likelihood_r<-c(logLik=ll_r, npar=18, rsq=rsq_r)
  
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
    
    out<-mmem(X=c(33.3*0.75, 35.06*0.7, 22.6*0.75), pars = pars, t=seq(0,130))
    cost<-modCost(model = out, obs = m2)
    
    return(cost)
    
  }
  
  for(i in 1:3){
    obs_Cmic<-append(obs_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$obs)
    mod_Cmic<-append(mod_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$mod)
  }
  
  OvP_Cmic<-data.frame(obs_Cmic, mod_Cmic)
  
  #logLik calculation
  mu_Cmic<-mean(obs_Cmic)
  variance_Cmic<-sd(obs_Cmic)^2
  
  
  ll_Cmic<--1*sum((obs_Cmic-mod_Cmic)^2)/2/variance_Cmic
  
  
  SSmodel_Cmic<-sum((obs_Cmic-mod_Cmic)^2)
  SSdata_Cmic<-sum((obs_Cmic-mean(obs_Cmic))^2)
  
  rsq_Cmic=1-(SSmodel_Cmic/SSdata_Cmic)
  
  likelihood_Cmic<-c(logLik=ll_Cmic, npar=18, rsq=rsq_Cmic)
  
  #OvP for enzymes
  obs_E<-numeric()
  mod_E<-numeric()
  
  cost_E<-function(pars, data){
    
    Obs_dat<-data[,c(2, 5, 9)]
    
    m1<-merge(Obs_dat[Obs_dat$Substrate=="Glucose", c(1,3)], 
              Obs_dat[Obs_dat$Substrate=="Cellobiose", c(1,3)],all = T, by="Time")
    
    m2<-merge(m1, 
              Obs_dat[Obs_dat$Substrate=="Mix", c(1,3)],all = T, by="Time")
    
    colnames(m2)<-c("time", "E", "E1", "E2")
    
    out<-mmem(X=c(33.3*0.75, 35.06*0.7, 22.6*0.75), pars = pars, t=seq(0,130))
    cost<-modCost(model = out, obs = m2)
    
    return(cost)
    
  }
  
  for(i in 1:3){
    obs_E<-append(obs_E, cost_E(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
    mod_E<-append(mod_E, cost_E(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
  }
  
  OvP_E<-data.frame(obs_E, mod_E)
  
  #logLik calculation
  mu_E<-mean(obs_E)
  variance_E<-sd(obs_E)^2
  
  
  ll_E<--1*sum((obs_E-mod_E)^2)/2/variance_E
  
  
  SSmodel_E<-sum((obs_E-mod_E)^2)
  SSdata_E<-sum((obs_E-mean(obs_E))^2)
  
  rsq_E=1-(SSmodel_E/SSdata_E)
  
  likelihood_E<-c(logLik=ll_E, npar=18, rsq=rsq_E)
  
  al<-list()
  
  al$parameters<-parameters.all
  al$OvP_r<-OvP_r
  al$ll_r<-likelihood_r
  al$OvP_Cmic<-OvP_Cmic
  al$ll_Cmic<-likelihood_Cmic
  al$OvP_E<-OvP_E
  al$ll_E<-likelihood_E
  
  return(al)
}

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

mmem_structures_results<-mmem_structures(dat)

stopImplicitCluster()

#showing results
mmem_structures_results$ll_r
mmem_structures_results$ll_Cmic
mmem_structures_results$ll_E
mmem_structures_results$parameters

#storing results
mmem_structures_ll<-mmem_structures_results$ll_r


#Figures
plot_mmem_structures<-mmem_structures_results$OvP_r
ggplot(plot_mmem_structures, aes(obs_r, mod_r))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  ylab(expression(paste("Predicted respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Observed respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  labs(title="Metabolic Microbial - Enzyme model \n for different structures")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

mmem_structures_Cmic<-mmem_structures_results$OvP_Cmic
ggplot(mmem_structures_Cmic, aes(obs_Cmic, mod_Cmic))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,30))+
  scale_y_continuous(limits = c(0,30))+
  ylab(expression(paste("Predicted ", C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Metabolic Microbial - Enzyme model \n for different structures")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

mmem_structures_E<-mmem_structures_results$OvP_E
ggplot(mmem_structures_E, aes(obs_E, mod_E))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,8))+
  scale_y_continuous(limits = c(0,8))+
  ylab(expression(paste("Predicted ", E, " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  E, " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Metabolic Microbial - Enzyme model \n for different structures")+
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

mem_all_ll
mem_structures_ll
mem_substrates_ll
mem_unique_ll

mmem_all_ll
mmem_structures_ll

##############################################################################################
#for different substrates
mmem_substrates<-function(data){
  
  dat<-setDT(data)[, id := .GRP, by = .(Substrate)]
  
  #this is the mmem model tobe fitted across different substrates 
  #with different initial carbon concentrations
  
  
  #mmem model
  deriv<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      #Carbon uptake
      
      Cu=Vmax*C*Cmic/(Km+C)
      
      #defining stochiometric coefficients
      #PO ratio
      #PO=1.5*(1-((Cu/Cmic)/((Cu/Cmic)+Sover)))
      y=2*1.5*Yatp
      
      #x
      #x=(Cu/Cmic)/((Cu/Cmic)+Sprod)
      
      #coefficients
      psi.r=1/((y*x)+(y*0.64*(1-x))+1)
      psi.g=y*x/((y*x)+(y*0.64*(1-x))+1)
      psi.e=0.64*y*(1-x)/((y*x)+(y*0.64*(1-x))+1)
      
      
      dCmic<--kmic*Cmic+Cu*psi.g
      dC<-kmic*Cmic-Cu+ke*E
      dE<-Cu*psi.e-ke*E
      
      return(list(c(dCmic, dC, dE), r=Cu*psi.r))
      
    })
  }
  
  
  
  #create cost function
  estim<-function(data){
    
    Obs_dat<-data[,c(2, 8, 12, 9)]
    
    colnames(Obs_dat)<-c("time", "Cmic", "r","E")
    Ci<-as.numeric(data[1, "DOCinit"])*0.75
    
    cost_function<-function(pars){
      
      
      out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
      cost<-modCost(model = out, obs = Obs_dat, weight = "mean")
      
      return(cost)
      
    }
    
    res<-modMCMC(f=cost_function, p=c(Vmax=0.1, Km=3, kmic=0.01, ke=0.01, x=0.7, Yatp=1),
                 lower=c(Vmax=1e-4, Km=1e-4, kmic=1e-5, ke=1e-6, x=0,  Yatp=0.01),
                 upper=c(Vmax=1e4, Km=1e4, kmic=1e5, ke=1e6, x=1, Yatp=30),niter=10000)
    
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
  colnames(parameters.l)<-c("Vmax.l", "Km.l", "kmic.l", "ke.l", "x.l", "Yatp.l")
  
  parameters.u<-rbind(as.data.frame(summary(res[[1]])[7,]), 
                      as.data.frame(summary(res[[2]])[7,]),
                      as.data.frame(summary(res[[3]])[7,]))
  colnames(parameters.u)<-c("Vmax.u", "Km.u", "kmic.u", "ke.u", "x.u", "Yatp.u")
  
  
  parameters.all<-cbind(parameters, parameters.l, parameters.u)
  
  
  
  
  #OvP for respiration
  obs_r<-numeric()
  mod_r<-numeric()
  
  cost_r<-function(pars, data){
    
    Obs_dat<-data[,c(2, 12)]
    
    colnames(Obs_dat)<-c("time", "r")
    Ci<-as.numeric(data[1, "DOCinit"])*0.75
    
    out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
  }
  
  for(i in 1:3){
    
    obs_r<-append(obs_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$obs)
    mod_r<-append(mod_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$mod)
    
  }
  
  OvP_r<-data.frame(obs_r, mod_r)
  
  #logLik calculation
  mu_r<-mean(obs_r)
  variance_r<-sd(obs_r)^2
  
  
  ll_r<--1*sum((obs_r-mod_r)^2)/2/variance_r
  
  
  SSmodel_r<-sum((obs_r-mod_r)^2)
  SSdata_r<-sum((obs_r-mean(obs_r))^2)
  
  rsq_r=1-(SSmodel_r/SSdata_r)
  
  likelihood_r<-c(logLik=ll_r, npar=18, rsq=rsq_r)
  
  #OvP for biomass
  obs_Cmic<-numeric()
  mod_Cmic<-numeric()
  
  cost_Cmic<-function(pars, data){
    
    Obs_dat<-data[,c(2, 8)]
    
    colnames(Obs_dat)<-c("time", "Cmic")
    Ci<-as.numeric(data[1, "DOCinit"])*0.75
    
    out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
    
    
  }
  
  for(i in 1:3){
    obs_Cmic<-append(obs_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$obs)
    mod_Cmic<-append(mod_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$mod)
  }
  
  OvP_Cmic<-data.frame(obs_Cmic, mod_Cmic)
  
  #logLik calculation
  mu_Cmic<-mean(obs_Cmic)
  variance_Cmic<-sd(obs_Cmic)^2
  
  
  ll_Cmic<--1*sum((obs_Cmic-mod_Cmic)^2)/2/variance_Cmic
  
  
  SSmodel_Cmic<-sum((obs_Cmic-mod_Cmic)^2)
  SSdata_Cmic<-sum((obs_Cmic-mean(obs_Cmic))^2)
  
  rsq_Cmic=1-(SSmodel_Cmic/SSdata_Cmic)
  
  likelihood_Cmic<-c(logLik=ll_Cmic, npar=18, rsq=rsq_Cmic)
  
  #OvP for enzymes
  obs_E<-numeric()
  mod_E<-numeric()
  
  cost_E<-function(pars, data){
    
    Obs_dat<-data[,c(2, 9)]
    
    colnames(Obs_dat)<-c("time", "E")
    Ci<-as.numeric(data[1, "DOCinit"])*0.75
    
    out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
  }
  
  for(i in 1:3){
    obs_E<-append(obs_E, cost_E(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
    mod_E<-append(mod_E, cost_E(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
  }
  
  OvP_E<-data.frame(obs_E, mod_E)
  
  #logLik calculation
  mu_E<-mean(obs_E)
  variance_E<-sd(obs_E)^2
  
  
  ll_E<--1*sum((obs_E-mod_E)^2)/2/variance_E
  
  
  SSmodel_E<-sum((obs_E-mod_E)^2)
  SSdata_E<-sum((obs_E-mean(obs_E))^2)
  
  rsq_E=1-(SSmodel_E/SSdata_E)
  
  likelihood_E<-c(logLik=ll_E, npar=18, rsq=rsq_E)
  
  al<-list()
  
  al$parameters<-parameters.all
  al$OvP_r<-OvP_r
  al$ll_r<-likelihood_r
  al$OvP_Cmic<-OvP_Cmic
  al$ll_Cmic<-likelihood_Cmic
  al$OvP_E<-OvP_E
  al$ll_E<-likelihood_E
  
  return(al)
}

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

mmem_substrates_results<-mmem_substrates(dat)

stopImplicitCluster()

#showing results
mmem_substrates_results$ll_r
mmem_substrates_results$ll_Cmic
mmem_substrates_results$ll_E
mmem_substrates_results$parameters

#storing results
mmem_substrates_ll<-mmem_substrates_results$ll_r


#Figures
plot_mmem_substrates<-mmem_substrates_results$OvP_r
ggplot(plot_mmem_substrates, aes(obs_r, mod_r))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  ylab(expression(paste("Predicted respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Observed respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  labs(title="Metabolic Microbial - Enzyme model \n for different substrates")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

mmem_substrates_Cmic<-mmem_substrates_results$OvP_Cmic
ggplot(mmem_substrates_Cmic, aes(obs_Cmic, mod_Cmic))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,30))+
  scale_y_continuous(limits = c(0,30))+
  ylab(expression(paste("Predicted ", C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Metabolic Microbial - Enzyme model \n for different substrates")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

mmem_substrates_E<-mmem_substrates_results$OvP_E
ggplot(mmem_substrates_E, aes(obs_E, mod_E))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,8))+
  scale_y_continuous(limits = c(0,8))+
  ylab(expression(paste("Predicted ", E, " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  E, " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Metabolic Microbial - Enzyme model \n for different substrates")+
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

mem_all_ll
mem_structures_ll
mem_substrates_ll
mem_unique_ll

mmem_all_ll
mmem_structures_ll
mmem_substrates_ll


#############################################################################################
#for all combinations of substrates and structures
mmem_unique<-function(data){
  
  dat<-setDT(data)[, id := .GRP, by = .(Substrate, Structure)]
  
  
  #mmem model
  deriv<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      #Carbon uptake
      
      Cu=Vmax*C*Cmic/(Km+C)
      
      #defining stochiometric coefficients
      #PO ratio
      #PO=2.5*(1-((Cu/Cmic)/((Cu/Cmic)+Sover)))
      y=2*1.5*Yatp
      
      #x
      #x=(Cu/Cmic)/((Cu/Cmic)+Sprod)
      
      #coefficients
      psi.r=1/((y*x)+(y*0.64*(1-x))+1)
      psi.g=y*x/((y*x)+(y*0.64*(1-x))+1)
      psi.e=0.64*y*(1-x)/((y*x)+(y*0.64*(1-x))+1)
      
      
      dCmic<--kmic*Cmic+Cu*psi.g
      dC<-kmic*Cmic-Cu+ke*E
      dE<-Cu*psi.e-ke*E
      
      return(list(c(dCmic, dC, dE), r=Cu*psi.r))
      
    })
  }
  
  
  #create cost function
  estim<-function(data){
    
    Obs_dat<-data[,c(2, 8, 12, 9)]
    
    colnames(Obs_dat)<-c("time", "Cmic", "r","E")
    Ci<-as.numeric(data[1, "DOCinit"])*0.75
    
    cost_function<-function(pars){
      
      
      out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
      cost<-modCost(model = out, obs = Obs_dat, weight = "mean")
      
      return(cost)
      
    }
    
    res<-modMCMC(f=cost_function, p=c(Vmax=0.1, Km=3, kmic=0.01, ke=0.01, x=0.7, Yatp=1),
                 lower=c(Vmax=1e-4, Km=1e-4, kmic=1e-5, ke=1e-6, x=0, Yatp=0.01),
                 upper=c(Vmax=1e4, Km=1e4, kmic=1e5, ke=1e6, x=1, Yatp=30),niter=10000)
    
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
  parameters$Structure<-c(res[[1]]$Structure,
                          res[[2]]$Structure,
                          res[[3]]$Structure,
                          res[[4]]$Structure,
                          res[[5]]$Structure,
                          res[[6]]$Structure,
                          res[[7]]$Structure,
                          res[[8]]$Structure,
                          res[[9]]$Structure)
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
  colnames(parameters.l)<-c("Vmax.l", "Km.l", "kmic.l", "ke.l", "x.l","Yatp.l")
  
  parameters.u<-rbind(as.data.frame(summary(res[[1]])[7,]), 
                      as.data.frame(summary(res[[2]])[7,]),
                      as.data.frame(summary(res[[3]])[7,]), 
                      as.data.frame(summary(res[[4]])[7,]), 
                      as.data.frame(summary(res[[5]])[7,]),
                      as.data.frame(summary(res[[6]])[7,]),
                      as.data.frame(summary(res[[7]])[7,]), 
                      as.data.frame(summary(res[[8]])[7,]),
                      as.data.frame(summary(res[[9]])[7,]))
  colnames(parameters.u)<-c("Vmax.u", "Km.u", "kmic.u", "ke.u", "x.u", "Yatp.l")
  
  
  parameters.all<-cbind(parameters, parameters.l, parameters.u)
  
  
  
  
  #OvP for respiration
  obs_r<-numeric()
  mod_r<-numeric()
  
  cost_r<-function(pars, data){
    
    Obs_dat<-data[,c(2, 12)]
    
    colnames(Obs_dat)<-c("time", "r")
    Ci<-as.numeric(data[1, "DOCinit"])*0.75
    
    out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
  }
  
  for(i in 1:9){
    
    obs_r<-append(obs_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$obs)
    mod_r<-append(mod_r, cost_r(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$mod)
    
  }
  
  OvP_r<-data.frame(obs_r, mod_r)
  
  #logLik calculation
  mu_r<-mean(obs_r)
  variance_r<-sd(obs_r)^2
  
  
  ll_r<--1*sum((obs_r-mod_r)^2)/2/variance_r
  
  
  SSmodel_r<-sum((obs_r-mod_r)^2)
  SSdata_r<-sum((obs_r-mean(obs_r))^2)
  
  rsq_r=1-(SSmodel_r/SSdata_r)
  
  likelihood_r<-c(logLik=ll_r, npar=54, rsq=rsq_r)
  
  #OvP for biomass
  obs_Cmic<-numeric()
  mod_Cmic<-numeric()
  
  cost_Cmic<-function(pars, data){
    
    Obs_dat<-data[,c(2, 8)]
    
    colnames(Obs_dat)<-c("time", "Cmic")
    Ci<-as.numeric(data[1, "DOCinit"])*0.75
    
    out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
    
    
  }
  
  for(i in 1:9){
    obs_Cmic<-append(obs_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$obs)
    mod_Cmic<-append(mod_Cmic, cost_Cmic(pars=summary(res[[i]])[6,], data=dat[dat$id==i,])$residuals$mod)
  }
  
  OvP_Cmic<-data.frame(obs_Cmic, mod_Cmic)
  
  #logLik calculation
  mu_Cmic<-mean(obs_Cmic)
  variance_Cmic<-sd(obs_Cmic)^2
  
  
  ll_Cmic<--1*sum((obs_Cmic-mod_Cmic)^2)/2/variance_Cmic
  
  
  SSmodel_Cmic<-sum((obs_Cmic-mod_Cmic)^2)
  SSdata_Cmic<-sum((obs_Cmic-mean(obs_Cmic))^2)
  
  rsq_Cmic=1-(SSmodel_Cmic/SSdata_Cmic)
  
  likelihood_Cmic<-c(logLik=ll_Cmic, npar=54, rsq=rsq_Cmic)
  
  #OvP for enzymes
  obs_E<-numeric()
  mod_E<-numeric()
  
  cost_E<-function(pars, data){
    
    Obs_dat<-data[,c(2, 9)]
    
    colnames(Obs_dat)<-c("time", "E")
    Ci<-as.numeric(data[1, "DOCinit"])*0.75
    
    out<-ode(y=c(Cmic=6.98, C=Ci, E=0), parms=pars, times=seq(0,130), func=deriv)
    cost<-modCost(model = out, obs = Obs_dat)
    
    return(cost)
    
  }
  
  for(i in 1:9){
    obs_E<-append(obs_E, cost_E(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$obs)
    mod_E<-append(mod_E, cost_E(pars=summary(res[[i]])[6,], data=dat[dat$id==i, ])$residuals$mod)
  }
  
  OvP_E<-data.frame(obs_E, mod_E)
  
  #logLik calculation
  mu_E<-mean(obs_E)
  variance_E<-sd(obs_E)^2
  
  
  ll_E<--1*sum((obs_E-mod_E)^2)/2/variance_E
  
  
  SSmodel_E<-sum((obs_E-mod_E)^2)
  SSdata_E<-sum((obs_E-mean(obs_E))^2)
  
  rsq_E=1-(SSmodel_E/SSdata_E)
  
  likelihood_E<-c(logLik=ll_E, npar=54, rsq=rsq_E)
  
  al<-list()
  
  al$parameters<-parameters.all
  al$OvP_r<-OvP_r
  al$ll_r<-likelihood_r
  al$OvP_Cmic<-OvP_Cmic
  al$ll_Cmic<-likelihood_Cmic
  al$OvP_E<-OvP_E
  al$ll_E<-likelihood_E
  
  return(al)
}

no_cors<-detectCores()-1
cl<-makeCluster(no_cors)
registerDoParallel(cl)

mmem_unique_results<-mmem_unique(dat)

stopImplicitCluster()

#showing results
mmem_unique_results$ll_r
mmem_unique_results$ll_Cmic
mmem_unique_results$ll_E
mmem_unique_results$parameters

#storing results
mmem_unique_ll<-mmem_unique_results$ll_r


#Figures
plot_mmem_unique<-mmem_unique_results$OvP_r
ggplot(plot_mmem_unique, aes(obs_r, mod_r))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  ylab(expression(paste("Predicted respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  xlab(expression(paste("Observed respiration rate (", mu, "mol ", ml^{-1}~h^{-1}, ")")))+
  labs(title="Metabolic Microbial - Enzyme model \n for all combinations of different substrates and structures")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

mmem_unique_Cmic<-mmem_unique_results$OvP_Cmic
ggplot(mmem_unique_Cmic, aes(obs_Cmic, mod_Cmic))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,30))+
  scale_y_continuous(limits = c(0,30))+
  ylab(expression(paste("Predicted ", C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  C[MIC], " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Metabolic Microbial - Enzyme model \n for all combinations of different substrates and structures")+
  theme(plot.title=element_text(size=14, face="bold.italic", hjust=0.5))

mmem_unique_E<-mmem_unique_results$OvP_E
ggplot(mmem_unique_E, aes(obs_E, mod_E))+theme_min+geom_point(cex=6, pch=1)+
  geom_abline(intercept = 0, slope=1, lwd=1.2)+
  scale_x_continuous(limits = c(0,8))+
  scale_y_continuous(limits = c(0,8))+
  ylab(expression(paste("Predicted ", E, " (", mu, "mol ", ml^{-1}, ")")))+
  xlab(expression(paste("Observed ",  E, " (", mu, "mol ", ml^{-1}, ")")))+
  labs(title="Metabolic Microbial - Enzyme model \n for all combinations of different substrates and structures")+
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

mem_all_ll
mem_structures_ll
mem_substrates_ll
mem_unique_ll

mmem_all_ll
mmem_structures_ll
mmem_substrates_ll
mmem_unique_ll

#######AIC calculation 2*npar - 2*loglik
(decay_all_ll[4]<-2*decay_all_ll[2]-2*decay_all_ll[1])
(decay_substrates_ll[4]<-c(AIC=2*decay_substrates_ll[2]-2*decay_substrates_ll[1]))
(decay_structures_ll[4]<-2*decay_structures_ll[2]-2*decay_structures_ll[1])
(decay_unique_ll[4]<-2*decay_unique_ll[2]-2*decay_unique_ll[1])

(monod_all_ll[4]<-2*monod_all_ll[2]-2*monod_all_ll[1])
(monod_structures_ll[4]<-2*monod_structures_ll[2]-2*monod_structures_ll[1])
(monod_substrates_ll[4]<-2*monod_substrates_ll[2]-2*monod_substrates_ll[1])
(monod_unique_ll[4]<-2*monod_unique_ll[2]-2*monod_unique_ll[1])

(mem_all_ll[4]<-2*mem_all_ll[2]-2*mem_all_ll[1])
(mem_structures_ll[4]<-2*mem_structures_ll[2]-2*mem_structures_ll[1])
(mem_substrates_ll[4]<-2*mem_substrates_ll[2]-2*mem_substrates_ll[1])
(mem_unique_ll[4]<-2*mem_unique_ll[2]-2*mem_unique_ll[1])

(mmem_all_ll[4]<-2*mmem_all_ll[2]-2*mmem_all_ll[1])
(mmem_structures_ll[4]<-2*mmem_structures_ll[2]-2*mmem_structures_ll[1])
(mmem_substrates_ll[4]<-2*mmem_substrates_ll[2]-2*mmem_substrates_ll[1])
(mmem_unique_ll[4]<-2*mmem_unique_ll[2]-2*mmem_unique_ll[1])


#comparing model fits based on AIC 
#monod growth for different substrates vs mem for different substrates 
exp((monod_substrates_ll[4]-mem_substrates_ll[4])/2)

#monod growth for different substrates vs mmem for different substrates 
exp((monod_substrates_ll[4]-mmem_substrates_ll[4])/2)

#Biomass AIC
(AICmonod_all_Cmic<-2*monod_all_results$ll_Cmic[2]-2*monod_all_results$ll_Cmic[1])
(AICmonod_structures_Cmic<-2*monod_structures_results$ll_Cmic[2]-2*monod_structures_results$ll_Cmic[1])
(AICmonod_substrates_Cmic<-2*monod_substrates_results$ll_Cmic[2]-2*monod_substrates_results$ll_Cmic[1])
(AICmonod_unique_Cmic<-2*monod_unique_results$ll_Cmic[2]-2*monod_unique_results$ll_Cmic[1])

(AICmem_all_Cmic<-2*mem_all_results$ll_Cmic[2]-2*mem_all_results$ll_Cmic[1])
(AICmem_structures_Cmic<-2*mem_structures_results$ll_Cmic[2]-2*mem_structures_results$ll_Cmic[1])
(AICmem_substrates_Cmic<-2*mem_substrates_results$ll_Cmic[2]-2*mem_substrates_results$ll_Cmic[1])
(AICmem_unique_Cmic<-2*mem_unique_results$ll_Cmic[2]-2*mem_unique_results$ll_Cmic[1])

(AICmmem_all_Cmic<-2*mmem_all_results$ll_Cmic[2]-2*mmem_all_results$ll_Cmic[1])
(AICmmem_structures_Cmic<-2*mmem_structures_results$ll_Cmic[2]-2*mmem_structures_results$ll_Cmic[1])
(AICmmem_substrates_Cmic<-2*mmem_substrates_results$ll_Cmic[2]-2*mmem_substrates_results$ll_Cmic[1])
(AICmmem_unique_Cmic<-2*mmem_unique_results$ll_Cmic[2]-2*mmem_unique_results$ll_Cmic[1])


#monod growth for different substrates vs mem for different substrates 
exp((AICmonod_substrates_Cmic-AICmem_substrates_Cmic)/2)

#monod growth for different substrates vs mmem for different substrates 
exp((AICmonod_substrates_Cmic-AICmmem_substrates_Cmic)/2)


#Enzyme AIC
(AICmem_all_E<-2*mem_all_results$ll_E[2]-2*mem_all_results$ll_E[1])
(AICmem_structures_E<-2*mem_structures_results$ll_E[2]-2*mem_structures_results$ll_E[1])
(AICmem_substrates_E<-2*mem_substrates_results$ll_E[2]-2*mem_substrates_results$ll_E[1])
(AICmem_unique_E<-2*mem_unique_results$ll_E[2]-2*mem_unique_results$ll_E[1])

(AICmmem_all_E<-2*mmem_all_results$ll_E[2]-2*mmem_all_results$ll_E[1])
(AICmmem_structures_E<-2*mmem_structures_results$ll_E[2]-2*mmem_structures_results$ll_E[1])
(AICmmem_substrates_E<-2*mmem_substrates_results$ll_E[2]-2*mmem_substrates_results$ll_E[1])
(AICmmem_unique_E<-2*mmem_unique_results$ll_E[2]-2*mmem_unique_results$ll_E[1])

#mem for different substrates vs mmem for different substrates 
exp((AICmmem_substrates_E-AICmem_substrates_Cmic)/2)
