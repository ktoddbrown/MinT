Mr_lin_est<-function(odeset, par_const){
  #Cost function
  cost<-function(x){
    p<-x
    #p<-out_const$pars
    names(p)<-c("Ecat", "f", "chi", "Yue", "j", "k", "fpr", "fps", "fds", "RSinit", "I", "B")
    
    #Initial conditions
    Sinit<-mean(odeset$DNA.initc, na.rm=T)/p[["fds"]]
    #Sinit<-mean(m0CB$DNA.initc, na.rm=T)/p[["fds"]]
    Rinit<-p[["RSinit"]]*Sinit
    
    #Simulation
    yhat_all <- as.data.frame(ode(y=c(Cs=25, DOC=0, R=Rinit, S=Sinit, Enz=0), func = Mr_lin, times = seq(0,120),
                                  parms = p[-which(names(p)=="RSinit")]))
    
    #Filter the output
    yhat <- yhat_all[, c("Protinc", "DNAc", "Enz", "r", "time")]
    colnames(yhat)<-c("Protinc", "DNAc", "Protoutc", "r", "time")
    Yhat <- melt(yhat, id.vars=c("time"))
    colnames(Yhat)<-c("time", "variable", "yhat")
    
    #Observations
    obs <- odeset[, c("Protinc", "DNAc", "Protoutc", "r", "Time2")]
    Obs <- melt(obs, id.vars=c("Time2"))
    colnames(Obs)<-c("time", "variable", "obs")
    
    #Merge both
    Yhat <- merge(Yhat, Obs, by = c("time", "variable"))
    
    #Add weighting factor
    yweights<-Yhat %>% group_by(variable) %>% summarise(weights=mean(obs, na.rm = T))
    
    #match with the Yhat data frame
    Yhat<-merge(Yhat, yweights, by=c("variable"))
    
    #now, the root mean square error is calculated
    NRMSE<-as.numeric(Yhat %>% group_by(variable) %>% summarise(NRMSE=sum((((yhat-obs)/mean(weights))^2), na.rm = T)) %>%
                        summarise(NRMSE=sum(NRMSE)))
    
    return(NRMSE)
  }
  
  #Goodness of fit function
  good<-function(x){
    p<-x
    names(p)<-c("Ecat", "f", "chi", "Yue", "j", "k", "fpr", "fps", "fds", "RSinit", "I", "B")
    
    #Initial conditions
    Sinit<-mean(odeset$DNA.initc, na.rm=T)/p[["fds"]]
    Rinit<-p[["RSinit"]]*Sinit
    
    #Simulation
    yhat_all <- as.data.frame(ode(y=c(Cs=25, DOC=0, R=Rinit, S=Sinit, Enz=0), func = Mr_lin, times = seq(0,120),
                                  parms = p[-which(names(p)=="RSinit")]))
    
    #Filter the output
    yhat <- yhat_all[, c("Protinc", "DNAc", "Enz", "r", "time")]
    colnames(yhat)<-c("Protinc", "DNAc", "Protoutc", "r", "time")
    Yhat <- melt(yhat, id.vars=c("time"))
    colnames(Yhat)<-c("time", "variable", "yhat")
    
    #Observations
    obs <- odeset[, c("Protinc", "DNAc", "Protoutc", "r", "Time2")]
    Obs <- melt(obs, id.vars=c("Time2"))
    colnames(Obs)<-c("time", "variable", "obs")
    
    #Merge both
    Yhat <- merge(Yhat, Obs, by = c("time", "variable"))
    
    #rsquared calculation for each variable
    Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-yhat)^2), na.rm = T),
                                                    SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                    ll=-sum(((obs-yhat)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
    Gfit$R2<-with(Gfit, 1-SSres/SStot)
    Gfit$N<-length(p)
    Gfit$AIC<-with(Gfit, 2*N-2*ll)
    
    rsq_out<-list(Yhat=Yhat, Gfit=Gfit)
    
    return(rsq_out)
  }
  
  #Estimate the parameters
  #approximate parameter estimation is done by MCMC method
  par_c<-par_const[-which(names(par_const)=="Mr")]
  par_t<-par_const[["Mr"]]
  par_mcmc<-modMCMC(f=cost, p=c(par_c, I=par_t, B=1),
                    lower=c(par_c*0.5, I=par_t*0.5, B=-100),
                    upper=c(Ecat=10, f=3e-1, chi=1, Yue=0.9, j=1, 
                            k=1e-1, fpr=1, fps=1, fds=1, RSinit=50, 
                            I=par_t*2, B=100), niter=10000)
  
  #lower and upper limits for parameters are extracted
  pl<-summary(par_mcmc)["min",]
  pu<-summary(par_mcmc)["max",]
  
  #these limits are used to find global optimum by DEoptim
  opt_par<-DEoptim(fn=cost, lower=pl, upper=pu,
                   control = c(itermax = 10000, steptol = 50, reltol = 1e-8,
                               trace=FALSE, strategy=3, NP=250))
  
  #goodness of fit
  fit<-good(opt_par$optim$bestmem)
  
  #best parameters
  p<-opt_par$optim$bestmem
  names(p)<-parnames
  
  #return list with opt_par and par_prof
  estim_out<-list(pars=p, par_mcmc=par_mcmc, fit=fit)
  
  return(estim_out)
}